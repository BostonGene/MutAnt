import csv
from typing import Literal

import pandas as pd
import pysam
import ray
from loguru import logger
from tqdm import tqdm

from dbnsfp_annotation.utils import (
    DATABASE_PATH,
    GENE_BASE_PATH,
    annotate_chunk,
    convert_chromosome_id,
    init_ray,
    liftover_hg19_hg38,
)


def annotate_maf_dbnsfp(
    maf: [str, pd.DataFrame],
    genome_version: str = "hg38",
    database_path: str = DATABASE_PATH,
    gene_base_path: str = GENE_BASE_PATH,
    n_jobs: int = 6,
    host: str = "0.0.0.0",
) -> pd.DataFrame:
    """
    Annotates MAF-file SNPs using dbNSFP. DB reference version must be compatible with MAF reference version.
    :param database_path: path to dbNSFP joined tabix indexed and GZ-compressed database file, string
    :gene_base_path: path to dbNSFP gene_complete annotation file
    :param maf: path to tab-separated MAF-file, string
    :param genome_version: 'hg19' or 'hg38' version of genome
    :param n_jobs: number of concurrent workers for annotation
    :param host: ip adress of machine to create multi-processing. default - current machine
    :return annotated pd.DataFrame
    """

    if isinstance(maf, str):
        sample_maf = pd.read_csv(maf, sep="\t", comment="#")

    if isinstance(maf, pd.DataFrame):
        sample_maf = maf

    if genome_version.lower() not in ["hg19", "hg38"]:
        logger.warning(
            f"Supported genome versions are hg19 or hg38, you passed {genome_version}. Setting hg38 as default"
        )
        genome_version = "hg38"

    if genome_version.lower() == "hg19":
        sample_maf = liftover_hg19_hg38(sample_maf, file_type="maf")

    logger.info(
        "The quality of annotation improves if the following columns are present in maf file:\
                RS_ID (rs_dbsnp), ensembl transcript ID, HGVSp (any version) and HGVSc."
    )

    # Valid names.
    rsid_columns = ["rs_dbsnp", "dbsnp_rs", "rs_id", "rsid"]
    enstid_columns = ["transcript_id", "enst_id", "enstid"]

    # Prepare specific columns names.
    maf_rsid_col = any([i for i in sample_maf.columns if i.lower() in rsid_columns])
    if maf_rsid_col:
        maf_rsid_col = [i for i in sample_maf.columns if i.lower() in rsid_columns][0]

    enstid_columns = ["transcript_id", "enst_id", "enstid"]

    maf_enstid_col = any([i for i in sample_maf.columns if i.lower() in enstid_columns])
    if maf_enstid_col:
        maf_enstid_col = [i for i in sample_maf.columns if i.lower() in enstid_columns][0]

    hgvsp_in_maf = [i for i in sample_maf.columns if "hgvsp" in i.lower()]

    hgvsc_in_maf = [i for i in sample_maf.columns if "hgvsc" in i.lower()]

    # Preloading indexed file
    init_ray(host=host, n_jobs=n_jobs)
    try:
        dbNSFP_file = pysam.TabixFile(database_path, "r", encoding="UTF-8", parser=pysam.asTuple())

        dbNSFP_columns = dbNSFP_file.header[0].split()
        sample_maf_annotated = []

        sample_maf_annotated = ray.get(
            [
                annotate_chunk.remote(
                    chromosome_id,
                    chromosome_data,
                    database_path,
                    hgvsp_in_maf,
                    hgvsc_in_maf,
                    dbNSFP_columns,
                    maf_rsid_col,
                    maf_enstid_col,
                )
                for chromosome_id, chromosome_data in sample_maf.groupby("Chromosome")
            ]
        )

    except Exception as e:
        print(e)
        ray.shutdown()
    finally:
        ray.shutdown()

    sample_maf_annotated = pd.concat(sample_maf_annotated)
    sample_maf_annotated = sample_maf_annotated.drop(["#chr", "pos(1-based)", "ref", "alt"], axis=1)

    sample_maf_annotated = sample_maf_annotated.reset_index(drop=True)
    target_col = "Hugo_Symbol" if "Hugo_Symbol" in sample_maf_annotated else "genename"

    # 3 sec loading
    genes_info_df = pd.read_csv(gene_base_path, sep="\t")

    sample_maf_annotated[f"solo_{target_col}"] = sample_maf_annotated[target_col].map(lambda x: str(x).split(";")[0])

    sample_maf_annotated = sample_maf_annotated.merge(
        genes_info_df, how="left", left_on=f"solo_{target_col}", right_on="Gene_name"
    )

    return sample_maf_annotated


def annotate_vcf_dbnsfp(
    vcf: [str, pd.DataFrame],
    genome_version: Literal["hg38", "hg19"],
    database_path: str = DATABASE_PATH,
    gene_base_path: str = GENE_BASE_PATH,
    n_jobs: int = 6,
) -> pd.DataFrame:
    """
    Annotates VCF-file SNPs using dbNSFP. DB reference version must be compatible with VCF reference version.
    :param database_path: path to dbNSFP joined tabix indexed and GZ-compressed database file, string
    :gene_base_path: path to dbNSFP gene_complete annotation file
    :param vcf: path to tab-separated VCF-file, string
    :param genome_version: 'hg19' or 'hg38' version of genome
    :param n_jobs: number of concurrent workers for annotation
    :return annotated pd.DataFrame
    """

    if isinstance(vcf, str):
        with open(vcf) as vcf_file:
            vcf_reader = csv.reader(vcf_file, delimiter="\t")

            vcf_filtered = (i for i in vcf_reader if not i[0].startswith("##"))

            df_columns = next(vcf_filtered)

            sample_vcf = pd.DataFrame(vcf_filtered, columns=df_columns)

    else:
        sample_vcf = vcf

    sample_vcf.columns = [
        {
            "#chr": "#CHROM",
            "#CHR": "#CHROM",
            "pos(1-based)": "POS",
            "ref": "REF",
            "alt": "ALT",
        }.get(i, i)
        for i in sample_vcf.columns
    ]

    if genome_version.lower() not in ["hg19", "hg38"]:
        raise ValueError("Unsupported genome version (hg19 or h38 are supported)")

    if genome_version.lower() == "hg19":
        sample_vcf = liftover_hg19_hg38(sample_vcf, file_type="vcf")

    sample_vcf["POS"] = sample_vcf["POS"].astype(int)

    # Preloading indexed file
    with pysam.TabixFile(database_path, parser=pysam.asTuple(), encoding="UTF-8", threads=n_jobs) as dbNSFP_file:

        dbNSFP_columns = dbNSFP_file.header[0].split()

        sample_vcf_annotated = []

        for chromosome_id, chromosome_data in tqdm(sample_vcf.groupby("#CHROM")):
            try:
                chromosome_id = convert_chromosome_id(chromosome_id)
            except AttributeError:
                pass

            intersected_variants = []

            for _, variant in chromosome_data.iterrows():

                # Converting SNP coordinate from 1 based to pysam-compatible 0-based
                # TODO: fix multiple AA subs for different transcripts for VCF
                try:
                    for tabix_fetched_row in dbNSFP_file.fetch(chromosome_id, variant["POS"] - 1, variant["POS"]):
                        intersected_variants.append(tabix_fetched_row)
                except (ValueError, KeyError) as e:
                    logger.info(e)
                    continue

            intersected_variants = pd.DataFrame(intersected_variants, columns=dbNSFP_columns)

            intersected_variants["pos(1-based)"] = intersected_variants["pos(1-based)"].apply(int)

            intersected_variants = pd.merge(
                chromosome_data,
                intersected_variants,
                left_on=["POS", "REF", "ALT"],
                right_on=["pos(1-based)", "ref", "alt"],
                how="left",
                suffixes=(None, "_from_dbNSFP"),
            )

            sample_vcf_annotated.append(intersected_variants)

    sample_vcf_annotated = pd.concat(sample_vcf_annotated)
    sample_vcf_annotated = sample_vcf_annotated.drop(["pos(1-based)", "ref", "alt"], axis=1)

    sample_vcf_annotated = sample_vcf_annotated.reset_index(drop=True)

    target_col = "Hugo_Symbol" if "Hugo_Symbol" in sample_vcf_annotated else "genename"

    genes_info_df = pd.read_csv(gene_base_path, sep="\t")

    sample_vcf_annotated[f"solo_{target_col}"] = sample_vcf_annotated[target_col].map(lambda x: str(x).split(";")[0])

    sample_vcf_annotated = sample_vcf_annotated.merge(
        genes_info_df, how="left", left_on=f"solo_{target_col}", right_on="Gene_name"
    )

    return sample_vcf_annotated
