import warnings
from typing import Literal, Optional, Tuple

import numpy as np
import pandas as pd
import pysam
import ray
from liftover import ChainFile, get_lifter

warnings.filterwarnings("ignore")
DATABASE_PATH = "/uftp/Databases/dbNSFP/dbNSFP4.5a/dbNSFPv4.5a_custombuild.gz"
GENE_BASE_PATH = "/uftp/Databases/dbNSFP/dbNSFP4.5a/dbNSFP4.5_gene.complete.gz"


def liftover_hg19_hg38(
    maf: pd.DataFrame,
    file_type: Literal["maf", "vcf"],
    chain_file: str = None,
):
    """
    Liftover hg19 genome position to hg38 position

    :param maf: input VCF or MAF file
    :param file_type: 'maf' or 'vcf'
    :param chain_file: path to chain file if necessary
    :return: Lifted over VCF/MAF file
    """

    if chain_file:
        converter = ChainFile(chain_file, "hg19", "hg38")
    else:
        converter = get_lifter("hg19", "hg38")
    if file_type.lower() == "maf":
        chr_arr, start_pos_arr, end_pos_arr = [], [], []
        for chr_, start_pos, end_pos in zip(
            maf["Chromosome"],
            maf["Start_Position"].astype(int),
            maf["End_Position"].astype(int),
        ):
            try:
                chrom, start_pos, strand = converter.query(chr_, start_pos)[0]
                chrom, end_pos, strand = converter.query(chr_, end_pos)[0]
                chr_arr.append(chrom)
                start_pos_arr.append(start_pos)
                end_pos_arr.append(end_pos)
            except (KeyError, IndexError):
                chr_arr.append(chr_)
                start_pos_arr.append(start_pos)
                end_pos_arr.append(end_pos)
        maf["Chromosome"] = chr_arr
        maf["Start_Position"] = start_pos_arr
        maf["End_Position"] = end_pos_arr
        maf["NCBI_Build"] = "GRCh38"
    elif file_type.lower() == "vcf":
        chr_arr, pos_arr = [], []
        for chr_, pos in zip(maf["#CHROM"], maf["POS"].astype(int)):
            try:
                chr_, pos, strand = converter.query(chr_, pos)[0]
                chr_arr.append(chr_)
                pos_arr.append(pos)
            except (KeyError, IndexError):
                chr_arr.append(chr_)
                pos_arr.append(pos)
        maf["#CHROM"] = chr_arr
        maf["POS"] = pos_arr
    return maf


def convert_chromosome_id(chrom_id: str):
    """
    Converts chromosome id string to format compatible with dbNSFP
    :param chrom_id: chromosome id string
    :return: converted chromosome id string
    """
    converted_chromosome_id = chrom_id.replace("chr", "")
    if converted_chromosome_id.isdigit():
        return int(converted_chromosome_id)
    return converted_chromosome_id


def fix_chromosome_prefix(x):
    """
    fix chromosome prefix for different types of maf/vcf

    :param x: string with chromosome name/number
    :return: fixed chromosome name
    """
    if str(x).startswith("chrMT"):
        return "chrM"

    if str(x).startswith("chr"):
        return x

    if x == "MT":
        x = "M"

    return "chr" + str(x)


def match_rsid(
    variant: pd.Series, tabix_fetched_row: list, maf_rsid_col: Optional[list], dbNSFP_columns: list
) -> Tuple[list, bool]:
    """
    Filter out rows from dbNSFP database fetch by rsid.

    Returns: filtered list and bool, representing whether even one element were filtered out.
    """
    matched = []
    maf_rsid = (
        variant[maf_rsid_col].tolist() if isinstance(variant[maf_rsid_col], pd.Series) else [variant[maf_rsid_col]]
    )
    for row in tabix_fetched_row:
        if row[dbNSFP_columns.index("rs_dbSNP")] in maf_rsid:
            matched.append(row)
    if not matched:
        return tabix_fetched_row, True
    return matched, False


def match_hgvsp(
    variant: pd.Series, tabix_fetched_row: list, hgvsp_in_maf: Optional[list], dbNSFP_columns: list
) -> Tuple[list, bool]:
    """
    Filter out rows from dbNSFP database fetch by HGVSp.

    Returns: filtered list and bool, representing whether even one element were filtered out.
    """
    matched = []
    maf_hgvsp = (
        variant[hgvsp_in_maf].tolist() if isinstance(variant[hgvsp_in_maf], pd.Series) else [variant[hgvsp_in_maf]]
    )
    for row in tabix_fetched_row:
        if any(
            item in maf_hgvsp
            for item in [
                elem
                for item in [row[i] for i in [num for (num, col) in enumerate(dbNSFP_columns) if "hgvsp" in col.lower()]]
                for elem in item.split(";")
            ]
        ):
            matched.append(row)
    if not matched:
        return tabix_fetched_row, True
    return matched, False


def match_hgvsc(
    variant: pd.Series, tabix_fetched_row: list, hgvsc_in_maf: Optional[list], dbNSFP_columns: list
) -> Tuple[list, bool]:
    """
    Filter out rows from dbNSFP database fetch by HGVSc.

    Returns: filtered list and bool, representing whether even one element were filtered out.
    """
    matched = []
    maf_hgvsc = (
        variant[hgvsc_in_maf].tolist() if isinstance(variant[hgvsc_in_maf], pd.Series) else [variant[hgvsc_in_maf]]
    )
    for row in tabix_fetched_row:
        if any(
            item in maf_hgvsc
            for item in [
                elem
                for item in [row[i] for i in [num for (num, col) in enumerate(dbNSFP_columns) if "hgvsc" in col.lower()]]
                for elem in item.split(";")
            ]
        ):
            matched.append(row)
    if not matched:
        return tabix_fetched_row, True
    return matched, False


def match_enstid(
    variant: pd.Series, tabix_fetched_row: list, maf_enstid_col: list, dbNSFP_columns: list
) -> Tuple[list, bool]:
    """
    Filter out rows from dbNSFP database fetch by ENSTID.

    Returns: filtered list and bool, representing whether even one element were filtered out.
    """
    matched = []
    maf_enstid = (
        variant[maf_enstid_col].tolist() if isinstance(variant[maf_enstid_col], pd.Series) else [variant[maf_enstid_col]]
    )
    for row in tabix_fetched_row:
        if any(
            item in maf_enstid
            for item in [
                elem
                for item in [
                    row[i]
                    for i in [num for (num, col) in enumerate(dbNSFP_columns) if "ensembl_transcriptid" in col.lower()]
                ]
                for elem in item.split(";")
            ]
        ):
            matched.append(row)
    if not matched:
        return tabix_fetched_row, True
    return matched, False


def match_ref_alt(variant: pd.Series, tabix_fetched_row: list) -> Tuple[list, bool]:
    """
    Filter out rows from dbNSFP database fetch by Reference and Alternative alleles.

    Returns: filtered list and bool, representing whether even one element were filtered out.
    """
    matched = []
    for row in tabix_fetched_row:
        if variant["Reference_Allele"] in row[2] and variant["Tumor_Seq_Allele2"] in row[3]:
            matched.append(row)
    if not matched:
        return tabix_fetched_row, True
    return matched, False


@ray.remote(num_cpus=1)
def annotate_chunk(
    chromosome_id: str,
    chromosome_data: pd.DataFrame,
    database_path: str,
    hgvsp_in_maf: list,
    hgvsc_in_maf: list,
    dbNSFP_columns: list,
    maf_rsid_col: Optional[list] = None,
    maf_enstid_col: Optional[list] = None,
) -> pd.DataFrame:
    """
    Annotate chromosome specific variants.

    Args:
        chromosome_id: number of chromosome (X and Y and MT kept untouched).
        chromosome_data: pd.DataFrame with mutations from chromosome = chromosome_id.
        database_path: path to dbNSFP database.
        hgvsp_in_maf: list of column names with hgvsp names.
        hgvsc_in_maf: list of column names with hgvsc names.
        dbNSFP_columns: list with names of columns from dbNSFP database file.
        maf_rsid_col: list of column names with rsid names.
        maf_enstid_col: list of column names with enstid names.

    Returns:
        annotated pd.DataFrame by dbNSFP database
    """
    chromosome_id = convert_chromosome_id(chromosome_id)

    intersected_variants = []

    # Filtering out non SNP events
    snp_data = chromosome_data[chromosome_data["Variant_Type"] == "SNP"]

    dbNSFP_file = pysam.TabixFile(database_path, "r", encoding="UTF-8", parser=pysam.asTuple())

    for _, variant in snp_data.iterrows():
        fetched_row = []
        # Converting SNP coordinate from 1 based to pysam-compatible 0-based
        try:
            fetched_row = list(
                dbNSFP_file.fetch(
                    chromosome_id, variant["Start_Position"] - 1, variant["End_Position"], multiple_iterators=True
                )
            )

        except ValueError as e:
            print(e)
            continue

        if not fetched_row:
            continue

        if len(fetched_row) == 1:
            intersected_variants.extend(fetched_row)
            continue

        # Workings with multiple variants found.
        if (
            maf_rsid_col
            and (not pd.isna(variant[maf_rsid_col]))
            and (variant[maf_rsid_col] not in [None, np.nan, "novel", "."])
        ):
            fetched_row, failed_to_filter = match_rsid(variant, fetched_row, maf_rsid_col, dbNSFP_columns)
            if not failed_to_filter and len(fetched_row) == 1:
                intersected_variants.extend(fetched_row)
                continue

        if hgvsp_in_maf and not pd.isna(variant[hgvsp_in_maf].all()):
            fetched_row, failed_to_filter = match_hgvsp(variant, fetched_row, hgvsp_in_maf, dbNSFP_columns)
            if not failed_to_filter and len(fetched_row) == 1:
                intersected_variants.extend(fetched_row)
                continue

        if hgvsc_in_maf and not pd.isna(variant[hgvsc_in_maf].all()):
            fetched_row, failed_to_filter = match_hgvsc(variant, fetched_row, hgvsc_in_maf, dbNSFP_columns)
            if not failed_to_filter and len(fetched_row) == 1:
                intersected_variants.extend(fetched_row)
                continue

        if maf_enstid_col and variant[maf_enstid_col] not in [None, np.nan, "."]:
            fetched_row, failed_to_filter = match_enstid(variant, fetched_row, maf_enstid_col, dbNSFP_columns)
            if not failed_to_filter and len(fetched_row) == 1:
                intersected_variants.extend(fetched_row)
                continue

        if (variant["Reference_Allele"] in ["A", "G", "C", "T"]) and (
            variant["Tumor_Seq_Allele2"] in ["A", "G", "C", "T"]
        ):
            fetched_row, failed_to_filter = match_ref_alt(variant, fetched_row)
            if not failed_to_filter and len(fetched_row) == 1:
                intersected_variants.extend(fetched_row)
                continue

        intersected_variants.append(fetched_row[0])
        continue

    intersected_variants = pd.DataFrame(intersected_variants, columns=dbNSFP_columns)

    intersected_variants["pos(1-based)"] = intersected_variants["pos(1-based)"].apply(int)

    intersected_variants = pd.merge(
        chromosome_data,
        intersected_variants,
        left_on=["Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"],
        right_on=["pos(1-based)", "ref", "alt"],
        how="left",
        suffixes=[None, "_from_dbNSFP"],
    )

    intersected_variants = intersected_variants.drop_duplicates()

    return intersected_variants


def init_ray(host="0.0.0.0", n_jobs=6):
    ray.init(num_cpus=n_jobs, dashboard_host=host, ignore_reinit_error=True)
