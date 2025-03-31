import csv
import os
import pickle
from typing import Literal, Union

import numpy as np
import pandas as pd
from dbnsfp_annotation.dbnsfp_annotation import annotate_maf_dbnsfp, annotate_vcf_dbnsfp
from dbnsfp_annotation.utils import DATABASE_PATH, GENE_BASE_PATH
from dotenv import load_dotenv

from mutant.utils import (
    add_features_info,
    dnp_columns_manipulator,
    dnp_to_snp_formater,
    process_dbnsfp_annotated_file_for_mutant,
)

load_dotenv()

MODEL_PATH = "v45a_academic_mutant.pickle"


def add_mutant_annotation(
    obj: Union[pd.DataFrame, str],
    file_type: Literal["maf", "vcf"],
    genome_version: Literal["hg19", "hg38"],
    path_to_database: str = DATABASE_PATH,
    gene_base_path: str = GENE_BASE_PATH,
    development_mode: bool = False,
    dnp_convert_info: bool = False,
    n_jobs: int = 2,
    model_path: str = os.environ.get("MODEL_PATH") or MODEL_PATH,
    method: str = "predict_proba",
    annotated: bool = False,
) -> pd.DataFrame:
    """
    Adds MutAnt_score to maf/vcf file as a new column.

    :param obj: Path to a maf/vcf file, or a maf/vcf file loaded into pd.DataFrame.
    :param file_type: Either 'maf' or 'vcf'.
    :param genome_version: Either 'hg19' or 'hg38'.
    :param path_to_database: Path to the gzipped database file.
    :param development_mode: If True, it will return a DataFrame with all columns from dbNSFP.
                                If False, only the MutAnt score will be added.
    :param dnp_convert_info: If True, it will return information about DNP converted SNP.
                                If False, it will silently replace fields.
    :param n_jobs: Maximum number of cores for processing.
    :param model_path: Absolute path to MutAnt model_academic.lgb.
    :param cohort: Cohort in TCGA format, used to filter out genes with low expression.
    :param percentile: If cohort is specified, it is used to set the threshold.
                        Genes with expression lower than this percentile will be considered non-pathogenic.
    :param method: Use model.predict_proba() if 'predict_proba', otherwise use model.predict().
    :param annotated: True means that the maf file is already annotated with dbNSFP database,
                        so annotation will not be repeated. Genome_version should be hg38.
    :return: DataFrame with new columns 'MutAnt_score' and 'MutAnt_oncogene_supressor_classification'.
    :rtype: pd.DataFrame

    """
    if file_type == "maf":

        if isinstance(obj, str):
            maf = pd.read_csv(obj, sep="\t", comment="#")
            initial_cols = maf.columns.tolist()

        else:
            maf = obj
            initial_cols = maf.columns.tolist()

        initial_cols = [
            i
            for i in initial_cols
            if i
            not in [
                "MutAnt_score",
                "MutAnt_feature_availability",
            ]
        ]

        # DNP to SNP.
        dnp_triggered = False
        if "DNP" in maf["Variant_Type"].unique():
            dnp_triggered = True
            maf = dnp_to_snp_formater(maf)

        if not annotated:
            df = annotate_maf_dbnsfp(
                database_path=path_to_database,
                gene_base_path=gene_base_path,
                maf=maf,
                genome_version=genome_version,
                n_jobs=n_jobs,
            )
        else:
            df = maf

        X = process_dbnsfp_annotated_file_for_mutant(df)

        with open(model_path, "rb") as f:
            model_features = pickle.load(f)
        model, features = model_features.values()

        if method == "predict_proba":
            df["MutAnt_score"] = model.predict_proba(X[features].astype(float))[:, 1]
        if method == "predict":
            df["MutAnt_score"] = model.predict(X[features].astype(float))

        if dnp_triggered:
            df = dnp_columns_manipulator(df, dnp_convert_info)

        df = add_features_info(df, features)

        if development_mode:
            return df

        output_columns = [
            "MutAnt_score",
            "MutAnt_feature_availability",
        ]
        if dnp_convert_info:
            output_columns.append("DNP_Converted")

        return df.loc[
            :,
            initial_cols + output_columns,
        ]

    elif file_type == "vcf":

        if isinstance(obj, str):
            with open(obj, "r") as vcf_file:
                vcf_reader = csv.reader(vcf_file, delimiter="\t")
                vcf_filtered = (i for i in vcf_reader if not i[0].startswith("##"))

                df_columns = next(vcf_filtered)

                vcf = pd.DataFrame(vcf_filtered, columns=df_columns)
                vcf_cols = vcf.columns.tolist()
        else:
            vcf = obj
            vcf_cols = vcf.columns.tolist()

        vcf_cols = [
            i
            for i in vcf_cols
            if i
            not in [
                "MutAnt_score",
                "MutAnt_feature_availability",
            ]
        ]

        df = annotate_vcf_dbnsfp(
            database_path=path_to_database,
            gene_base_path=gene_base_path,
            vcf=vcf,
            genome_version=genome_version,
            n_jobs=n_jobs,
        )

        X = process_dbnsfp_annotated_file_for_mutant(df)

        with open(model_path, "rb") as f:
            model_features = pickle.load(f)
        model, features = model_features.values()

        if method == "predict_proba":
            df["MutAnt_score"] = model.predict_proba(X[features].astype(float))[:, 1]
        if method == "predict":
            df["MutAnt_score"] = model.predict(X[features].astype(float))

        df = add_features_info(df, features)

        if development_mode:
            return df

        return df.loc[
            :,
            vcf_cols
            + [
                "MutAnt_score",
                "MutAnt_feature_availability",
            ],
        ]