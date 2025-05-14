import glob
import importlib.resources as pkg_resources
import pickle
import re
from collections import defaultdict
from typing import Union

import numpy as np
import pandas as pd
import yaml
from Bio.Seq import reverse_complement
from sklearn.preprocessing import minmax_scale

from mutant import mutant_resources

GENCODE = {
    "ATA": "I",
    "ATC": "I",
    "ATT": "I",
    "ATG": "M",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "AAC": "N",
    "AAT": "N",
    "AAA": "K",
    "AAG": "K",
    "AGC": "S",
    "AGT": "S",
    "AGA": "R",
    "AGG": "R",
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CAC": "H",
    "CAT": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "GAC": "D",
    "GAT": "D",
    "GAA": "E",
    "GAG": "E",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    "TTC": "F",
    "TTT": "F",
    "TTA": "L",
    "TTG": "L",
    "TAC": "Y",
    "TAT": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGC": "C",
    "TGT": "C",
    "TGA": "*",
    "TGG": "W",
}
AA_TO_CODONS = defaultdict(list)
{AA_TO_CODONS[v].append(k) for k, v in GENCODE.items()}

GENCODE_REVERSED = {reverse_complement(k): v for k, v in GENCODE.items()}
AA_TO_CODONS_REVERSED = defaultdict(list)
{AA_TO_CODONS_REVERSED[v].append(k) for k, v in GENCODE_REVERSED.items()}


VARIANT_INTERPRETATION_MAP = {
    "T": 0,
    "D": 1,
    "A": 1,
    "N": 0,
    "P": 0,
    "H": 1,
    "M": 1,
    "L": 0,
}


def process_multiple_transcript(s):
    """
    Split multiple transcripts stored in one line in dbNSFP annotated clinvar VCF
    :param s: row of pd.DataFrame
    :return: single transcript
    """
    first_transcript_data = re.search(r"^(.*?);", str(s) + ";")
    if first_transcript_data is not None:
        return first_transcript_data.group(1)
    else:
        return np.nan


def nan_replace(s) -> Union[str, None]:
    """
    Replace dbNSFP nan symbol to np.nan
    :param s: row of pd.DataFrame
    :return: string with replaced . to np.nan
    """
    try:
        if re.search(r"^\.\B|^-\B", str(s)):
            return None
        return s
    except AttributeError:
        return s


def process_dbnsfp_annotated_file_for_mutant(df):
    """
    Transform dbNSFP annotated MAF/VCF file for MutAnt annotation
    :param df: pd.DataFrame - annotated with dbNSFP MAF or VCF file
    :return: df ready for predictions
    """

    def blosum_score(x, blosum):
        try:
            return blosum[(x["aaref"], x["aaalt"])]
        except KeyError:
            try:
                return blosum[(x["aaalt"], x["aaref"])]
            except Exception:
                return None

    # SiPhy_29way_pi
    def process_siphy_score(s) -> str:
        """
        Fuction to process SiPhy_29way_pi column in annotated clinvar VCF
        :param s: row of pd.DataFrame
        :return: single string
        """
        siphy_first_el = re.search(r"^(.*?):", s + ":")
        if siphy_first_el is not None:
            return siphy_first_el.group(1)
        else:
            return s

    df = df.applymap(process_multiple_transcript)
    df = df.applymap(process_siphy_score)
    df = df.applymap(nan_replace)

    str_cols = []
    for col in df.columns:
        try:
            df[col].sample(n=int(df.shape[0] / 2), replace=True).fillna(-10).astype(
                float
            )
        except Exception:
            str_cols.append(col)

    df = df.drop(columns=str_cols)

    return df


def snp_finder(row: pd.Series) -> dict:
    """
    Change variant info from DNP to SNP if available.

    :param row: pd.Series, one variant info from maf.
    :return: dict, summary for variant.
    """
    # Set defaults.
    ref = row["Reference_Allele"]
    alt = row["Tumor_Seq_Allele2"]
    ref_1_codon = row["ref_1"]
    ref_2_codon = row["ref_2"]
    alt_1_codon = row["alt_1"]
    alt_2_codon = row["alt_2"]
    start_position = row["Start_Position"]
    variant_type = row["Variant_Type"]
    dnp_converted = False
    ref_saved, alt_saved, start_pos_saved, end_pos_saved = (
        np.nan,
        np.nan,
        np.nan,
        np.nan,
    )

    # Check strand for correct gencode.
    gencode = GENCODE
    aa_to_codons = AA_TO_CODONS
    if row["STRAND_VEP"] == -1.0:
        gencode = GENCODE_REVERSED
        aa_to_codons = AA_TO_CODONS_REVERSED

    # if 1 Codon affected.
    if (ref_2_codon == "" and alt_2_codon == "") or (
        pd.isna(ref_2_codon) and pd.isna(alt_2_codon)
    ):

        # Generate all possible AA codons.
        alt_codons_1 = aa_to_codons[gencode[alt_1_codon]]

        # Find SNP codon and position.
        snp_codon = None
        snp_codon_position = None
        save = 3
        for case in alt_codons_1:
            diffs = [
                ref_1_codon[0] != case[0],
                ref_1_codon[1] != case[1],
                ref_1_codon[2] != case[2],
            ]
            if diffs.count(True) < save:
                save = diffs.count(True)
            if diffs.count(True) == 1:
                snp_codon = case
                snp_codon_position = diffs.index(True)
                break

        # Get start DNP position in codon.
        first_dnp_nucleotide_position = ref_1_codon.find(row["Reference_Allele"])

        # Set if found snp codon.
        if snp_codon is not None:
            ref = ref_1_codon[snp_codon_position]
            alt = snp_codon[snp_codon_position]
            start_position = row["Start_Position"] - (
                first_dnp_nucleotide_position - snp_codon_position
            )
            variant_type = "SNP"
            dnp_converted = True
        if variant_type == "DNP" and save < 2:
            raise ValueError("DNP -> SNP missmatch codons difference.")

    # if 2 Codons affected.
    elif (
        ref_2_codon != ""
        and alt_2_codon != ""
        and pd.notna(ref_2_codon)
        and pd.notna(alt_2_codon)
    ):
        # Because of 2 alt nucleotides are on the bondary of two codons.
        # Example: agCCcg
        if gencode[ref_1_codon] != gencode[alt_1_codon]:
            ref = ref_1_codon[-1:]
            alt = alt_1_codon[-1:]
            start_position = row["Start_Position"]
            variant_type = "SNP"
            dnp_converted = True
        elif gencode[ref_2_codon] != gencode[alt_2_codon]:
            ref = ref_2_codon[:1]
            alt = alt_2_codon[:1]
            start_position = row["End_Position"]
            variant_type = "SNP"
            dnp_converted = True

    if dnp_converted:
        ref_saved = row["Reference_Allele"]
        alt_saved = row["Tumor_Seq_Allele2"]
        start_pos_saved = row["Start_Position"]
        end_pos_saved = row["End_Position"]

    return {
        "Reference_Allele": ref,
        "Tumor_Seq_Allele2": alt,
        "Start_Position": start_position,
        "End_Position": start_position,
        "Variant_Type": variant_type,
        "DNP_Converted": dnp_converted,
        "Reference_Allele_DNP_Saved": ref_saved,
        "Tumor_Seq_Allele2_DNP_Saved": alt_saved,
        "Start_Position_DNP_Saved": start_pos_saved,
        "End_Position_DNP_Saved": end_pos_saved,
    }


def dnp_to_snp_formater(maf: pd.DataFrame) -> pd.DataFrame:
    """
    Change DNP variants to SNP if available.
    """

    # Filter DNP and non empty codons.
    maf_slice_with_dnp = maf[
        (maf["Variant_Type"] == "DNP")
        & (~maf["Codons"].empty)
        & (maf["Codons"].notna())
    ]

    maf["Reference_Allele_DNP_Saved"] = np.nan
    maf["Tumor_Seq_Allele2_DNP_Saved"] = np.nan
    maf["Start_Position_DNP_Saved"] = np.nan
    maf["End_Position_DNP_Saved"] = np.nan
    maf["DNP_Converted"] = False

    if maf_slice_with_dnp.empty:
        return maf

    # Extraction of ref and alt codon(s).
    maf_slice_with_dnp["ref_all"] = (
        maf_slice_with_dnp["Codons"].str.upper().str.extract(r"([ACGTacgt]+)")
    )

    maf_slice_with_dnp["alt_all"] = (
        maf_slice_with_dnp["Codons"].str.upper().str.extract(r"([^/]*$)")
    )

    # Extraction of target codons (sometimes we can see 4 codons - two ref and two alt.
    def extract_codons(row):
        found_ref, found_alt = re.findall(
            r"([ACGTacgt]{3,6})\/([ACGTacgt]{3,6})", row["Codons"]
        )[0]
        found_alt = found_alt.upper()
        found_ref = found_ref.upper()
        if row["STRAND_VEP"] == -1.0:
            found_alt, found_ref = map(reverse_complement, [found_alt, found_ref])
        return {
            "alt_1": found_alt[:3],
            "alt_2": found_alt[3:],
            "ref_1": found_ref[:3],
            "ref_2": found_ref[3:],
        }

    codons_info = maf_slice_with_dnp.apply(extract_codons, axis=1, result_type="expand")
    maf_slice_with_dnp[codons_info.columns] = codons_info

    snp_info = maf_slice_with_dnp.apply(snp_finder, axis=1, result_type="expand")
    maf.loc[
        maf_slice_with_dnp.index,
        snp_info.columns,
    ] = snp_info

    return maf


def dnp_revert(row: pd.Series) -> dict:
    """
    Back saved DNP information.
    """
    if row["DNP_Converted"]:
        return {
            "Reference_Allele": row["Reference_Allele_DNP_Saved"],
            "Tumor_Seq_Allele2": row["Tumor_Seq_Allele2_DNP_Saved"],
            "Start_Position": int(row["Start_Position_DNP_Saved"]),
            "End_Position": int(row["End_Position_DNP_Saved"]),
            "Variant_Type": "DNP",
        }
    return {
        "Reference_Allele": row["Reference_Allele"],
        "Tumor_Seq_Allele2": row["Tumor_Seq_Allele2"],
        "Start_Position": row["Start_Position"],
        "End_Position": row["End_Position"],
        "Variant_Type": row["Variant_Type"],
    }


def dnp_columns_manipulator(maf: pd.DataFrame, dnp_convert_info: bool) -> pd.DataFrame:
    revert_info = maf.apply(dnp_revert, axis=1, result_type="expand")
    maf[revert_info.columns] = revert_info

    maf = maf.drop(
        columns=[
            "Reference_Allele_DNP_Saved",
            "Tumor_Seq_Allele2_DNP_Saved",
            "Start_Position_DNP_Saved",
            "End_Position_DNP_Saved",
        ]
    )
    if not dnp_convert_info:
        maf = maf.drop(columns=["DNP_Converted"])

    return maf


def add_features_info(df: pd.DataFrame, features_from_model: list) -> pd.DataFrame:
    """
    Adding features info.

    :param df: Annotated df with MutAnt score.
    :return: df with features inference info.
    """

    # Open yaml file with MutAnt's features.
    with pkg_resources.path(mutant_resources, "features.yaml") as p:
        with open(p, "r") as f:
            features = yaml.safe_load(f)

    freq_features = [
        c for c in features_from_model if c in features["frequency_features"]
    ]
    tools_features = [c for c in features_from_model if c in features["tools_features"]]
    evolution_features = [
        c for c in features_from_model if c in features["evolution_features"]
    ]

    def basis_extract(row):
        result = []
        if not row[freq_features].isnull().values.all():
            result.append("frequencies")
        if not row[tools_features].isnull().values.all():
            result.append("tools")
        if not row[evolution_features].isnull().values.all():
            result.append("evolution")

        if not result:
            return "no info"

        return ", ".join(result)

    df = df.replace(".", None)
    df["MutAnt_feature_availability"] = df.apply(
        basis_extract, axis=1, result_type="expand"
    )
    df.loc[
        (df["MutAnt_feature_availability"] == "no info")
        | (df["MutAnt_feature_availability"] == "evolution"),
        "MutAnt_score",
    ] = None
    return df