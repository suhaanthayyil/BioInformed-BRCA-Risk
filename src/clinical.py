"""Clinical feature lists and stage parsing."""

from __future__ import annotations

import re

import numpy as np
import pandas as pd

PATHWAY_FEATURES = [
    "Pathway_Immune",
    "Pathway_Proliferation",
    "Pathway_DNA_Repair",
    "Pathway_Metabolism",
    "Pathway_Stromal_EMT",
    "Pathway_Apoptosis_Stress",
    "Pathway_Hormone",
]

PRIMARY_FEATURES = list(PATHWAY_FEATURES)

ROMAN = {"I": 1.0, "II": 2.0, "III": 3.0, "IV": 4.0}


def tnm_to_stage(t: str, n: str, m: str) -> float:
    if m == "1":
        return 4.0
    if t == "4" or n == "3":
        return 3.0
    if n == "2" or (t == "3" and n == "1"):
        return 3.0
    if (t == "3" and n == "0") or (t == "2" and n == "1") or (t == "2" and n == "0") or (t in {"0", "1"} and n == "1"):
        return 2.0
    if t in {"1", "IS"} and n == "0":
        return 1.0
    return np.nan


def stage_to_ordinal(value) -> float:
    if value is None or pd.isna(value):
        return np.nan
    text = str(value).strip().upper()
    hit = re.fullmatch(r"(?:STAGE\s*)?(IV|III|II|I)[ABC]?", text)
    if hit:
        return ROMAN[hit.group(1)]
    hit = re.fullmatch(r"([1-4])(?:\.0)?", text)
    if hit:
        return float(hit.group(1))
    hit = re.fullmatch(r"T(IS|\d)(?:\.0)?N(\d)(?:\.0)?M(\d)(?:\.0)?", text)
    if hit:
        return tnm_to_stage(*hit.groups())
    return np.nan
