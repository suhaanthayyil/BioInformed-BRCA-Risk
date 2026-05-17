#!/usr/bin/env python3
"""Stage 6: patient characteristics table and BMC table renumbering."""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from scripts.phase3_common import (  # noqa: E402
    ALL_COHORTS,
    PAM50_MODEL,
    PROCESSED,
    RESULTS,
    ensure_dirs,
    load_samples_survival,
    log_phase3,
)


def median_iqr(series: pd.Series, digits: int = 1) -> str:
    values = pd.to_numeric(series, errors="coerce").dropna()
    if values.empty:
        return "NA"
    q1, median, q3 = np.quantile(values, [0.25, 0.5, 0.75])
    return f"{median:.{digits}f} ({q1:.{digits}f}-{q3:.{digits}f})"


def positive_percent(series: pd.Series) -> str:
    clean = series.dropna().astype(str).str.lower()
    known = clean[clean.ne("") & ~clean.isin(["none", "nan", "unknown", "not_available"])]
    if known.empty:
        return "NA"
    positives = known.str.contains("positive|pos|amplified|yes", regex=True).sum()
    return f"{positives} ({100 * positives / len(known):.1f}%)"


def distribution(series: pd.Series, order: list[str] | None = None) -> str:
    clean = series.dropna().astype(str)
    clean = clean[~clean.str.lower().isin(["none", "nan", "unknown", "not_available", ""])]
    if clean.empty:
        return "NA"
    counts = clean.value_counts()
    if order:
        ordered = [item for item in order if item in counts.index]
        ordered += [item for item in counts.index if item not in ordered]
        counts = counts.loc[ordered]
    total = counts.sum()
    return "; ".join(f"{idx}: {int(val)} ({100 * val / total:.1f}%)" for idx, val in counts.items())


def build_table() -> pd.DataFrame:
    frame = load_samples_survival()
    pam = pd.read_parquet(PROCESSED / "baselines_pam50.parquet")
    pam = pam[pam["baseline"].eq(PAM50_MODEL) & pam["status"].eq("ok")][
        ["sample_id", "pam50_subtype_genefu", "n_pam50_genes_input"]
    ]
    frame = frame.merge(pam, on="sample_id", how="left")
    rows = []
    for cohort in ALL_COHORTS:
        sub = frame[frame["cohort"].eq(cohort)].copy()
        survival_valid = sub["os_days"].notna() & (sub["os_days"] > 0)
        pam50_valid = sub["pam50_subtype_genefu"].notna()
        analysis = sub[survival_valid & pam50_valid].copy()
        followup = pd.to_numeric(analysis["last_followup_days"], errors="coerce").fillna(
            pd.to_numeric(analysis["os_days"], errors="coerce")
        )
        rows.append(
            {
                "cohort": cohort,
                "initial_samples_in_harmonized_database": len(sub),
                "final_analysis_n": len(analysis),
                "age_at_diagnosis_median_iqr_years": median_iqr(analysis["age_at_dx"]),
                "tumor_size_median_iqr": "NA",
                "er_positive": positive_percent(analysis["er_status"]),
                "pr_positive": positive_percent(analysis["pr_status"]),
                "her2_positive": positive_percent(analysis["her2_status"]),
                "lymph_node_positive": "NA",
                "stage_distribution": distribution(
                    analysis["stage_tnm"],
                    ["Stage I", "Stage IA", "Stage IB", "Stage II", "Stage IIA", "Stage IIB", "Stage III", "Stage IV"],
                ),
                "pam50_subtype_distribution_genefu": distribution(
                    analysis["pam50_subtype_genefu"], ["LumA", "LumB", "Her2", "Basal", "Normal"]
                ),
                "median_followup_months_iqr": median_iqr(followup / 30.4375),
                "os_events": int(pd.to_numeric(analysis["os_event"], errors="coerce").fillna(0).sum()),
                "treatment_received": distribution(analysis["treatment_summary"]),
                "notes": "Treatment, tumor size, and lymph node fields are unavailable for some public cohorts.",
            }
        )
    return pd.DataFrame(rows)


def renumber_tables() -> None:
    legacy = RESULTS / "Table_1_LEGACY_external_validation.csv"
    old = RESULTS / "Table_1_ml_external_validation.csv"
    new_external = RESULTS / "Table_2_external_validation.csv"
    if old.exists():
        shutil.copy2(old, legacy)
        shutil.copy2(old, new_external)
    old_h2h = RESULTS / "Table_2_head_to_head.csv"
    new_h2h = RESULTS / "Table_3_head_to_head.csv"
    if old_h2h.exists():
        shutil.copy2(old_h2h, new_h2h)


def main() -> None:
    ensure_dirs()
    log_phase3("Stage 6 patient characteristics table started.")
    table = build_table()
    table.to_csv(RESULTS / "Table_1_patient_characteristics.csv", index=False)
    renumber_tables()
    log_phase3("Stage 6 patient characteristics table completed.")


if __name__ == "__main__":
    main()
