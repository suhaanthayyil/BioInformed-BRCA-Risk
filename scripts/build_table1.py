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
    HEADLINE_MODEL,
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
    scored = set(pd.read_parquet(PROCESSED / "04_pathway_features.parquet", columns=["sample_id"])["sample_id"])
    rows = []
    for cohort in ALL_COHORTS:
        sub = frame[frame["cohort"].eq(cohort)].copy()
        analysis = sub[sub["sample_id"].isin(scored) & sub["os_days"].notna()].copy()
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
                    analysis["stage_ordinal"].map({1.0: "I", 2.0: "II", 3.0: "III", 4.0: "IV"}),
                    ["I", "II", "III", "IV"],
                ),
                "stage_recorded_n": int(analysis["stage_ordinal"].notna().sum()),
                "pam50_subtype_distribution_genefu": distribution(
                    analysis["pam50_subtype_genefu"], ["LumA", "LumB", "Her2", "Basal", "Normal"]
                ),
                "median_followup_months_iqr": median_iqr(followup / 30.4375),
                "os_events": int(pd.to_numeric(analysis["os_event"], errors="coerce").fillna(0).sum()),
                "notes": "Stage is AJCC anatomic stage I-IV (TCGA pathologic stage, METABRIC tumour stage, GSE20685 derived from TNM); GSE96058 releases no stage.",
            }
        )
    return pd.DataFrame(rows)


def renumber_tables() -> None:
    old = RESULTS / "Table_1_ml_external_validation.csv"
    new_external = RESULTS / "Table_2_external_validation.csv"
    if old.exists():
        shutil.copy2(old, new_external)
    old_h2h = RESULTS / "Table_2_head_to_head.csv"
    new_h2h = RESULTS / "Table_3_head_to_head.csv"
    if old_h2h.exists():
        shutil.copy2(old_h2h, new_h2h)


def denominator_crosswalk(table: pd.DataFrame) -> pd.DataFrame:
    pred = pd.read_csv(RESULTS / "06_external_model_predictions.csv")
    pred = pred[pred[HEADLINE_MODEL].notna() & pred[PAM50_MODEL].notna() & pred["os_days"].gt(0)]
    rows = []
    for _, r in table.iterrows():
        ev = pred[pred["cohort"].eq(r["cohort"])]
        rows.append(
            {
                "cohort": r["cohort"],
                "harmonized_n": r["initial_samples_in_harmonized_database"],
                "patient_characteristics_final_n": r["final_analysis_n"],
                "metric_evaluable_n": len(ev),
                "patient_characteristics_events": r["os_events"],
                "metric_evaluable_events": int(ev["os_event"].sum()),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    ensure_dirs()
    log_phase3("Stage 6 patient characteristics table started.")
    table = build_table()
    table.to_csv(RESULTS / "Table_1_patient_characteristics.csv", index=False)
    denominator_crosswalk(table).to_csv(RESULTS / "Table_1_denominator_crosswalk.csv", index=False)
    renumber_tables()
    log_phase3("Stage 6 patient characteristics table completed.")


if __name__ == "__main__":
    main()
