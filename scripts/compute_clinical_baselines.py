#!/usr/bin/env python3
"""Compute Phase 3 clinical baseline scores and survival metrics."""

from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from src.baselines import (
    ONCOTYPE_GENES,
    PROLIFERATION_GENES,
    commercial_formula_status,
    expression_wide_to_samples_by_gene,
    mammaprint_status,
    oncotype_dx_surrogate,
    pam50_ror_surrogate,
    percentile_0_100,
    samples_by_gene_from_tcga,
)
from src.survival import (
    brier_score_at,
    bootstrap_c_index_ci,
    cox_hr_high_low,
    harrell_c_index,
    time_dependent_auc,
    uno_c_index,
)


DATA = REPO_ROOT / "data"
PROCESSED = DATA / "processed"
RESULTS = REPO_ROOT / "results"
REPORTS = RESULTS / "reports"
LOGS = REPO_ROOT / "logs"
DB_PATH = PROCESSED / "unified_cohorts.duckdb"

BOOTSTRAP_ITER = 1000
SEED = 42


def timestamp() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def log(message: str) -> None:
    LOGS.mkdir(exist_ok=True)
    with (LOGS / "agent_log.md").open("a") as fh:
        fh.write(f"\n- {timestamp()} {message}\n")


def load_harmonized_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    con = duckdb.connect(str(DB_PATH), read_only=True)
    try:
        samples = con.execute("select * from samples").fetchdf()
        survival = con.execute("select * from survival").fetchdf()
    finally:
        con.close()
    samples["sample_id"] = samples["sample_id"].astype(str)
    survival["sample_id"] = survival["sample_id"].astype(str)
    return samples, survival


def needed_genes() -> list[str]:
    genes = set(PROLIFERATION_GENES)
    for items in ONCOTYPE_GENES.values():
        genes.update(items)
    return sorted(genes)


def load_baseline_expression(cohort: str, sample_ids: list[str]) -> pd.DataFrame:
    genes = needed_genes()
    if cohort == "TCGA-BRCA":
        return samples_by_gene_from_tcga(
            str(DATA / "01_tcga_expression_normalized.csv"),
            str(PROCESSED / "02_tcga_feature_matrix.csv"),
            keep_genes=genes,
        ).reindex(sample_ids)
    if cohort == "GSE96058":
        expr = pd.read_parquet(PROCESSED / "02_gse96058_expression.parquet")
        return expression_wide_to_samples_by_gene(expr, sample_ids=sample_ids, keep_genes=genes)
    if cohort == "METABRIC":
        expr = pd.read_parquet(PROCESSED / "03_metabric_expression.parquet")
        return expression_wide_to_samples_by_gene(expr, sample_ids=sample_ids, keep_genes=genes)
    if cohort == "GSE20685":
        expr = pd.read_parquet(PROCESSED / "04_gse20685_expression.parquet")
        return expression_wide_to_samples_by_gene(expr, sample_ids=sample_ids, keep_genes=genes)
    return pd.DataFrame(index=pd.Index([], name="sample_id"))


def tumor_size_series(cohort: str, sample_ids: list[str]) -> pd.Series | None:
    if cohort == "GSE96058":
        clinical = pd.read_csv(DATA / "clinical" / "01_gse96058_clinical.csv")
        return clinical.set_index("sample_id").reindex(sample_ids)["tumor_size"]
    if cohort == "METABRIC":
        clinical = pd.read_parquet(PROCESSED / "03_metabric_clinical.parquet")
        return clinical.set_index("sample_id").reindex(sample_ids)["TUMOR_SIZE"]
    return None


def score_to_frame(result, cohort: str, samples: pd.DataFrame, score_0_100: pd.Series | None = None) -> pd.DataFrame:
    if result.score.empty:
        return pd.DataFrame(
            [
                {
                    "cohort": cohort,
                    "sample_id": pd.NA,
                    "baseline": result.baseline,
                    "score": np.nan,
                    "score_0_100": np.nan,
                    "method_label": result.method_label,
                    "status": result.status,
                    "note": result.note,
                }
            ]
        )
    frame = pd.DataFrame(
        {
            "cohort": cohort,
            "sample_id": result.score.index.astype(str),
            "baseline": result.baseline,
            "score": result.score.astype(float).to_numpy(),
            "score_0_100": (score_0_100 if score_0_100 is not None else percentile_0_100(result.score)).reindex(result.score.index).to_numpy(),
            "method_label": result.method_label,
            "status": result.status,
            "note": result.note,
        }
    )
    return frame.merge(samples[["sample_id", "er_status", "tnbc_flag", "intrinsic_subtype_pam50_published"]], on="sample_id", how="left")


def compute_scores(samples: pd.DataFrame) -> pd.DataFrame:
    all_scores: list[pd.DataFrame] = []
    for cohort, cohort_samples in samples.groupby("cohort", sort=True):
        sample_ids = cohort_samples["sample_id"].astype(str).tolist()
        expr = load_baseline_expression(cohort, sample_ids)
        if not expr.empty:
            expr = expr.reindex(sample_ids)

        subtype = cohort_samples.set_index("sample_id")["intrinsic_subtype_pam50_published"].reindex(sample_ids)
        pam50 = pam50_ror_surrogate(subtype, expr, tumor_size_series(cohort, sample_ids))
        all_scores.append(score_to_frame(pam50, cohort, cohort_samples))

        if expr.empty or expr.shape[1] < 8:
            note = f"Oncotype surrogate unavailable in {cohort}: expression genes missing after harmonization."
            result = commercial_formula_status("Oncotype_DX_21_gene_surrogate")
            result = type(result)(result.score, result.baseline, result.method_label, "unavailable", note)
            all_scores.append(score_to_frame(result, cohort, cohort_samples))
        else:
            oncotype = oncotype_dx_surrogate(expr)
            frame = score_to_frame(oncotype, cohort, cohort_samples)
            if "er_status" in frame:
                missing_er = frame["er_status"].isna().all()
                frame.loc[frame["er_status"] != "Positive", ["score", "score_0_100"]] = np.nan
                frame["note"] = frame["note"] + " Restricted to ER-positive samples for performance evaluation."
                if missing_er:
                    frame["status"] = "not_evaluable_missing_er_status"
                    frame["note"] = (
                        frame["note"]
                        + " No ER status is available for this cohort, so the ER-positive Oncotype analysis is not evaluable."
                    )
            all_scores.append(frame)

        all_scores.append(score_to_frame(mammaprint_status(), cohort, cohort_samples))
        all_scores.append(score_to_frame(commercial_formula_status("EndoPredict"), cohort, cohort_samples))
        all_scores.append(score_to_frame(commercial_formula_status("Breast_Cancer_Index"), cohort, cohort_samples))

    scores = pd.concat(all_scores, ignore_index=True)
    return scores


def evaluate_scores(scores: pd.DataFrame, samples: pd.DataFrame, survival: pd.DataFrame) -> pd.DataFrame:
    merged = (
        scores[scores["sample_id"].notna()]
        .merge(survival, on="sample_id", how="left")
        .merge(samples[["sample_id", "cohort", "tnbc_flag"]], on=["sample_id", "cohort"], how="left", suffixes=("", "_sample"))
    )
    rows = []
    for (cohort, baseline), group in merged.groupby(["cohort", "baseline"], sort=True):
        tnbc_flag = group["tnbc_flag"].astype("boolean").fillna(False).astype(bool)
        for subgroup_name, mask in {
            "overall": pd.Series(True, index=group.index),
            "tnbc": tnbc_flag,
            "non_tnbc": ~tnbc_flag,
        }.items():
            sub = group.loc[mask].copy()
            sub = sub[sub["score"].notna() & sub["os_days"].notna() & sub["os_event"].notna()]
            n = len(sub)
            events = int(pd.to_numeric(sub["os_event"], errors="coerce").fillna(0).sum()) if n else 0
            if n < 10 or events < 2:
                row = {
                    "cohort": cohort,
                    "baseline": baseline,
                    "endpoint": "OS",
                    "subgroup": subgroup_name,
                    "n": n,
                    "events": events,
                    "harrell_c": np.nan,
                    "harrell_c_ci_low": np.nan,
                    "harrell_c_ci_high": np.nan,
                    "uno_c": np.nan,
                    "auc_3y": np.nan,
                    "auc_5y": np.nan,
                    "auc_10y": np.nan,
                    "brier_5y": np.nan,
                    "hr": np.nan,
                    "hr_ci_low": np.nan,
                    "hr_ci_high": np.nan,
                    "hr_p": np.nan,
                    "status": "too_few_events",
                    "note": "Fewer than 10 samples or 2 events after score/endpoint filtering.",
                }
            else:
                c = harrell_c_index(sub["os_days"], sub["os_event"], sub["score"])
                ci_low, ci_high = bootstrap_c_index_ci(
                    sub["os_days"],
                    sub["os_event"],
                    sub["score"],
                    n_bootstrap=BOOTSTRAP_ITER,
                    seed=SEED,
                )
                row = {
                    "cohort": cohort,
                    "baseline": baseline,
                    "endpoint": "OS",
                    "subgroup": subgroup_name,
                    "n": n,
                    "events": events,
                    "harrell_c": c,
                    "harrell_c_ci_low": ci_low,
                    "harrell_c_ci_high": ci_high,
                    "uno_c": uno_c_index(sub["os_days"], sub["os_event"], sub["score"]),
                    **time_dependent_auc(sub["os_days"], sub["os_event"], sub["score"], years=(3, 5, 10)),
                    **brier_score_at(sub["os_days"], sub["os_event"], sub["score"], years=(5,)),
                    **cox_hr_high_low(sub["os_days"], sub["os_event"], sub["score"]),
                    "status": "ok",
                    "note": str(sub["note"].dropna().iloc[0]) if sub["note"].notna().any() else "",
                }
            rows.append(row)
    return pd.DataFrame(rows)


def write_report(scores: pd.DataFrame, metrics: pd.DataFrame) -> None:
    REPORTS.mkdir(parents=True, exist_ok=True)
    ok_scores = scores[scores["status"] == "ok"]
    unavailable = scores[scores["status"] != "ok"][["cohort", "baseline", "status", "note"]].drop_duplicates()
    lines = [
        "# Phase 3 Clinical Baseline QC",
        "",
        f"Generated: {timestamp()}",
        "",
        "## Score Inventory",
        "",
        ok_scores.assign(scored=ok_scores["score"].notna())
        .groupby(["cohort", "baseline"])["scored"]
        .sum()
        .rename("n_scored")
        .reset_index()
        .to_markdown(index=False),
        "",
        "## Unavailable / Surrogate Notes",
        "",
        unavailable.to_markdown(index=False),
        "",
        "## Overall OS Metrics",
        "",
        metrics[(metrics["subgroup"] == "overall")][
            ["cohort", "baseline", "n", "events", "harrell_c", "harrell_c_ci_low", "harrell_c_ci_high", "uno_c", "status"]
        ].to_markdown(index=False, floatfmt=".3f"),
        "",
        "## TNBC OS Metrics",
        "",
        metrics[(metrics["subgroup"] == "tnbc")][
            ["cohort", "baseline", "n", "events", "harrell_c", "harrell_c_ci_low", "harrell_c_ci_high", "status"]
        ].to_markdown(index=False, floatfmt=".3f"),
        "",
        "PAM50 values are marked as `PAM50_ROR_surrogate` unless official genefu scoring is available. "
        "MammaPrint, EndoPredict, and BCI are not fabricated when exact public formulae are absent.",
        "",
    ]
    (REPORTS / "clinical_baselines_qc.md").write_text("\n".join(lines))


def main() -> None:
    np.random.seed(SEED)
    log("Phase 3 clinical baseline scoring started.")
    samples, survival = load_harmonized_tables()
    scores = compute_scores(samples)

    pam50 = scores[scores["baseline"].eq("PAM50_ROR_surrogate")]
    oncotype = scores[scores["baseline"].eq("Oncotype_DX_21_gene_surrogate")]
    mammaprint = scores[scores["baseline"].eq("MammaPrint_gene70")]
    other = scores[scores["baseline"].isin(["EndoPredict", "Breast_Cancer_Index"])]

    pam50.to_parquet(PROCESSED / "baselines_pam50.parquet", index=False)
    oncotype.to_parquet(PROCESSED / "baselines_oncotype.parquet", index=False)
    mammaprint.to_parquet(PROCESSED / "baselines_mammaprint.parquet", index=False)
    other.to_parquet(PROCESSED / "baselines_other.parquet", index=False)

    metrics = evaluate_scores(scores, samples, survival)
    RESULTS.mkdir(exist_ok=True)
    metrics.to_csv(RESULTS / "Table_S1_clinical_baselines.csv", index=False)
    metrics.to_latex(RESULTS / "Table_S1_clinical_baselines.tex", index=False, float_format="%.3f")

    metadata = {
        "generated_at": timestamp(),
        "bootstrap_iterations": BOOTSTRAP_ITER,
        "seed": SEED,
        "score_files": {
            "pam50": "data/processed/baselines_pam50.parquet",
            "oncotype": "data/processed/baselines_oncotype.parquet",
            "mammaprint": "data/processed/baselines_mammaprint.parquet",
            "other": "data/processed/baselines_other.parquet",
        },
        "metric_file": "results/Table_S1_clinical_baselines.csv",
        "surrogate_policy": "No proprietary assay score was fabricated; unavailable exact formulas are labeled unavailable.",
    }
    (PROCESSED / "clinical_baselines.metadata.json").write_text(json.dumps(metadata, indent=2))
    write_report(scores, metrics)
    log("Phase 3 clinical baseline scoring completed.")


if __name__ == "__main__":
    main()
