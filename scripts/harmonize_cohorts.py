#!/usr/bin/env python3
"""Build the Phase 2 harmonized cohort database.

The database stores compact sample/survival/QC tables and expression-long
views over processed wide expression files. The views avoid materializing
hundreds of millions of expression rows while preserving the requested schema.
"""

from __future__ import annotations

import json
from pathlib import Path

import duckdb
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = REPO_ROOT / "data"
PROCESSED = DATA_ROOT / "processed"
CLINICAL = DATA_ROOT / "clinical"
RESULTS = REPO_ROOT / "results" / "reports"
DB_PATH = PROCESSED / "unified_cohorts.duckdb"


SAMPLE_COLUMNS = [
    "sample_id",
    "cohort",
    "patient_id",
    "age_at_dx",
    "sex",
    "race",
    "ethnicity",
    "stage_tnm",
    "grade",
    "er_status",
    "pr_status",
    "her2_status",
    "intrinsic_subtype_pam50_published",
    "tnbc_flag",
    "primary_site",
    "treatment_summary",
]

SURVIVAL_COLUMNS = [
    "sample_id",
    "os_days",
    "os_event",
    "dfs_days",
    "dfs_event",
    "last_followup_days",
]


def _norm_status(value: object) -> str | None:
    if pd.isna(value):
        return None
    text = str(value).strip().lower()
    if text in {"1", "1.0", "positive", "pos", "positve", "yes"}:
        return "Positive"
    if text in {"0", "0.0", "negative", "neg", "no"}:
        return "Negative"
    if text in {"equivocal", "borderline"}:
        return "Equivocal"
    if text in {"nan", "none", "", "unknown", "not available"}:
        return None
    return str(value).strip()


def _tnbc(er: pd.Series, pr: pd.Series, her2: pd.Series) -> pd.Series:
    return (er == "Negative") & (pr == "Negative") & (her2 == "Negative")


def _empty_frame(columns: list[str]) -> pd.DataFrame:
    return pd.DataFrame(columns=columns)


def load_tcga() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    clinical = pd.read_csv(CLINICAL / "01_tcga_clinical.csv")
    sample_id = clinical.get("sample_id", clinical.get("sampleID")).astype(str)
    er = clinical.get("ER_Status_nature2012", pd.Series(index=clinical.index, dtype=object)).map(_norm_status)
    pr = clinical.get("PR_Status_nature2012", pd.Series(index=clinical.index, dtype=object)).map(_norm_status)
    her2 = clinical.get("HER2_Final_Status_nature2012", pd.Series(index=clinical.index, dtype=object)).map(_norm_status)

    samples = pd.DataFrame(
        {
            "sample_id": sample_id,
            "cohort": "TCGA-BRCA",
            "patient_id": clinical.get("_PATIENT", sample_id.str[:12]),
            "age_at_dx": pd.to_numeric(
                clinical.get(
                    "age_at_initial_pathologic_diagnosis",
                    clinical.get("Age_at_Initial_Pathologic_Diagnosis_nature2012"),
                ),
                errors="coerce",
            ),
            "sex": clinical.get("Gender_nature2012", clinical.get("gender")),
            "race": clinical.get("race"),
            "ethnicity": clinical.get("ethnicity"),
            "stage_tnm": clinical.get("Converted_Stage_nature2012", clinical.get("pathologic_stage")),
            "grade": clinical.get("neoplasm_histologic_grade"),
            "er_status": er,
            "pr_status": pr,
            "her2_status": her2,
            "intrinsic_subtype_pam50_published": clinical.get("pam50_subtype", clinical.get("PAM50Call_RNAseq")),
            "tnbc_flag": _tnbc(er, pr, her2),
            "primary_site": clinical.get("_primary_site", "Breast"),
            "treatment_summary": None,
        }
    )

    os_days = pd.to_numeric(clinical.get("time_to_event"), errors="coerce") * 30.4375
    if os_days.isna().all() and "OS_Time_nature2012" in clinical:
        os_days = pd.to_numeric(clinical["OS_Time_nature2012"], errors="coerce")
    survival = pd.DataFrame(
        {
            "sample_id": sample_id,
            "os_days": os_days,
            "os_event": pd.to_numeric(clinical.get("event_status", clinical.get("OS_event_nature2012")), errors="coerce"),
            "dfs_days": pd.Series(pd.NA, index=clinical.index, dtype="Float64"),
            "dfs_event": pd.Series(pd.NA, index=clinical.index, dtype="Float64"),
            "last_followup_days": os_days,
        }
    )

    qc = pd.DataFrame(
        {
            "sample_id": sample_id,
            "expression_qc_flag": False,
            "clinical_qc_flag": survival["os_days"].notna() & survival["os_event"].notna(),
            "retained_in_analysis_flag": survival["os_days"].notna() & survival["os_event"].notna(),
            "qc_note": "TCGA clinical harmonized; expression CSV lacks explicit sample_id column in current repo.",
        }
    )
    return samples[SAMPLE_COLUMNS], survival[SURVIVAL_COLUMNS], qc


def load_gse96058() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    clinical = pd.read_csv(CLINICAL / "01_gse96058_clinical.csv")
    er = clinical["er_status"].map(_norm_status)
    pr = clinical["pgr_status"].map(_norm_status)
    her2 = clinical["her2_status"].map(_norm_status)
    samples = pd.DataFrame(
        {
            "sample_id": clinical["sample_id"].astype(str),
            "cohort": "GSE96058",
            "patient_id": clinical["sample_id"].astype(str),
            "age_at_dx": pd.to_numeric(clinical["age_at_diagnosis"], errors="coerce"),
            "sex": "Female",
            "race": None,
            "ethnicity": None,
            "stage_tnm": None,
            "grade": clinical.get("nhg"),
            "er_status": er,
            "pr_status": pr,
            "her2_status": her2,
            "intrinsic_subtype_pam50_published": clinical.get("pam50_subtype"),
            "tnbc_flag": _tnbc(er, pr, her2),
            "primary_site": "Breast",
            "treatment_summary": clinical[["endocrine_treated", "chemo_treated"]].astype(str).agg("; ".join, axis=1),
        }
    )
    survival = pd.DataFrame(
        {
            "sample_id": clinical["sample_id"].astype(str),
            "os_days": pd.to_numeric(clinical["overall_survival_days"], errors="coerce"),
            "os_event": pd.to_numeric(clinical["overall_survival_event"], errors="coerce"),
            "dfs_days": pd.Series(pd.NA, index=clinical.index, dtype="Float64"),
            "dfs_event": pd.Series(pd.NA, index=clinical.index, dtype="Float64"),
            "last_followup_days": pd.to_numeric(clinical["overall_survival_days"], errors="coerce"),
        }
    )
    qc = pd.DataFrame(
        {
            "sample_id": clinical["sample_id"].astype(str),
            "expression_qc_flag": True,
            "clinical_qc_flag": survival["os_days"].notna() & survival["os_event"].notna(),
            "retained_in_analysis_flag": survival["os_days"].notna() & survival["os_event"].notna(),
            "qc_note": "SCAN-B expression and clinical metadata available.",
        }
    )
    return samples[SAMPLE_COLUMNS], survival[SURVIVAL_COLUMNS], qc


def load_metabric() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    clinical = pd.read_parquet(PROCESSED / "03_metabric_clinical.parquet")
    samples = pd.DataFrame(
        {
            "sample_id": clinical["sample_id"].astype(str),
            "cohort": "METABRIC",
            "patient_id": clinical["patient_id"].astype(str),
            "age_at_dx": pd.to_numeric(clinical["age_at_dx"], errors="coerce"),
            "sex": clinical.get("SEX"),
            "race": None,
            "ethnicity": None,
            "stage_tnm": clinical.get("TUMOR_STAGE"),
            "grade": clinical.get("GRADE"),
            "er_status": clinical.get("er_status_norm"),
            "pr_status": clinical.get("pr_status_norm"),
            "her2_status": clinical.get("her2_status_norm"),
            "intrinsic_subtype_pam50_published": clinical.get("CLAUDIN_SUBTYPE"),
            "tnbc_flag": clinical.get("tnbc_flag"),
            "primary_site": "Breast",
            "treatment_summary": clinical[["CHEMOTHERAPY", "HORMONE_THERAPY", "RADIO_THERAPY"]].astype(str).agg("; ".join, axis=1),
        }
    )
    survival = pd.DataFrame(
        {
            "sample_id": clinical["sample_id"].astype(str),
            "os_days": pd.to_numeric(clinical["os_days"], errors="coerce"),
            "os_event": pd.to_numeric(clinical["os_event"], errors="coerce"),
            "dfs_days": pd.to_numeric(clinical["dfs_days"], errors="coerce"),
            "dfs_event": pd.to_numeric(clinical["dfs_event"], errors="coerce"),
            "last_followup_days": pd.to_numeric(clinical["os_days"], errors="coerce"),
        }
    )
    qc = pd.DataFrame(
        {
            "sample_id": clinical["sample_id"].astype(str),
            "expression_qc_flag": True,
            "clinical_qc_flag": survival["os_days"].notna() & survival["os_event"].notna(),
            "retained_in_analysis_flag": survival["os_days"].notna() & survival["os_event"].notna(),
            "qc_note": "METABRIC expression z-scores and clinical metadata available.",
        }
    )
    return samples[SAMPLE_COLUMNS], survival[SURVIVAL_COLUMNS], qc


def load_gse20685() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    clinical = pd.read_parquet(PROCESSED / "04_gse20685_clinical.parquet")
    samples = pd.DataFrame(
        {
            "sample_id": clinical["sample_id"].astype(str),
            "cohort": "GSE20685",
            "patient_id": clinical["patient_id"].astype(str),
            "age_at_dx": pd.to_numeric(clinical["age_at_diagnosis"], errors="coerce"),
            "sex": clinical.get("gender"),
            "race": None,
            "ethnicity": None,
            "stage_tnm": clinical.get("stage_tnm"),
            "grade": None,
            "er_status": None,
            "pr_status": None,
            "her2_status": None,
            "intrinsic_subtype_pam50_published": clinical.get("subtype"),
            "tnbc_flag": pd.Series(pd.NA, index=clinical.index, dtype="boolean"),
            "primary_site": "Breast",
            "treatment_summary": clinical[["adjuvant_chemotherapy", "neoadjuvant_chemotherapy"]].astype(str).agg("; ".join, axis=1),
        }
    )
    survival = pd.DataFrame(
        {
            "sample_id": clinical["sample_id"].astype(str),
            "os_days": pd.to_numeric(clinical["os_days"], errors="coerce"),
            "os_event": pd.to_numeric(clinical["os_event"], errors="coerce"),
            "dfs_days": pd.to_numeric(clinical.get("time_to_relapse_years"), errors="coerce") * 365.25,
            "dfs_event": pd.to_numeric(clinical.get("event_metastasis"), errors="coerce"),
            "last_followup_days": pd.to_numeric(clinical["os_days"], errors="coerce"),
        }
    )
    qc = pd.DataFrame(
        {
            "sample_id": clinical["sample_id"].astype(str),
            "expression_qc_flag": True,
            "clinical_qc_flag": survival["os_days"].notna() & survival["os_event"].notna(),
            "retained_in_analysis_flag": survival["os_days"].notna() & survival["os_event"].notna(),
            "qc_note": "GSE20685 has expression and OS metadata; receptor status unavailable, TNBC excluded.",
        }
    )
    return samples[SAMPLE_COLUMNS], survival[SURVIVAL_COLUMNS], qc


def _quote_ident(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def _expression_view_sql(view_name: str, path: Path, cohort: str, id_columns: list[str], sample_columns: list[str]) -> str:
    id_select = ", ".join(_quote_ident(c) for c in id_columns)
    sample_list = ", ".join(_quote_ident(c) for c in sample_columns)
    gene_expr = _quote_ident("gene_symbol") if "gene_symbol" in id_columns else _quote_ident("Hugo_Symbol")
    return f"""
CREATE OR REPLACE VIEW {view_name} AS
SELECT
  sample_id,
  '{cohort}' AS cohort,
  {gene_expr} AS gene_symbol_hgnc,
  CAST(NULL AS VARCHAR) AS ensembl_id,
  CAST(log2_normalized_expression AS DOUBLE) AS log2_normalized_expression
FROM (
  SELECT {id_select}, sample_id, log2_normalized_expression
  FROM read_parquet('{path.as_posix()}')
  UNPIVOT (log2_normalized_expression FOR sample_id IN ({sample_list}))
);
"""


def build_database() -> dict[str, object]:
    cohort_loaders = [load_tcga, load_gse96058, load_metabric, load_gse20685]
    samples_list: list[pd.DataFrame] = []
    survival_list: list[pd.DataFrame] = []
    qc_list: list[pd.DataFrame] = []
    for loader in cohort_loaders:
        s, surv, qc = loader()
        samples_list.append(s)
        survival_list.append(surv)
        qc_list.append(qc)

    samples = pd.concat(samples_list, ignore_index=True)
    survival = pd.concat(survival_list, ignore_index=True)
    qc = pd.concat(qc_list, ignore_index=True)
    for col in [
        "sample_id",
        "cohort",
        "patient_id",
        "sex",
        "race",
        "ethnicity",
        "stage_tnm",
        "grade",
        "er_status",
        "pr_status",
        "her2_status",
        "intrinsic_subtype_pam50_published",
        "primary_site",
        "treatment_summary",
    ]:
        samples[col] = samples[col].astype("string")

    con = duckdb.connect(str(DB_PATH))
    con.execute("DROP TABLE IF EXISTS samples")
    con.execute("DROP TABLE IF EXISTS survival")
    con.execute("DROP TABLE IF EXISTS qc")
    con.register("samples_df", samples)
    con.register("survival_df", survival)
    con.register("qc_df", qc)
    con.execute("CREATE TABLE samples AS SELECT * FROM samples_df")
    con.execute("CREATE TABLE survival AS SELECT * FROM survival_df")
    con.execute("CREATE TABLE qc AS SELECT * FROM qc_df")

    for view in [
        "gse96058_expression_long",
        "metabric_expression_long",
        "gse20685_expression_long",
        "expression_long",
    ]:
        con.execute(f"DROP VIEW IF EXISTS {view}")

    gse_cols = pd.read_parquet(PROCESSED / "02_gse96058_expression.parquet").columns.tolist()
    met_cols = pd.read_parquet(PROCESSED / "03_metabric_expression.parquet").columns.tolist()
    gse20685_cols = pd.read_parquet(PROCESSED / "04_gse20685_expression.parquet").columns.tolist()
    con.execute(
        _expression_view_sql(
            "gse96058_expression_long",
            PROCESSED / "02_gse96058_expression.parquet",
            "GSE96058",
            ["gene_symbol"],
            [c for c in gse_cols if c != "gene_symbol"],
        )
    )
    con.execute(
        _expression_view_sql(
            "metabric_expression_long",
            PROCESSED / "03_metabric_expression.parquet",
            "METABRIC",
            ["gene_symbol", "entrez_gene_id"],
            [c for c in met_cols if c not in {"gene_symbol", "entrez_gene_id"}],
        )
    )
    con.execute(
        _expression_view_sql(
            "gse20685_expression_long",
            PROCESSED / "04_gse20685_expression.parquet",
            "GSE20685",
            ["gene_symbol"],
            [c for c in gse20685_cols if c != "gene_symbol"],
        )
    )
    con.execute(
        """
CREATE OR REPLACE VIEW expression_long AS
SELECT * FROM gse96058_expression_long
UNION ALL
SELECT * FROM metabric_expression_long
UNION ALL
SELECT * FROM gse20685_expression_long;
"""
    )

    cohort_counts = con.execute(
        """
        SELECT s.cohort,
               COUNT(*) AS n_samples,
               SUM(CASE WHEN q.retained_in_analysis_flag THEN 1 ELSE 0 END) AS n_retained,
               SUM(CASE WHEN s.tnbc_flag THEN 1 ELSE 0 END) AS n_tnbc,
               SUM(CASE WHEN v.os_event = 1 THEN 1 ELSE 0 END) AS n_os_events,
               median(v.os_days) AS median_os_days
        FROM samples s
        LEFT JOIN survival v USING (sample_id)
        LEFT JOIN qc q USING (sample_id)
        GROUP BY s.cohort
        ORDER BY s.cohort
        """
    ).fetchdf()
    con.close()

    return {
        "samples_rows": int(len(samples)),
        "survival_rows": int(len(survival)),
        "qc_rows": int(len(qc)),
        "cohort_counts": cohort_counts,
        "db_path": str(DB_PATH),
    }


def write_report(summary: dict[str, object]) -> None:
    cohort_counts: pd.DataFrame = summary["cohort_counts"]  # type: ignore[assignment]
    lines = [
        "# Harmonization QC",
        "",
        "Timestamp: generated by `scripts/harmonize_cohorts.py`.",
        "",
        f"DuckDB: `{summary['db_path']}`",
        "",
        "## Tables and Views",
        "",
        "- `samples`: harmonized sample-level clinical metadata.",
        "- `survival`: harmonized OS/DFS time-to-event fields.",
        "- `qc`: sample-level expression and clinical QC flags.",
        "- `expression_long`: DuckDB view over external expression parquet files; TCGA expression is excluded because the current repo CSV lacks explicit sample IDs.",
        "",
        "## Cohort Summary",
        "",
        cohort_counts.to_markdown(index=False),
        "",
        "## Normalization Notes",
        "",
        "- TCGA: existing repo-normalized expression file has no explicit sample ID column; TCGA pathway features remain available in `data/processed/02_tcga_feature_matrix.csv`.",
        "- GSE96058: GEO transformed RNA-seq expression, gene rows by sample columns.",
        "- METABRIC: cBioPortal z-scores relative to diploid samples; not directly commensurate with RNA-seq expression scale.",
        "- GSE20685: Affymetrix GPL570 processed expression, probe-level values collapsed to gene symbols by mean.",
        "",
        "## Analysis Flags",
        "",
        "- GSE20685 receptor status was unavailable from parsed GEO metadata, so TNBC-specific primary endpoint analyses must exclude it unless additional metadata is obtained.",
        "- Cross-cohort modeling must avoid fitting normalization parameters on external validation cohorts.",
        "",
    ]
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "harmonization_qc.md").write_text("\n".join(lines))
    (RESULTS / "harmonization_summary.json").write_text(
        json.dumps(
            {
                "db_path": summary["db_path"],
                "samples_rows": summary["samples_rows"],
                "survival_rows": summary["survival_rows"],
                "qc_rows": summary["qc_rows"],
                "cohort_counts": cohort_counts.to_dict(orient="records"),
            },
            indent=2,
        )
        + "\n"
    )


def main() -> None:
    summary = build_database()
    write_report(summary)
    print(json.dumps({k: v for k, v in summary.items() if k != "cohort_counts"}, indent=2))
    print(summary["cohort_counts"])


if __name__ == "__main__":
    main()
