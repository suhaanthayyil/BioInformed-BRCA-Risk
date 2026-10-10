#!/usr/bin/env python3
"""Parse raw TCGA-BRCA, GSE96058, METABRIC and GSE20685 downloads into analysis inputs."""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.csv as pacsv

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
CLINICAL = ROOT / "data" / "clinical"
PROCESSED = ROOT / "data" / "processed"

TCGA_EXPR = RAW / "tcga" / "HiSeqV2.gz"
TCGA_CLIN = RAW / "tcga" / "BRCA_clinicalMatrix"
TCGA_SURV = RAW / "tcga" / "BRCA_survival.txt"
GSE96058_EXPR = RAW / "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"
GSE96058_SERIES = [
    RAW / "gse96058" / "GSE96058-GPL11154_series_matrix.txt.gz",
    RAW / "gse96058" / "GSE96058-GPL18573_series_matrix.txt.gz",
]
METABRIC_TAR = RAW / "brca_metabric.tar.gz"
METABRIC_EXPR = RAW / "metabric" / "data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt"
METABRIC_PATIENT = RAW / "metabric" / "data_clinical_patient.txt"
METABRIC_SAMPLE = RAW / "metabric" / "data_clinical_sample.txt"
GSE20685_SERIES = RAW / "GSE20685_series_matrix.txt.gz"
GPL570 = RAW / "GPL570.annot.gz"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def receptor(value) -> str | None:
    if pd.isna(value):
        return None
    text = str(value).strip().lower()
    if text in {"positive", "1"}:
        return "Positive"
    if text in {"negative", "0"}:
        return "Negative"
    if text == "equivocal":
        return "Equivocal"
    return None


def tcga_her2(ihc, fish) -> str | None:
    ihc, fish = receptor(ihc), receptor(fish)
    if ihc == "Positive" or fish == "Positive":
        return "Positive"
    if fish == "Negative":
        return "Negative"
    if ihc == "Negative":
        return "Negative"
    return None


def ingest_tcga() -> dict:
    expr = pd.read_csv(TCGA_EXPR, sep="\t", index_col=0)
    expr = expr[~expr.index.duplicated()]
    clin = pd.read_csv(TCGA_CLIN, sep="\t", low_memory=False).rename(columns={"sampleID": "sample_id"})
    surv = pd.read_csv(TCGA_SURV, sep="\t").rename(columns={"sample": "sample_id"})

    samples = pd.DataFrame({"sample_id": [c for c in expr.columns if c[13:15] == "01"]})
    flow = {"source": "UCSC Xena TCGA.BRCA.sampleMap/HiSeqV2 + TCGA-CDR survival", "rnaseq_samples": int(expr.shape[1])}
    flow["primary_tumour_samples"] = len(samples)

    df = samples.merge(clin, on="sample_id", how="left").merge(surv, on="sample_id", how="left")
    df["patient_id"] = df["sample_id"].str[:12]
    df = df.drop_duplicates("patient_id")
    flow["unique_patients"] = len(df)

    male = df["gender"].ne("FEMALE")
    flow["excluded_male_or_sex_missing"] = int(male.sum())
    df = df[~male]
    redacted = df["Redaction"].eq("Redacted")
    flow["excluded_redacted"] = int(redacted.sum())
    df = df[~redacted]
    no_os = df["OS"].isna() | df["OS.time"].isna() | (pd.to_numeric(df["OS.time"], errors="coerce") <= 0)
    flow["excluded_missing_or_zero_os"] = int(no_os.sum())
    df = df[~no_os]
    flow["included"] = len(df)
    flow["events"] = int(df["OS"].sum())

    out = pd.DataFrame(
        {
            "sample_id": df["sample_id"],
            "patient_id": df["patient_id"],
            "age_at_dx": pd.to_numeric(df["age_at_initial_pathologic_diagnosis"], errors="coerce"),
            "sex": "Female",
            "stage_tnm": df["pathologic_stage"],
            "er_status": df["breast_carcinoma_estrogen_receptor_status"].map(receptor),
            "pr_status": df["breast_carcinoma_progesterone_receptor_status"].map(receptor),
            "her2_status": [
                tcga_her2(a, b)
                for a, b in zip(
                    df["lab_proc_her2_neu_immunohistochemistry_receptor_status"],
                    df["lab_procedure_her2_neu_in_situ_hybrid_outcome_type"],
                    strict=True,
                )
            ],
            "pam50_published": df["PAM50Call_RNAseq"],
            "os_days": pd.to_numeric(df["OS.time"], errors="coerce"),
            "os_event": pd.to_numeric(df["OS"], errors="coerce").astype(int),
            "pfi_days": pd.to_numeric(df["PFI.time"], errors="coerce"),
            "pfi_event": pd.to_numeric(df["PFI"], errors="coerce"),
        }
    )
    out.to_csv(CLINICAL / "tcga_brca_clinical.csv", index=False)

    keep = expr[out["sample_id"].tolist()].astype("float32")
    keep.index.name = "gene_symbol"
    keep.reset_index().to_parquet(PROCESSED / "01_tcga_expression.parquet", index=False)
    flow["expression_genes"] = int(keep.shape[0])
    flow["sha256"] = {p.name: sha256(p) for p in [TCGA_EXPR, TCGA_CLIN, TCGA_SURV]}
    return flow


def read_series_table(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt") as fh:
        lines = fh.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("!series_matrix_table_begin"))
    end = next(i for i, line in enumerate(lines) if line.startswith("!series_matrix_table_end"))
    return pd.read_csv(io.StringIO("".join(lines[start + 1 : end])), sep="\t", index_col=0)


def read_series_matrix(path: Path) -> pd.DataFrame:
    titles, accessions, fields, sources, descriptions, platforms = [], [], [], [], [], []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                accessions = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_source_name_ch1"):
                sources = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_description"):
                descriptions = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_platform_id"):
                platforms = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                fields.append([x.strip('"') for x in line.rstrip("\n").split("\t")[1:]])
            elif line.startswith("!series_matrix_table_begin"):
                break
    rows = []
    for i, title in enumerate(titles):
        row = {"sample_id": title, "geo_accession": accessions[i]}
        for key, values in [("source_name_ch1", sources), ("description", descriptions), ("platform_id", platforms)]:
            if values:
                row[key] = values[i]
        for field in fields:
            if ": " in field[i]:
                key, value = field[i].split(": ", 1)
                key = key.strip().replace(" (", "_").replace(")", "").replace(" ", "_").replace("-", "_")
                row[key] = None if value == "NA" else value
        rows.append(row)
    return pd.DataFrame(rows)


def ingest_gse96058() -> dict:
    clin = pd.concat([read_series_matrix(p) for p in GSE96058_SERIES], ignore_index=True)
    flow = {"source": "GEO GSE96058 series matrices (GPL11154, GPL18573)", "geo_samples": len(clin)}
    replicate = clin["sample_id"].str.contains("repl", case=False)
    flow["excluded_technical_replicates"] = int(replicate.sum())
    clin = clin[~replicate].copy()
    flow["unique_patients"] = len(clin)

    clin["os_days"] = pd.to_numeric(clin["overall_survival_days"], errors="coerce")
    clin["os_event"] = pd.to_numeric(clin["overall_survival_event"], errors="coerce")
    no_os = clin["os_days"].isna() | clin["os_event"].isna() | (clin["os_days"] <= 0)
    flow["excluded_missing_or_zero_os"] = int(no_os.sum())
    clin = clin[~no_os]

    table = pacsv.read_csv(GSE96058_EXPR).to_pandas()
    gene_col = table.columns[0]
    table = table.rename(columns={gene_col: "gene_symbol"})
    no_expr = ~clin["sample_id"].isin(table.columns)
    flow["excluded_no_expression"] = int(no_expr.sum())
    clin = clin[~no_expr]
    flow["included"] = len(clin)
    flow["events"] = int(clin["os_event"].sum())

    clin.to_csv(CLINICAL / "gse96058_clinical.csv", index=False)
    cols = ["gene_symbol", *clin["sample_id"].tolist()]
    expr = table[cols].copy()
    expr[cols[1:]] = expr[cols[1:]].astype("float32")
    expr.to_parquet(PROCESSED / "02_gse96058_expression.parquet", index=False)
    flow["expression_genes"] = int(len(expr))
    flow["sha256"] = {p.name: sha256(p) for p in [GSE96058_EXPR, *GSE96058_SERIES]}
    return flow


def ingest_metabric() -> dict:
    if not METABRIC_EXPR.exists():
        with tarfile.open(METABRIC_TAR) as tar:
            member = next(m for m in tar.getmembers() if m.name.endswith(METABRIC_EXPR.name))
            member.name = METABRIC_EXPR.name
            tar.extract(member, METABRIC_EXPR.parent, filter="data")
    expr = pd.read_csv(METABRIC_EXPR, sep="\t")
    expr = expr.rename(columns={"Hugo_Symbol": "gene_symbol", "Entrez_Gene_Id": "entrez_gene_id"})
    samples = [c for c in expr.columns if c not in {"gene_symbol", "entrez_gene_id"}]
    expr[samples] = expr[samples].astype("float32")
    expr.to_parquet(PROCESSED / "03_metabric_expression.parquet", index=False)
    return {"expression_genes": int(len(expr)), "expression_samples": len(samples), "sha256": {METABRIC_TAR.name: sha256(METABRIC_TAR)}}


def read_cbioportal(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", comment="#", dtype=str)


def receptor_series(values: pd.Series) -> pd.Series:
    return values.map(receptor)


def ingest_metabric_clinical() -> dict:
    patient = read_cbioportal(METABRIC_PATIENT)
    sample = read_cbioportal(METABRIC_SAMPLE)
    df = sample.merge(patient, on="PATIENT_ID", how="left")
    df["sample_id"] = df["SAMPLE_ID"].astype(str)
    df["patient_id"] = df["PATIENT_ID"].astype(str)
    df["cohort"] = "METABRIC"
    df["age_at_dx"] = pd.to_numeric(df["AGE_AT_DIAGNOSIS"], errors="coerce")
    df["os_days"] = pd.to_numeric(df["OS_MONTHS"], errors="coerce") * 30.4375
    df["os_event"] = df["OS_STATUS"].astype(str).str.startswith("1").astype(int)
    df["dfs_days"] = pd.to_numeric(df["RFS_MONTHS"], errors="coerce") * 30.4375
    df["dfs_event"] = df["RFS_STATUS"].astype(str).str.startswith("1").astype(int)
    df["er_status_norm"] = receptor_series(df["ER_STATUS"])
    df["pr_status_norm"] = receptor_series(df["PR_STATUS"])
    df["her2_status_norm"] = receptor_series(df["HER2_STATUS"])
    df["tnbc_flag"] = df[["er_status_norm", "pr_status_norm", "her2_status_norm"]].eq("Negative").all(axis=1)
    df.to_parquet(PROCESSED / "03_metabric_clinical.parquet", index=False)
    return {"clinical_rows": len(df), "sha256": {p.name: sha256(p) for p in [METABRIC_PATIENT, METABRIC_SAMPLE]}}


def ingest_gse20685() -> dict:
    clin = read_series_matrix(GSE20685_SERIES).drop(columns=["sample_id"]).rename(columns={"geo_accession": "sample_id"})
    clin["cohort"] = "GSE20685"
    clin["patient_id"] = clin["sample_id"]
    clin["age_at_diagnosis"] = pd.to_numeric(clin["age_at_diagnosis"], errors="coerce")
    clin["os_days"] = pd.to_numeric(clin["follow_up_duration_years"], errors="coerce") * 365.25
    clin["follow_up_duration_years"] = pd.to_numeric(clin["follow_up_duration_years"], errors="coerce")
    clin["regional_relapse"] = pd.to_numeric(clin["regional_relapse"], errors="coerce")
    clin["os_event"] = pd.to_numeric(clin["event_death"], errors="coerce")
    clin["t_stage"] = pd.to_numeric(clin["t_stage"], errors="coerce")
    clin["stage_tnm"] = [
        f"T{t}N{n}M{m}" if pd.notna(t) and pd.notna(n) and pd.notna(m) else None
        for t, n, m in zip(clin["t_stage"], clin["n_stage"], clin["m_stage"], strict=True)
    ]
    clin["tnbc_flag"] = None
    clin.to_parquet(PROCESSED / "04_gse20685_clinical.parquet", index=False)

    expr = read_series_table(GSE20685_SERIES)
    annot = pd.read_csv(GPL570, sep="\t", skiprows=27, usecols=["ID", "Gene symbol"], dtype=str)
    symbol = annot.set_index("ID")["Gene symbol"].reindex(expr.index)
    expr = expr[symbol.notna()].groupby(symbol[symbol.notna()].to_numpy()).mean()
    expr.index.name = "gene_symbol"
    expr.reset_index().to_parquet(PROCESSED / "04_gse20685_expression.parquet", index=False)
    return {"clinical_rows": len(clin), "expression_genes": len(expr), "sha256": {p.name: sha256(p) for p in [GSE20685_SERIES, GPL570]}}


def merge_flow(a: dict, b: dict) -> dict:
    return {**a, **b, "sha256": {**a.get("sha256", {}), **b.get("sha256", {})}}


def main() -> None:
    CLINICAL.mkdir(parents=True, exist_ok=True)
    flow = {
        "TCGA-BRCA": ingest_tcga(),
        "GSE96058": ingest_gse96058(),
        "METABRIC": merge_flow(ingest_metabric(), ingest_metabric_clinical()),
        "GSE20685": ingest_gse20685(),
    }
    (PROCESSED / "ingest_flow.json").write_text(json.dumps(flow, indent=2, default=lambda x: x.item() if isinstance(x, np.generic) else str(x)))
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != "sha256"} for k, v in flow.items()}, indent=2))


if __name__ == "__main__":
    main()
