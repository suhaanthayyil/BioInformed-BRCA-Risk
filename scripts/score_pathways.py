#!/usr/bin/env python3
"""Phase 4 pathway scoring and seven-pathway aggregation."""

from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from src.pathways import (  # noqa: E402
    KEGG_KEYWORDS,
    REACTOME_KEYWORDS,
    SEVEN_PATHWAY_COMPONENTS,
    aggregate_seven_pathways,
    expression_wide_to_samples_by_gene,
    merge_gene_sets,
    rank_ssgsea_scores,
    read_gmt,
    select_by_keywords,
    zscore_by_gene,
)


DATA = REPO_ROOT / "data"
PROCESSED = DATA / "processed"
GENE_SETS = DATA / "raw" / "gene_sets"
RESULTS = REPO_ROOT / "results" / "reports"
LOGS = REPO_ROOT / "logs"
SEED = 42
MSIGDB_VERSION = "v2024.1.Hs"


def timestamp() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def log(message: str) -> None:
    LOGS.mkdir(exist_ok=True)
    with (LOGS / "agent_log.md").open("a") as fh:
        fh.write(f"\n- {timestamp()} {message}\n")


def load_tcga_expression() -> pd.DataFrame:
    expr = pd.read_csv(DATA / "01_tcga_expression_normalized.csv")
    feature_matrix = pd.read_csv(PROCESSED / "02_tcga_feature_matrix.csv")
    if "sample_id" in expr.columns:
        expr = expr.set_index(expr["sample_id"].astype(str)).drop(columns=["sample_id"])
        expr.index.name = "sample_id"
        expr = expr.reindex(feature_matrix["sample_id"].astype(str))
        expr.columns = [str(c).upper() for c in expr.columns]
        return expr.apply(pd.to_numeric, errors="coerce").astype("float32")
    if len(expr) != len(feature_matrix):
        raise ValueError("TCGA expression rows do not match TCGA feature matrix rows")
    expr.index = feature_matrix["sample_id"].astype(str)
    expr.index.name = "sample_id"
    expr.columns = [str(c).upper() for c in expr.columns]
    return expr.astype("float32")


def load_expression(cohort: str) -> pd.DataFrame:
    if cohort == "TCGA-BRCA":
        return load_tcga_expression()
    if cohort == "GSE96058":
        return expression_wide_to_samples_by_gene(pd.read_parquet(PROCESSED / "02_gse96058_expression.parquet"))
    if cohort == "METABRIC":
        return expression_wide_to_samples_by_gene(pd.read_parquet(PROCESSED / "03_metabric_expression.parquet"))
    if cohort == "GSE20685":
        return expression_wide_to_samples_by_gene(pd.read_parquet(PROCESSED / "04_gse20685_expression.parquet"))
    raise ValueError(f"Unknown cohort: {cohort}")


def build_gene_sets() -> tuple[dict[str, list[str]], dict[str, int]]:
    hallmark = read_gmt(GENE_SETS / f"h.all.{MSIGDB_VERSION}.symbols.gmt")
    reactome_all = read_gmt(GENE_SETS / f"c2.cp.reactome.{MSIGDB_VERSION}.symbols.gmt")
    kegg_all = read_gmt(GENE_SETS / f"c2.cp.kegg_medicus.{MSIGDB_VERSION}.symbols.gmt")
    reactome = select_by_keywords(reactome_all, REACTOME_KEYWORDS, limit=300)
    kegg = select_by_keywords(kegg_all, KEGG_KEYWORDS, limit=50)
    mandatory = {}
    for components in SEVEN_PATHWAY_COMPONENTS.values():
        for name in components:
            if name in hallmark:
                mandatory[name] = hallmark[name]
            elif name in reactome_all:
                mandatory[name] = reactome_all[name]
            elif name in kegg_all:
                mandatory[name] = kegg_all[name]
    return merge_gene_sets(hallmark, reactome, kegg, mandatory), {
        "hallmark": len(hallmark),
        "reactome_curated": len(reactome),
        "kegg_medicus_curated": len(kegg),
        "mandatory_aggregation_sets": len(mandatory),
    }


def score_one_cohort(cohort: str, gene_sets: dict[str, list[str]]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    expr = load_expression(cohort)
    expr = zscore_by_gene(expr)
    scores, coverage = rank_ssgsea_scores(expr, gene_sets)
    features, agg_qc = aggregate_seven_pathways(scores)

    scores = scores.reset_index()
    scores.insert(1, "cohort", cohort)
    features = features.reset_index()
    features.insert(1, "cohort", cohort)
    coverage.insert(0, "cohort", cohort)
    agg_qc.insert(0, "cohort", cohort)
    return scores, features, pd.concat([coverage.assign(qc_type="gene_coverage"), agg_qc.assign(qc_type="aggregation")], ignore_index=True)


def write_report(summary: dict, coverage: pd.DataFrame) -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    coverage_part = coverage[coverage["qc_type"].eq("gene_coverage")].copy()
    agg_part = coverage[coverage["qc_type"].eq("aggregation")].copy()
    lines = [
        "# Phase 4 Pathway Scoring QC",
        "",
        f"Generated: {timestamp()}",
        f"MSigDB version: {MSIGDB_VERSION}",
        "",
        "## Gene Set Inventory",
        "",
        pd.DataFrame([summary["gene_set_counts"]]).to_markdown(index=False),
        "",
        "## Cohort Output Shapes",
        "",
        pd.DataFrame(summary["cohort_shapes"]).to_markdown(index=False),
        "",
        "## Seven-Pathway Aggregation Coverage",
        "",
        agg_part[
            ["cohort", "feature", "n_components_present", "components_present"]
        ].to_markdown(index=False),
        "",
        "## Gene Coverage Summary",
        "",
        coverage_part.groupby("cohort")
        .agg(
            pathways_scored=("scored", "sum"),
            median_coverage=("coverage", "median"),
            min_coverage=("coverage", "min"),
        )
        .reset_index()
        .to_markdown(index=False, floatfmt=".3f"),
        "",
        "Scores use a deterministic rank-percentile ssGSEA-style implementation. "
        "The R/GSVA entrypoint remains available for package-based cross-checking.",
        "",
    ]
    (RESULTS / "pathway_scoring_qc.md").write_text("\n".join(lines))


def main() -> None:
    np.random.seed(SEED)
    log("Phase 4 pathway scoring started.")
    gene_sets, counts = build_gene_sets()
    all_scores = []
    all_features = []
    all_coverage = []
    shapes = []

    for cohort in ["TCGA-BRCA", "GSE96058", "METABRIC", "GSE20685"]:
        log(f"Phase 4 scoring cohort {cohort}.")
        scores, features, coverage = score_one_cohort(cohort, gene_sets)
        all_scores.append(scores)
        all_features.append(features)
        all_coverage.append(coverage)
        shapes.append(
            {
                "cohort": cohort,
                "n_samples_pathway_scores": scores.shape[0],
                "n_pathway_sets": scores.shape[1] - 2,
                "n_samples_features": features.shape[0],
                "n_features": features.shape[1] - 2,
            }
        )

    pathway_scores = pd.concat(all_scores, ignore_index=True)
    pathway_features = pd.concat(all_features, ignore_index=True)
    coverage = pd.concat(all_coverage, ignore_index=True)

    pathway_scores.to_parquet(PROCESSED / "pathway_scores_all.parquet", index=False)
    pathway_features.to_parquet(PROCESSED / "04_pathway_features.parquet", index=False)
    coverage.to_csv(RESULTS / "pathway_scoring_coverage.csv", index=False)

    summary = {
        "generated_at": timestamp(),
        "seed": SEED,
        "msigdb_version": MSIGDB_VERSION,
        "gene_set_counts": counts,
        "total_gene_sets_scored": len(gene_sets),
        "cohort_shapes": shapes,
        "score_file": "data/processed/pathway_scores_all.parquet",
        "feature_file": "data/processed/04_pathway_features.parquet",
        "coverage_file": "results/reports/pathway_scoring_coverage.csv",
    }
    (PROCESSED / "pathway_scoring.metadata.json").write_text(json.dumps(summary, indent=2))
    write_report(summary, coverage)
    log("Phase 4 pathway scoring completed.")


if __name__ == "__main__":
    main()
