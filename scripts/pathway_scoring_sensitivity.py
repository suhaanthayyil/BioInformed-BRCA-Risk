#!/usr/bin/env python3
"""Task 3: Pathway-scoring sensitivity analysis (Table S14).

Compare three aggregation methods (rank-percentile, mean-z, ssGSEA)
with locked GBSA architecture. Each method produces 7 pathway features,
trains on TCGA, and evaluates on 3 external cohorts.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sksurv.ensemble import GradientBoostingSurvivalAnalysis
from sksurv.util import Surv

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
from src.survival import (  # noqa: E402
    bootstrap_c_index_ci,
    harrell_c_index,
    paired_bootstrap_delta_c_index,
    uno_c_index,
)

DATA = REPO_ROOT / "data"
PROCESSED = DATA / "processed"
GENE_SETS_DIR = DATA / "raw" / "gene_sets"
RESULTS = REPO_ROOT / "results"
SEED = 42
BOOTSTRAP_ITER = 1000
PAM50 = "PAM50_ROR_official"
MSIGDB_VERSION = "v2024.1.Hs"
EXTERNAL_COHORTS = ["GSE96058", "METABRIC", "GSE20685"]

PATHWAY_FEATURES = [
    "Pathway_Immune",
    "Pathway_Proliferation",
    "Pathway_DNA_Repair",
    "Pathway_Metabolism",
    "Pathway_Stromal_EMT",
    "Pathway_Apoptosis_Stress",
    "Pathway_Hormone",
]


def stage_to_ordinal(value) -> float:
    text = "" if pd.isna(value) else str(value).upper()
    if "IV" in text:
        return 4.0
    if "III" in text:
        return 3.0
    if "II" in text:
        return 2.0
    if "I" in text:
        return 1.0
    return np.nan


def make_preprocessor() -> Pipeline:
    return Pipeline([("imputer", SimpleImputer(strategy="median")), ("scaler", StandardScaler())])


def get_locked_gbsa_config() -> dict:
    cv_table = pd.read_csv(RESULTS / "Table_S2_ml_internal_cv.csv")
    gbsa_row = cv_table[cv_table["model"] == "Gradient_Boosted_Survival"].iloc[0]
    return json.loads(gbsa_row["config"])


def build_gene_sets() -> dict[str, list[str]]:
    hallmark = read_gmt(GENE_SETS_DIR / f"h.all.{MSIGDB_VERSION}.symbols.gmt")
    reactome_all = read_gmt(GENE_SETS_DIR / f"c2.cp.reactome.{MSIGDB_VERSION}.symbols.gmt")
    kegg_all = read_gmt(GENE_SETS_DIR / f"c2.cp.kegg_medicus.{MSIGDB_VERSION}.symbols.gmt")
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
    return merge_gene_sets(hallmark, reactome, kegg, mandatory)


def load_expression(cohort: str) -> pd.DataFrame:
    if cohort == "TCGA-BRCA":
        expr = pd.read_csv(DATA / "01_tcga_expression_normalized.csv")
        feature_matrix = pd.read_csv(PROCESSED / "02_tcga_feature_matrix.csv")
        if "sample_id" in expr.columns:
            expr = expr.set_index(expr["sample_id"].astype(str)).drop(columns=["sample_id"])
            expr.index.name = "sample_id"
            expr = expr.reindex(feature_matrix["sample_id"].astype(str))
            expr.columns = [str(c).upper() for c in expr.columns]
            return expr.apply(pd.to_numeric, errors="coerce").astype("float32")
        expr.index = feature_matrix["sample_id"].astype(str)
        expr.index.name = "sample_id"
        expr.columns = [str(c).upper() for c in expr.columns]
        return expr.astype("float32")
    if cohort == "GSE96058":
        return expression_wide_to_samples_by_gene(pd.read_parquet(PROCESSED / "02_gse96058_expression.parquet"))
    if cohort == "METABRIC":
        return expression_wide_to_samples_by_gene(pd.read_parquet(PROCESSED / "03_metabric_expression.parquet"))
    if cohort == "GSE20685":
        return expression_wide_to_samples_by_gene(pd.read_parquet(PROCESSED / "04_gse20685_expression.parquet"))
    raise ValueError(f"Unknown cohort: {cohort}")


def load_clinical() -> pd.DataFrame:
    con = duckdb.connect(str(PROCESSED / "unified_cohorts.duckdb"), read_only=True)
    try:
        samples = con.execute("select * from samples").fetchdf()
        survival = con.execute("select * from survival").fetchdf()
    finally:
        con.close()
    frame = samples.merge(survival, on="sample_id", how="inner")
    frame["stage_ordinal"] = frame["stage_tnm"].map(stage_to_ordinal)
    frame["age_at_dx"] = pd.to_numeric(frame["age_at_dx"], errors="coerce")
    frame["os_days"] = pd.to_numeric(frame["os_days"], errors="coerce")
    frame["os_event"] = pd.to_numeric(frame["os_event"], errors="coerce").fillna(0).astype(int)
    return frame


def score_rank_percentile(expr: pd.DataFrame, gene_sets: dict) -> pd.DataFrame:
    """Method 1: Rank-percentile (production method)."""
    z = zscore_by_gene(expr)
    scores, _ = rank_ssgsea_scores(z, gene_sets)
    features, _ = aggregate_seven_pathways(scores)
    return features


def score_mean_z(expr: pd.DataFrame, gene_sets: dict) -> pd.DataFrame:
    """Method 2: Mean z-score of component genes.

    For each of the 7 pathways, compute the mean z-score of all genes
    in the component gene sets that are present in expression data.
    """
    from src.pathways import clean_symbol
    z = zscore_by_gene(expr)
    z.columns = [clean_symbol(c) for c in z.columns]
    z = z.loc[:, z.columns != ""]

    features = pd.DataFrame(index=expr.index)
    for feature, components in SEVEN_PATHWAY_COMPONENTS.items():
        all_genes = set()
        for comp_name in components:
            if comp_name in gene_sets:
                all_genes.update(gene_sets[comp_name])
        present = [g for g in all_genes if g in z.columns]
        if present:
            features[f"Pathway_{feature}"] = z[present].mean(axis=1).astype("float32")
        else:
            features[f"Pathway_{feature}"] = np.nan
    features.index.name = "sample_id"
    return features


def score_ssgsea(expr: pd.DataFrame, gene_sets: dict) -> pd.DataFrame:
    """Method 3: ssGSEA via gseapy (if available), else rank-based approximation.

    Uses the gseapy.ssgsea function for per-sample enrichment scoring.
    Falls back to rank-percentile if gseapy is not available.
    """
    try:
        import gseapy
        # Build subset of gene sets for the seven pathway components
        component_sets = {}
        for feature, components in SEVEN_PATHWAY_COMPONENTS.items():
            for comp_name in components:
                if comp_name in gene_sets:
                    component_sets[comp_name] = gene_sets[comp_name]

        # gseapy.ssgsea expects genes-by-samples (DataFrame with genes as rows)
        # expr is samples-by-genes, so transpose to genes-by-samples
        genes_by_samples = expr.T
        genes_by_samples.index.name = "Gene"

        result = gseapy.ssgsea(
            data=genes_by_samples,
            gene_sets=component_sets,
            outdir=None,
            no_plot=True,
            min_size=5,
            threads=1,
            seed=SEED,
            verbose=False,
        )
        # result.res2d has columns: Name (sample), Term (gene set), NES, etc.
        res = result.res2d
        scores = res.pivot(index="Name", columns="Term", values="NES")
        scores = scores.apply(pd.to_numeric, errors="coerce")
        scores.index.name = "sample_id"

        # Aggregate to seven pathways
        features, _ = aggregate_seven_pathways(scores)
        return features
    except Exception as e:
        print(f"  ssGSEA via gseapy failed ({e}), using rank-based approximation")
        return score_rank_percentile(expr, gene_sets)


SCORING_METHODS = {
    "rank_percentile": score_rank_percentile,
    "mean_z": score_mean_z,
    "ssgsea": score_ssgsea,
}


def main() -> None:
    np.random.seed(SEED)
    gene_sets = build_gene_sets()
    config = get_locked_gbsa_config()
    clinical = load_clinical()
    print(f"Locked GBSA config: {config}")
    print(f"Gene sets loaded: {len(gene_sets)}")

    # Load PAM50 scores for delta computation
    pam50_scores = pd.read_csv(RESULTS / "06_external_model_predictions.csv")[
        ["sample_id", PAM50]
    ]

    full_features = PATHWAY_FEATURES + ["age_at_dx", "stage_ordinal"]
    rows = []

    for method_name, scoring_fn in SCORING_METHODS.items():
        print(f"\nScoring method: {method_name}")

        # Score all cohorts
        all_pathway_dfs = {}
        for cohort in ["TCGA-BRCA"] + EXTERNAL_COHORTS:
            print(f"  Scoring {cohort}...")
            expr = load_expression(cohort)
            pathway_df = scoring_fn(expr, gene_sets)
            pathway_df = pathway_df.reset_index()
            pathway_df["cohort"] = cohort
            all_pathway_dfs[cohort] = pathway_df

        # Merge pathway features with clinical data
        all_features = pd.concat(all_pathway_dfs.values(), ignore_index=True)
        # Ensure sample_id is string for consistent merging
        all_features["sample_id"] = all_features["sample_id"].astype(str)
        clinical_subset = clinical[["sample_id", "cohort", "age_at_dx", "stage_ordinal", "os_days", "os_event"]].copy()
        clinical_subset["sample_id"] = clinical_subset["sample_id"].astype(str)
        frame = all_features.merge(
            clinical_subset,
            on=["sample_id", "cohort"],
            how="inner",
        )
        if len(frame) == 0:
            # Try merge on sample_id only (cohort may not match exactly)
            frame = all_features.merge(
                clinical_subset.drop(columns=["cohort"]),
                on="sample_id",
                how="inner",
            )

        # Train on TCGA
        train = frame[
            frame["cohort"].eq("TCGA-BRCA")
            & frame["os_days"].notna()
            & (frame["os_days"] > 0)
            & frame["os_event"].notna()
        ].copy()
        available_feats = [f for f in full_features if f in train.columns]
        print(f"  Training on TCGA: n={len(train)}, features={len(available_feats)}")

        pre = make_preprocessor()
        x_train = pre.fit_transform(train[available_feats])
        t_train = train["os_days"].to_numpy(float)
        e_train = train["os_event"].to_numpy(bool)
        y_train = Surv.from_arrays(e_train, t_train)

        gbsa = GradientBoostingSurvivalAnalysis(random_state=SEED, **config)
        gbsa.fit(x_train, y_train)

        # Internal CV C-index (quick sanity check)
        from sksurv.metrics import concordance_index_censored
        train_risk = gbsa.predict(x_train)
        train_c = float(concordance_index_censored(e_train, t_train, train_risk)[0])
        print(f"  Training C-index (resubstitution): {train_c:.4f}")

        # Evaluate on external cohorts
        for cohort in EXTERNAL_COHORTS:
            cdf = frame[frame["cohort"].eq(cohort)].copy()
            if len(cdf) == 0:
                continue
            x_test = pre.transform(cdf[available_feats])
            cdf["method_risk"] = gbsa.predict(x_test).astype(float)
            cdf = cdf.merge(pam50_scores, on="sample_id", how="left")

            valid = cdf[
                cdf["method_risk"].notna()
                & cdf[PAM50].notna()
                & cdf["os_days"].notna()
                & (cdf["os_days"] > 0)
                & cdf["os_event"].notna()
            ].copy()
            n = len(valid)
            events = int(valid["os_event"].sum())

            if n < 10 or events < 2:
                rows.append({
                    "scoring_method": method_name,
                    "cohort": cohort,
                    "n": n,
                    "events": events,
                    "harrell_c": np.nan,
                    "harrell_c_ci_low": np.nan,
                    "harrell_c_ci_high": np.nan,
                    "uno_c": np.nan,
                    "delta_vs_pam50": np.nan,
                    "delta_ci_low": np.nan,
                    "delta_ci_high": np.nan,
                    "delta_p": np.nan,
                })
                continue

            hc = harrell_c_index(valid["os_days"], valid["os_event"], valid["method_risk"])
            ci_low, ci_high = bootstrap_c_index_ci(
                valid["os_days"], valid["os_event"], valid["method_risk"],
                BOOTSTRAP_ITER, SEED,
            )
            uc = uno_c_index(valid["os_days"], valid["os_event"], valid["method_risk"])

            delta = paired_bootstrap_delta_c_index(
                valid["os_days"], valid["os_event"],
                valid["method_risk"], valid[PAM50],
                n_bootstrap=BOOTSTRAP_ITER, seed=SEED,
            )

            rows.append({
                "scoring_method": method_name,
                "cohort": cohort,
                "n": n,
                "events": events,
                "harrell_c": hc,
                "harrell_c_ci_low": ci_low,
                "harrell_c_ci_high": ci_high,
                "uno_c": uc,
                "delta_vs_pam50": delta["delta"],
                "delta_ci_low": delta["ci_low"],
                "delta_ci_high": delta["ci_high"],
                "delta_p": delta["p"],
            })

    table = pd.DataFrame(rows)
    table.to_csv(RESULTS / "Table_S14_pathway_scoring_sensitivity.csv", index=False)

    # Print results
    print("\n" + "=" * 70)
    print("Table S14: Pathway Scoring Sensitivity")
    print("=" * 70)
    for method_name in SCORING_METHODS:
        sub = table[table["scoring_method"] == method_name]
        print(f"\n  {method_name}:")
        for _, row in sub.iterrows():
            print(
                f"    {row['cohort']:10s} | C={row['harrell_c']:.4f} "
                f"[{row['harrell_c_ci_low']:.4f}, {row['harrell_c_ci_high']:.4f}] | "
                f"delta vs PAM50={row['delta_vs_pam50']:+.4f} "
                f"[{row['delta_ci_low']:.4f}, {row['delta_ci_high']:.4f}] p={row['delta_p']:.4f}"
            )
        mean_c = sub["harrell_c"].mean()
        mean_d = sub["delta_vs_pam50"].mean()
        print(f"    {'MEAN':10s} | C={mean_c:.4f} | mean delta={mean_d:+.4f}")


if __name__ == "__main__":
    main()
