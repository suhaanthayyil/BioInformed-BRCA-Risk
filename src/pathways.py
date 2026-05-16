"""MSigDB pathway parsing and rank-based pathway scoring."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import numpy as np
import pandas as pd


SEVEN_PATHWAY_COMPONENTS = {
    "Immune": [
        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "HALLMARK_INTERFERON_ALPHA_RESPONSE",
        "HALLMARK_INFLAMMATORY_RESPONSE",
        "HALLMARK_ALLOGRAFT_REJECTION",
        "HALLMARK_COMPLEMENT",
    ],
    "Proliferation": [
        "HALLMARK_E2F_TARGETS",
        "HALLMARK_G2M_CHECKPOINT",
        "HALLMARK_MYC_TARGETS_V1",
        "HALLMARK_MYC_TARGETS_V2",
        "HALLMARK_MITOTIC_SPINDLE",
    ],
    "DNA_Repair": [
        "HALLMARK_DNA_REPAIR",
        "REACTOME_HDR_THROUGH_HOMOLOGOUS_RECOMBINATION_HRR",
        "REACTOME_NUCLEOTIDE_EXCISION_REPAIR",
        "REACTOME_MISMATCH_REPAIR",
        "KEGG_MEDICUS_REFERENCE_HOMOLOGOUS_RECOMBINATION",
        "KEGG_MEDICUS_REFERENCE_MISMATCH_REPAIR",
    ],
    "Metabolism": [
        "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
        "HALLMARK_GLYCOLYSIS",
        "HALLMARK_FATTY_ACID_METABOLISM",
        "KEGG_MEDICUS_REFERENCE_GLYCOLYSIS",
    ],
    "Stromal_EMT": [
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
        "HALLMARK_ANGIOGENESIS",
        "HALLMARK_HEDGEHOG_SIGNALING",
    ],
    "Apoptosis_Stress": [
        "HALLMARK_APOPTOSIS",
        "HALLMARK_P53_PATHWAY",
        "HALLMARK_HYPOXIA",
    ],
    "Hormone": [
        "HALLMARK_ESTROGEN_RESPONSE_EARLY",
        "HALLMARK_ESTROGEN_RESPONSE_LATE",
        "HALLMARK_ANDROGEN_RESPONSE",
    ],
}

REACTOME_KEYWORDS = [
    "IMMUNE",
    "INTERFERON",
    "CYTOKINE",
    "T_CELL",
    "B_CELL",
    "CELL_CYCLE",
    "MITOTIC",
    "DNA_REPAIR",
    "HOMOLOGOUS",
    "NUCLEOTIDE_EXCISION",
    "MISMATCH_REPAIR",
    "METABOLISM",
    "GLYCOLYSIS",
    "FATTY",
    "OXIDATIVE",
    "APOPTOSIS",
    "TP53",
    "HYPOXIA",
    "ESTROGEN",
    "ANDROGEN",
    "EXTRACELLULAR_MATRIX",
    "COLLAGEN",
    "ANGIOGENESIS",
    "ERBB",
    "PI3K",
    "MAPK",
]

KEGG_KEYWORDS = [
    "REFERENCE",
    "HOMOLOGOUS_RECOMBINATION",
    "MISMATCH_REPAIR",
    "GLYCOLYSIS",
    "ESTROGEN",
    "P53",
    "CELL_CYCLE",
    "APOPTOSIS",
    "PI3K",
    "MAPK",
]


def clean_symbol(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip().upper()


def read_gmt(path: str | Path) -> dict[str, list[str]]:
    """Read a GMT file into ``name -> genes``."""

    gene_sets: dict[str, list[str]] = {}
    with Path(path).open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name = parts[0]
            genes = [clean_symbol(g) for g in parts[2:] if clean_symbol(g)]
            gene_sets[name] = sorted(set(genes))
    return gene_sets


def select_by_keywords(gene_sets: dict[str, list[str]], keywords: Iterable[str], limit: int) -> dict[str, list[str]]:
    """Select a deterministic keyword-curated subset of a GMT dictionary."""

    keys = []
    upper_keywords = [k.upper() for k in keywords]
    for name in sorted(gene_sets):
        upper = name.upper()
        if any(k in upper for k in upper_keywords):
            keys.append(name)
    return {k: gene_sets[k] for k in keys[:limit]}


def merge_gene_sets(*items: dict[str, list[str]]) -> dict[str, list[str]]:
    merged: dict[str, list[str]] = {}
    for item in items:
        merged.update(item)
    return merged


def expression_wide_to_samples_by_gene(expression: pd.DataFrame, gene_col: str = "gene_symbol") -> pd.DataFrame:
    """Convert genes-by-samples expression data to samples-by-genes."""

    if gene_col not in expression.columns:
        raise ValueError(f"Missing gene column: {gene_col}")
    df = expression.copy()
    df[gene_col] = df[gene_col].map(clean_symbol)
    df = df[df[gene_col] != ""]
    id_cols = {gene_col, "entrez_gene_id", "Entrez_Gene_Id", "ENTREZ_GENE_ID"}
    sample_cols = [c for c in df.columns if c not in id_cols]
    numeric = df[[gene_col, *sample_cols]].set_index(gene_col).apply(pd.to_numeric, errors="coerce")
    numeric = numeric.groupby(level=0).mean()
    out = numeric.T.astype("float32")
    out.index.name = "sample_id"
    return out


def zscore_by_gene(samples_by_gene: pd.DataFrame) -> pd.DataFrame:
    values = samples_by_gene.apply(pd.to_numeric, errors="coerce").astype("float32")
    mean = values.mean(axis=0, skipna=True)
    std = values.std(axis=0, skipna=True, ddof=0).replace(0, np.nan)
    return ((values - mean) / std).fillna(0.0).astype("float32")


def rank_ssgsea_scores(samples_by_gene: pd.DataFrame, gene_sets: dict[str, list[str]], min_genes: int = 5) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute a rank-based ssGSEA-style score for each set and sample.

    The score is the mean percentile rank of genes in the set within each
    sample. It is deterministic and platform-local, and avoids materializing
    per-gene long tables.
    """

    data = samples_by_gene.copy()
    data.columns = [clean_symbol(c) for c in data.columns]
    data = data.loc[:, data.columns != ""]
    ranks = data.rank(axis=1, pct=True, method="average").astype("float32")

    scores = {}
    coverage_rows = []
    columns = set(ranks.columns)
    for name, genes in sorted(gene_sets.items()):
        present = [g for g in genes if g in columns]
        coverage_rows.append(
            {
                "pathway": name,
                "n_genes_defined": len(genes),
                "n_genes_present": len(present),
                "coverage": len(present) / len(genes) if genes else np.nan,
                "scored": len(present) >= min_genes,
            }
        )
        if len(present) >= min_genes:
            scores[name] = ranks[present].mean(axis=1).astype("float32")
    score_df = pd.DataFrame(scores, index=ranks.index)
    score_df.index.name = "sample_id"
    coverage = pd.DataFrame(coverage_rows)
    return score_df, coverage


def aggregate_seven_pathways(pathway_scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Aggregate pathway-set scores into the seven pre-registered features."""

    z = pathway_scores.apply(pd.to_numeric, errors="coerce")
    z = (z - z.mean(axis=0)) / z.std(axis=0, ddof=0).replace(0, np.nan)
    z = z.fillna(0.0)

    features = pd.DataFrame(index=pathway_scores.index)
    rows = []
    for feature, components in SEVEN_PATHWAY_COMPONENTS.items():
        present = [c for c in components if c in z.columns]
        features[f"Pathway_{feature}"] = z[present].mean(axis=1) if present else np.nan
        rows.append(
            {
                "feature": f"Pathway_{feature}",
                "components_requested": ";".join(components),
                "components_present": ";".join(present),
                "n_components_present": len(present),
            }
        )
    features.index.name = "sample_id"
    return features.astype("float32"), pd.DataFrame(rows)
