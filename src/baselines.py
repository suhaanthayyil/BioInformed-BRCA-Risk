"""Clinical baseline and pathway scoring utilities.

The commercial assays in this project are handled conservatively. Public,
citation-backed components are implemented directly; proprietary or incomplete
formulae are returned as unavailable rather than approximated silently.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd

from src.features import PATHWAY_DEFINITIONS
from src.models import get_classifiers, run_cv_evaluation
from src.survival import cox_cv


RANDOM_SEED = 42

PAM50_SUBTYPE_RISK = {
    "luma": 0.15,
    "luminal a": 0.15,
    "lumb": 0.45,
    "luminal b": 0.45,
    "normal": 0.25,
    "normal-like": 0.25,
    "her2": 0.70,
    "her2-enriched": 0.70,
    "her2 enriched": 0.70,
    "basal": 0.85,
    "basal-like": 0.85,
    "claudin-low": 0.65,
}

ONCOTYPE_GENES = {
    "er": ["ESR1", "PGR", "BCL2", "SCUBE2"],
    "proliferation": ["MKI67", "AURKA", "BIRC5", "CCNB1", "MYBL2"],
    "her2": ["GRB7", "ERBB2"],
    "invasion": ["MMP11", "CTSL2", "CTSV"],
    "other": ["CD68", "GSTM1", "BAG1"],
}

PROLIFERATION_GENES = ["MKI67", "AURKA", "BIRC5", "CCNB1", "MYBL2"]


@dataclass(frozen=True)
class ScoreResult:
    """A vector of baseline scores plus provenance metadata."""

    score: pd.Series
    baseline: str
    method_label: str
    status: str
    note: str


def _clean_gene_symbol(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip().upper()


def _coerce_numeric_frame(df: pd.DataFrame) -> pd.DataFrame:
    return df.apply(pd.to_numeric, errors="coerce")


def expression_wide_to_samples_by_gene(
    expression: pd.DataFrame,
    *,
    gene_col: str = "gene_symbol",
    sample_ids: Iterable[str] | None = None,
    keep_genes: Iterable[str] | None = None,
) -> pd.DataFrame:
    """Convert a genes-by-samples expression table to samples-by-genes.

    Duplicate gene symbols are averaged. When ``sample_ids`` is provided, only
    those columns are retained and ordered. The function keeps memory use small
    by subsetting genes before transposition.
    """

    if gene_col not in expression.columns:
        raise ValueError(f"Expression table is missing required gene column: {gene_col}")

    df = expression.copy()
    df[gene_col] = df[gene_col].map(_clean_gene_symbol)
    df = df[df[gene_col] != ""]
    if keep_genes is not None:
        wanted = {_clean_gene_symbol(g) for g in keep_genes}
        df = df[df[gene_col].isin(wanted)]

    id_cols = {gene_col, "entrez_gene_id", "Entrez_Gene_Id", "ENTREZ_GENE_ID"}
    available_samples = [c for c in df.columns if c not in id_cols]
    if sample_ids is not None:
        sample_ids = [str(s) for s in sample_ids]
        available = [s for s in sample_ids if s in df.columns]
    else:
        available = available_samples

    if not available:
        return pd.DataFrame(index=pd.Index([], name="sample_id"))

    numeric = _coerce_numeric_frame(df[[gene_col, *available]].set_index(gene_col))
    numeric = numeric.groupby(level=0).mean()
    out = numeric.T
    out.index.name = "sample_id"
    return out


def samples_by_gene_from_tcga(
    expression_path: str,
    feature_matrix_path: str,
    *,
    keep_genes: Iterable[str] | None = None,
) -> pd.DataFrame:
    """Load the TCGA expression matrix and recover sample IDs from v1 features.

    The current TCGA expression CSV has no sample-id column, but its row order
    exactly reproduces the committed v1 pathway feature matrix. The recovery is
    therefore explicit and auditable instead of inferred silently.
    """

    expression = pd.read_csv(expression_path)
    feature_matrix = pd.read_csv(feature_matrix_path)
    if len(expression) != len(feature_matrix):
        raise ValueError("TCGA expression row count does not match feature matrix row count")
    if keep_genes is not None:
        keep = [_clean_gene_symbol(g) for g in keep_genes]
        rename = {c: _clean_gene_symbol(c) for c in expression.columns}
        expression = expression.rename(columns=rename)
        keep_existing = [g for g in keep if g in expression.columns]
        expression = expression[keep_existing]
    expression.index = feature_matrix["sample_id"].astype(str)
    expression.index.name = "sample_id"
    return _coerce_numeric_frame(expression)


def zscore_by_gene(samples_by_gene: pd.DataFrame) -> pd.DataFrame:
    """Z-score each gene across samples with stable zero-variance handling."""

    values = _coerce_numeric_frame(samples_by_gene)
    mean = values.mean(axis=0, skipna=True)
    std = values.std(axis=0, ddof=0, skipna=True).replace(0, np.nan)
    z = (values - mean) / std
    return z.fillna(0.0)


def percentile_0_100(series: pd.Series) -> pd.Series:
    """Map a score vector to cohort-relative 0-100 percentiles."""

    valid = series.dropna()
    out = pd.Series(np.nan, index=series.index, dtype=float)
    if valid.empty:
        return out
    out.loc[valid.index] = valid.rank(method="average", pct=True) * 100.0
    return out.clip(0, 100)


def _mean_available(z: pd.DataFrame, genes: Iterable[str], min_genes: int = 1) -> tuple[pd.Series, list[str]]:
    wanted = [_clean_gene_symbol(g) for g in genes]
    present = [g for g in wanted if g in z.columns]
    if len(present) < min_genes:
        return pd.Series(np.nan, index=z.index, dtype=float), present
    return z[present].mean(axis=1), present


def oncotype_dx_surrogate(samples_by_gene: pd.DataFrame) -> ScoreResult:
    """Compute a public Oncotype DX 21-gene recurrence-score surrogate.

    This uses the published Paik-style recurrence score group structure and
    weights. Because commercial reference-gene normalization is not available
    for public cohorts, the input is z-scored within cohort and the final vector
    is reported both as RSU and as a cohort-relative 0-100 surrogate.
    """

    z = zscore_by_gene(samples_by_gene)
    er, er_present = _mean_available(z, ONCOTYPE_GENES["er"], min_genes=2)
    prolif, prolif_present = _mean_available(z, ONCOTYPE_GENES["proliferation"], min_genes=3)
    her2, her2_present = _mean_available(z, ONCOTYPE_GENES["her2"], min_genes=1)
    invasion, invasion_present = _mean_available(z, ONCOTYPE_GENES["invasion"], min_genes=1)

    cd68 = z["CD68"] if "CD68" in z else pd.Series(0.0, index=z.index)
    gstm1 = z["GSTM1"] if "GSTM1" in z else pd.Series(0.0, index=z.index)
    bag1 = z["BAG1"] if "BAG1" in z else pd.Series(0.0, index=z.index)

    rsu = (
        0.47 * her2
        - 0.34 * er
        + 1.04 * prolif
        + 0.10 * invasion
        + 0.05 * cd68
        - 0.08 * gstm1
        - 0.07 * bag1
    )
    note = (
        "Paik 2004 group-weight surrogate on within-cohort z-scored expression; "
        f"genes present ER={len(er_present)}, proliferation={len(prolif_present)}, "
        f"HER2={len(her2_present)}, invasion={len(invasion_present)}. "
        "Commercial reference normalization is unavailable."
    )
    status = "ok" if rsu.notna().sum() else "unavailable"
    return ScoreResult(rsu, "Oncotype_DX_21_gene_surrogate", "Paik_2004_RSU_surrogate", status, note)


def pam50_ror_surrogate(
    subtype: pd.Series,
    samples_by_gene: pd.DataFrame | None = None,
    tumor_size: pd.Series | None = None,
) -> ScoreResult:
    """Compute a PAM50-ROR-like risk surrogate from public cohort fields.

    Official PAM50 ROR calculation requires platform-specific centroid mapping
    and formulae generally supplied by packages such as ``genefu``. Until that
    path is available, this function uses published intrinsic subtype labels
    plus proliferation/tumor-size modifiers and labels the output as a surrogate.
    """

    subtype_norm = subtype.astype("string").str.strip().str.lower()
    base = subtype_norm.map(PAM50_SUBTYPE_RISK).astype(float)
    risk = base.copy()
    pieces = ["published intrinsic subtype risk mapping"]

    if samples_by_gene is not None and not samples_by_gene.empty:
        z = zscore_by_gene(samples_by_gene)
        proliferation, present = _mean_available(z, PROLIFERATION_GENES, min_genes=2)
        if proliferation.notna().any():
            risk = risk + 0.20 * percentile_0_100(proliferation).div(100.0)
            pieces.append(f"proliferation modifier ({len(present)} genes)")

    if tumor_size is not None:
        size = pd.to_numeric(tumor_size, errors="coerce")
        if size.notna().any():
            risk = risk + 0.10 * percentile_0_100(size).div(100.0)
            pieces.append("tumor-size modifier")

    scaled = percentile_0_100(risk)
    status = "ok" if scaled.notna().sum() else "unavailable"
    note = (
        "PAM50-ROR surrogate, not official ROR-S/ROR-P/ROR-PT; "
        + ", ".join(pieces)
        + ". Official genefu reproduction should replace this if package scoring succeeds."
    )
    return ScoreResult(scaled, "PAM50_ROR_surrogate", "published_subtype_plus_proliferation", status, note)


def mammaprint_status() -> ScoreResult:
    """Return an explicit unavailable marker for MammaPrint when coefficients are absent."""

    note = (
        "Exact MammaPrint scoring was not computed in Python because the public repo "
        "does not contain the validated 70-gene coefficients/centroid. Use genefu::gene70 "
        "if the R package installation succeeds; otherwise report as not reproducible "
        "from public inputs rather than fabricating a score."
    )
    return ScoreResult(pd.Series(dtype=float), "MammaPrint_gene70", "not_scored", "unavailable", note)


def commercial_formula_status(name: str) -> ScoreResult:
    """Return an explicit unavailable marker for proprietary assays."""

    return ScoreResult(
        pd.Series(dtype=float),
        name,
        "not_scored",
        "unavailable",
        f"{name} exact public formula/coefficients were not available in the repository; score not fabricated.",
    )


def compute_ssgsea_scores(expression_df: pd.DataFrame, pathway_dict: dict[str, list[str]] | None = None) -> pd.DataFrame:
    """Compute a lightweight rank-based ssGSEA-style pathway score."""

    if pathway_dict is None:
        pathway_dict = PATHWAY_DEFINITIONS

    data = expression_df.drop(columns=["sample_id"], errors="ignore")
    data = _coerce_numeric_frame(data)
    n_total = data.shape[1]
    scores = pd.DataFrame(index=expression_df.index)

    if n_total == 0:
        return scores

    ranks = data.rank(axis=1, ascending=True)
    for name, genes in pathway_dict.items():
        valid_genes = [g for g in genes if g in ranks.columns]
        if not valid_genes:
            scores[f"Pathway_{name}"] = 0.0
            continue
        n_set = len(valid_genes)
        expected_rank_sum = n_set * (n_total + 1) / 2
        normalizer = np.sqrt(n_total * n_set)
        scores[f"Pathway_{name}"] = (ranks[valid_genes].sum(axis=1) - expected_rank_sum) / normalizer

    return scores


def compare_scoring_methods(
    meanz_scores: pd.DataFrame,
    ssgsea_scores: pd.DataFrame,
    y: pd.Series,
    time: pd.Series | None = None,
    event: pd.Series | None = None,
    classifiers: dict | None = None,
) -> pd.DataFrame:
    """Compare mean-z and ssGSEA pathway scoring under identical CV."""

    if classifiers is None:
        classifiers = get_classifiers()

    meanz_results = run_cv_evaluation(meanz_scores, y, classifiers)
    ssgsea_results = run_cv_evaluation(ssgsea_scores, y, classifiers)

    comparison = pd.DataFrame(
        {
            "Model": meanz_results["Model"],
            "Mean_Z_AUC": meanz_results["AUC"],
            "ssGSEA_AUC": ssgsea_results["AUC"],
            "Delta": meanz_results["AUC"] - ssgsea_results["AUC"],
        }
    )

    if time is not None and event is not None:
        meanz_cindex, _, _ = cox_cv(meanz_scores, time, event)
        ssgsea_cindex, _, _ = cox_cv(ssgsea_scores, time, event)
        comparison = pd.concat(
            [
                comparison,
                pd.DataFrame(
                    [
                        {
                            "Model": "Cox PH",
                            "Mean_Z_AUC": meanz_cindex,
                            "ssGSEA_AUC": ssgsea_cindex,
                            "Delta": meanz_cindex - ssgsea_cindex,
                        }
                    ]
                ),
            ],
            ignore_index=True,
        )

    return comparison
