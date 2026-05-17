#!/usr/bin/env python3
"""Task 1: External-only TNBC sensitivity analysis (Table S15).

Meta-analyzed delta C-index of locked GBSA vs PAM50-ROR in TNBC,
restricted to external cohorts GSE96058 (n~55) and METABRIC (n~320).
Both DerSimonian-Laird random-effects and fixed-effect are reported.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from src.meta import dersimonian_laird  # noqa: E402
from src.survival import (  # noqa: E402
    bootstrap_c_index_ci,
    harrell_c_index,
    paired_bootstrap_delta_c_index,
    uno_c_index,
)

RESULTS = REPO_ROOT / "results"
SEED = 42
BOOTSTRAP_ITER = 1000
HEADLINE = "Gradient_Boosted_Survival"
PAM50 = "PAM50_ROR_official"
EXTERNAL_COHORTS = ["GSE96058", "METABRIC"]


def main() -> None:
    np.random.seed(SEED)
    preds = pd.read_csv(RESULTS / "06_external_model_predictions.csv")
    preds["os_days"] = pd.to_numeric(preds["os_days"], errors="coerce")
    preds["os_event"] = pd.to_numeric(preds["os_event"], errors="coerce").fillna(0).astype(int)
    preds["tnbc_flag"] = preds["tnbc_flag"].astype("boolean")

    # Filter to external TNBC only
    ext_tnbc = preds[
        preds["cohort"].isin(EXTERNAL_COHORTS)
        & preds["tnbc_flag"].fillna(False).astype(bool)
    ].copy()

    rows = []
    for cohort, cdf in ext_tnbc.groupby("cohort", sort=True):
        valid = cdf[
            cdf[HEADLINE].notna()
            & cdf[PAM50].notna()
            & cdf["os_days"].notna()
            & (cdf["os_days"] > 0)
            & cdf["os_event"].notna()
        ].copy()
        n = len(valid)
        events = int(valid["os_event"].sum())

        # Headline metrics
        hl_c = harrell_c_index(valid["os_days"], valid["os_event"], valid[HEADLINE])
        hl_ci_low, hl_ci_high = bootstrap_c_index_ci(
            valid["os_days"], valid["os_event"], valid[HEADLINE], BOOTSTRAP_ITER, SEED
        )
        hl_uno = uno_c_index(valid["os_days"], valid["os_event"], valid[HEADLINE])

        # PAM50 metrics
        p50_c = harrell_c_index(valid["os_days"], valid["os_event"], valid[PAM50])
        p50_ci_low, p50_ci_high = bootstrap_c_index_ci(
            valid["os_days"], valid["os_event"], valid[PAM50], BOOTSTRAP_ITER, SEED
        )
        p50_uno = uno_c_index(valid["os_days"], valid["os_event"], valid[PAM50])

        # Paired delta
        delta = paired_bootstrap_delta_c_index(
            valid["os_days"], valid["os_event"],
            valid[HEADLINE], valid[PAM50],
            n_bootstrap=BOOTSTRAP_ITER, seed=SEED,
        )

        rows.append({
            "cohort": cohort,
            "subgroup": "TNBC",
            "n": n,
            "events": events,
            "headline_harrell_c": hl_c,
            "headline_c_ci_low": hl_ci_low,
            "headline_c_ci_high": hl_ci_high,
            "headline_uno_c": hl_uno,
            "pam50_harrell_c": p50_c,
            "pam50_c_ci_low": p50_ci_low,
            "pam50_c_ci_high": p50_ci_high,
            "pam50_uno_c": p50_uno,
            "delta_cindex": delta["delta"],
            "delta_ci_low": delta["ci_low"],
            "delta_ci_high": delta["ci_high"],
            "delta_p": delta["p"],
        })

    table = pd.DataFrame(rows)

    # --- Meta-analysis ---
    # Recover SE from bootstrap 95% CI for each cohort
    effects = table["delta_cindex"].to_numpy()
    se_arr = (table["delta_ci_high"].to_numpy() - table["delta_ci_low"].to_numpy()) / (2 * 1.96)
    variances = se_arr ** 2

    # Random-effects (DerSimonian-Laird)
    re = dersimonian_laird(effects, variances)

    # Fixed-effect (use inverse-variance weights with tau2=0)
    wi = 1.0 / variances
    fe_effect = float(np.sum(wi * effects) / np.sum(wi))
    fe_se = float(np.sqrt(1.0 / np.sum(wi)))
    from scipy import stats as sp_stats
    fe_z = fe_effect / fe_se if fe_se > 0 else np.nan
    fe_p = float(2 * sp_stats.norm.sf(abs(fe_z))) if np.isfinite(fe_z) else np.nan

    # Add meta rows
    meta_rows = [
        {
            "cohort": "Meta_RE",
            "subgroup": "TNBC",
            "n": int(table["n"].sum()),
            "events": int(table["events"].sum()),
            "headline_harrell_c": np.nan,
            "headline_c_ci_low": np.nan,
            "headline_c_ci_high": np.nan,
            "headline_uno_c": np.nan,
            "pam50_harrell_c": np.nan,
            "pam50_c_ci_low": np.nan,
            "pam50_c_ci_high": np.nan,
            "pam50_uno_c": np.nan,
            "delta_cindex": re["effect"],
            "delta_ci_low": re["ci_low"],
            "delta_ci_high": re["ci_high"],
            "delta_p": re["p"],
            "tau2": re["tau2"],
            "i2": re["i2"],
        },
        {
            "cohort": "Meta_FE",
            "subgroup": "TNBC",
            "n": int(table["n"].sum()),
            "events": int(table["events"].sum()),
            "headline_harrell_c": np.nan,
            "headline_c_ci_low": np.nan,
            "headline_c_ci_high": np.nan,
            "headline_uno_c": np.nan,
            "pam50_harrell_c": np.nan,
            "pam50_c_ci_low": np.nan,
            "pam50_c_ci_high": np.nan,
            "pam50_uno_c": np.nan,
            "delta_cindex": fe_effect,
            "delta_ci_low": fe_effect - 1.96 * fe_se,
            "delta_ci_high": fe_effect + 1.96 * fe_se,
            "delta_p": fe_p,
            "tau2": 0.0,
            "i2": 0.0,
        },
    ]
    full = pd.concat([table, pd.DataFrame(meta_rows)], ignore_index=True)

    # Save outputs
    RESULTS.mkdir(exist_ok=True)
    full.to_csv(RESULTS / "Table_S15_external_only_tnbc.csv", index=False)

    # Plain-text summary
    lines = [
        "Table S15: External-only TNBC Sensitivity Analysis",
        "=" * 55,
        f"Cohorts: {', '.join(EXTERNAL_COHORTS)} (TCGA-BRCA excluded as training set)",
        "",
    ]
    for _, row in table.iterrows():
        lines.append(
            f"  {row['cohort']}: n={row['n']}, events={row['events']}, "
            f"GBSA C={row['headline_harrell_c']:.4f}, PAM50 C={row['pam50_harrell_c']:.4f}, "
            f"delta={row['delta_cindex']:+.4f} [{row['delta_ci_low']:.4f}, {row['delta_ci_high']:.4f}], "
            f"p={row['delta_p']:.4f}"
        )
    lines.extend([
        "",
        "Random-effects meta-analysis (DerSimonian-Laird):",
        f"  delta = {re['effect']:+.4f}, 95% CI [{re['ci_low']:.4f}, {re['ci_high']:.4f}], "
        f"p = {re['p']:.4f}, tau2 = {re['tau2']:.4f}, I2 = {re['i2']:.1%}",
        "",
        "Fixed-effect meta-analysis:",
        f"  delta = {fe_effect:+.4f}, 95% CI [{fe_effect - 1.96 * fe_se:.4f}, "
        f"{fe_effect + 1.96 * fe_se:.4f}], p = {fe_p:.4f}",
        "",
    ])
    txt = "\n".join(lines)
    (RESULTS / "Table_S15_external_only_tnbc.txt").write_text(txt)
    print(txt)


if __name__ == "__main__":
    main()
