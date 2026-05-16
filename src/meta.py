"""Meta-analysis helpers."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def dersimonian_laird(effects, variances):
    """DerSimonian-Laird random-effects meta-analysis."""

    yi = np.asarray(effects, dtype=float)
    vi = np.asarray(variances, dtype=float)
    valid = np.isfinite(yi) & np.isfinite(vi) & (vi > 0)
    yi = yi[valid]
    vi = vi[valid]
    k = len(yi)
    if k == 0:
        return {
            "n_studies": 0,
            "effect": np.nan,
            "se": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "z": np.nan,
            "p": np.nan,
            "tau2": np.nan,
            "i2": np.nan,
        }
    if k == 1:
        se = float(np.sqrt(vi[0]))
        effect = float(yi[0])
        z = effect / se if se > 0 else np.nan
        return {
            "n_studies": 1,
            "effect": effect,
            "se": se,
            "ci_low": effect - 1.96 * se,
            "ci_high": effect + 1.96 * se,
            "z": z,
            "p": float(2 * stats.norm.sf(abs(z))) if np.isfinite(z) else np.nan,
            "tau2": 0.0,
            "i2": 0.0,
        }

    wi = 1.0 / vi
    fixed = np.sum(wi * yi) / np.sum(wi)
    q = np.sum(wi * (yi - fixed) ** 2)
    c = np.sum(wi) - (np.sum(wi**2) / np.sum(wi))
    tau2 = max(0.0, (q - (k - 1)) / c) if c > 0 else 0.0
    w_star = 1.0 / (vi + tau2)
    effect = float(np.sum(w_star * yi) / np.sum(w_star))
    se = float(np.sqrt(1.0 / np.sum(w_star)))
    z = effect / se if se > 0 else np.nan
    i2 = max(0.0, (q - (k - 1)) / q) if q > 0 else 0.0
    return {
        "n_studies": int(k),
        "effect": effect,
        "se": se,
        "ci_low": effect - 1.96 * se,
        "ci_high": effect + 1.96 * se,
        "z": z,
        "p": float(2 * stats.norm.sf(abs(z))) if np.isfinite(z) else np.nan,
        "tau2": float(tau2),
        "i2": float(i2),
    }


def random_effects_from_ci(df: pd.DataFrame, effect_col: str, low_col: str, high_col: str):
    """Pool effects using standard errors recovered from 95% CIs."""

    tmp = df[[effect_col, low_col, high_col]].copy()
    tmp["se"] = (tmp[high_col] - tmp[low_col]) / (2 * 1.96)
    tmp["variance"] = tmp["se"] ** 2
    return dersimonian_laird(tmp[effect_col], tmp["variance"])
