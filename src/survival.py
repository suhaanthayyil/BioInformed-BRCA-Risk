"""
Cox proportional hazards survival analysis.

Implements Cox PH model fitting, cross-validated C-index evaluation,
and Kaplan-Meier visualization with log-rank testing.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.utils import concordance_index
from lifelines.statistics import logrank_test
from scipy import stats

try:
    from sksurv.metrics import (
        brier_score,
        concordance_index_ipcw,
        cumulative_dynamic_auc,
    )
    from sksurv.util import Surv
except Exception:  # pragma: no cover - optional dependency guard
    brier_score = None
    concordance_index_ipcw = None
    cumulative_dynamic_auc = None
    Surv = None


def fit_cox_model(X, time, event, penalizer=0.1):
    """
    Fit a Cox proportional hazards model.

    Args:
        X: Feature matrix (DataFrame).
        time: Time-to-event values (Series or array).
        event: Event indicator (1=event, 0=censored).
        penalizer: L2 regularization strength.

    Returns:
        Fitted CoxPHFitter instance.
    """
    df = X.copy()
    df["T"] = np.array(time)
    df["E"] = np.array(event)

    cph = CoxPHFitter(penalizer=penalizer)
    cph.fit(df, duration_col="T", event_col="E")
    return cph


def cox_cv(X, time, event, n_splits=5, penalizer=0.1):
    """
    Evaluate Cox model using stratified k-fold cross-validation.

    C-index is computed on each test fold. Stratification is on the
    binary event indicator.

    Args:
        X: Feature matrix (DataFrame or ndarray).
        time: Time-to-event values.
        event: Event indicator (binary).
        n_splits: Number of CV folds.
        penalizer: Cox PH regularization.

    Returns:
        Tuple of (mean_cindex, std_cindex, list_of_fold_cindexes).
    """
    X_arr = np.array(X)
    time_arr = np.array(time)
    event_arr = np.array(event).astype(int)
    feature_names = X.columns.tolist() if hasattr(X, "columns") else [f"f{i}" for i in range(X_arr.shape[1])]

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    fold_cindexes = []

    for train_idx, test_idx in skf.split(X_arr, event_arr):
        X_train, X_test = X_arr[train_idx], X_arr[test_idx]
        t_train, t_test = time_arr[train_idx], time_arr[test_idx]
        e_train, e_test = event_arr[train_idx], event_arr[test_idx]

        # Standardize within fold
        scaler = StandardScaler()
        X_train_s = scaler.fit_transform(X_train)
        X_test_s = scaler.transform(X_test)

        # Fit Cox on training fold
        train_df = pd.DataFrame(X_train_s, columns=feature_names)
        train_df["T"] = t_train
        train_df["E"] = e_train

        cph = CoxPHFitter(penalizer=penalizer)
        cph.fit(train_df, duration_col="T", event_col="E")

        # Evaluate on test fold
        test_df = pd.DataFrame(X_test_s, columns=feature_names)
        test_df["T"] = t_test
        test_df["E"] = e_test

        c_index = cph.score(test_df, scoring_method="concordance_index")
        fold_cindexes.append(c_index)

    return np.mean(fold_cindexes), np.std(fold_cindexes), fold_cindexes


def generate_km_curves(X, time, event, model, save_path=None):
    """
    Generate Kaplan-Meier survival curves for Cox-predicted risk groups.

    Patients are split at the median predicted partial hazard into
    high-risk and low-risk groups. Includes log-rank test p-value.

    Args:
        X: Feature matrix (DataFrame).
        time: Time-to-event values.
        event: Event indicator.
        model: Fitted CoxPHFitter.
        save_path: Optional path to save the figure.

    Returns:
        matplotlib Figure, log-rank test p-value.
    """
    # Predict partial hazards
    partial_hazards = model.predict_partial_hazard(X)
    median_hazard = partial_hazards.median()

    high_risk = partial_hazards >= median_hazard
    low_risk = ~high_risk

    time_arr = np.array(time)
    event_arr = np.array(event)

    # Filter to positive times
    valid = time_arr > 0
    time_arr = time_arr[valid]
    event_arr = event_arr[valid]
    high_risk = high_risk[valid]
    low_risk = low_risk[valid]

    fig, ax = plt.subplots(figsize=(8, 6))
    kmf = KaplanMeierFitter()

    # High risk group
    kmf.fit(
        time_arr[high_risk],
        event_observed=event_arr[high_risk],
        label="High Risk",
    )
    kmf.plot_survival_function(ax=ax, ci_show=True, color="red")

    # Low risk group
    kmf.fit(
        time_arr[low_risk],
        event_observed=event_arr[low_risk],
        label="Low Risk",
    )
    kmf.plot_survival_function(ax=ax, ci_show=True, color="blue")

    # Log-rank test
    lr_result = logrank_test(
        time_arr[high_risk], time_arr[low_risk],
        event_observed_A=event_arr[high_risk],
        event_observed_B=event_arr[low_risk],
    )

    ax.set_xlabel("Time (months)", fontsize=12)
    ax.set_ylabel("Survival Probability", fontsize=12)
    ax.set_title(
        f"Kaplan-Meier Survival Curves by Predicted Risk Group\n"
        f"Log-rank p = {lr_result.p_value:.2e}",
        fontsize=13,
    )
    ax.legend(fontsize=11)
    ax.set_ylim(0, 1.05)

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Saved KM curves to {save_path}")

    return fig, lr_result.p_value


def survival_arrays(time, event):
    """Return positive-time numeric survival arrays."""

    time_arr = pd.to_numeric(pd.Series(time), errors="coerce").to_numpy(float)
    event_arr = pd.to_numeric(pd.Series(event), errors="coerce").fillna(0).to_numpy(int)
    valid = np.isfinite(time_arr) & (time_arr > 0) & np.isfinite(event_arr)
    return time_arr[valid], event_arr[valid].astype(bool), valid


def harrell_c_index(time, event, risk):
    """Harrell C-index where larger risk means shorter survival."""

    time_arr, event_arr, valid_surv = survival_arrays(time, event)
    risk_arr = pd.to_numeric(pd.Series(risk), errors="coerce").to_numpy(float)
    risk_arr = risk_arr[valid_surv]
    valid = np.isfinite(risk_arr)
    if valid.sum() < 5 or event_arr[valid].sum() < 2:
        return np.nan
    return float(concordance_index(time_arr[valid], -risk_arr[valid], event_arr[valid].astype(int)))


def uno_c_index(time, event, risk, tau=None):
    """Uno IPCW C-index where larger risk means shorter survival."""

    if Surv is None or concordance_index_ipcw is None:
        return np.nan
    time_arr, event_arr, valid_surv = survival_arrays(time, event)
    risk_arr = pd.to_numeric(pd.Series(risk), errors="coerce").to_numpy(float)
    risk_arr = risk_arr[valid_surv]
    valid = np.isfinite(risk_arr)
    if valid.sum() < 10 or event_arr[valid].sum() < 2:
        return np.nan
    y = Surv.from_arrays(event_arr[valid], time_arr[valid])
    if tau is None:
        tau = np.nanpercentile(time_arr[valid], 80)
    try:
        return float(concordance_index_ipcw(y, y, risk_arr[valid], tau=tau)[0])
    except Exception:
        return np.nan


def bootstrap_c_index_ci(time, event, risk, n_bootstrap=1000, seed=42, alpha=0.05):
    """Bootstrap Harrell C-index CI for a single risk score."""

    rng = np.random.default_rng(seed)
    time_arr, event_arr, valid_surv = survival_arrays(time, event)
    risk_arr = pd.to_numeric(pd.Series(risk), errors="coerce").to_numpy(float)
    risk_arr = risk_arr[valid_surv]
    valid = np.isfinite(risk_arr)
    time_arr = time_arr[valid]
    event_arr = event_arr[valid]
    risk_arr = risk_arr[valid]
    n = len(time_arr)
    if n < 10 or event_arr.sum() < 2:
        return np.nan, np.nan

    values = []
    for _ in range(n_bootstrap):
        idx = rng.integers(0, n, n)
        if event_arr[idx].sum() < 2 or len(np.unique(time_arr[idx])) < 3:
            continue
        values.append(harrell_c_index(time_arr[idx], event_arr[idx], risk_arr[idx]))
    values = np.array([v for v in values if np.isfinite(v)])
    if len(values) < 25:
        return np.nan, np.nan
    return tuple(np.quantile(values, [alpha / 2, 1 - alpha / 2]).astype(float))


def time_dependent_auc(time, event, risk, years=(3, 5, 10)):
    """Compute cumulative dynamic AUC at requested year cutoffs."""

    out = {f"auc_{year}y": np.nan for year in years}
    if Surv is None or cumulative_dynamic_auc is None:
        return out
    time_arr, event_arr, valid_surv = survival_arrays(time, event)
    risk_arr = pd.to_numeric(pd.Series(risk), errors="coerce").to_numpy(float)
    risk_arr = risk_arr[valid_surv]
    valid = np.isfinite(risk_arr)
    time_arr = time_arr[valid]
    event_arr = event_arr[valid]
    risk_arr = risk_arr[valid]
    if len(time_arr) < 20 or event_arr.sum() < 3:
        return out
    y = Surv.from_arrays(event_arr, time_arr)
    eval_times = np.array([365.25 * year for year in years], dtype=float)
    allowed = eval_times[(eval_times > np.min(time_arr)) & (eval_times < np.max(time_arr))]
    if len(allowed) == 0:
        return out
    try:
        aucs, _ = cumulative_dynamic_auc(y, y, risk_arr, allowed)
    except Exception:
        return out
    for year, eval_time in zip(years, eval_times):
        matches = np.where(np.isclose(allowed, eval_time))[0]
        if len(matches):
            out[f"auc_{year}y"] = float(aucs[matches[0]])
    return out


def brier_score_at(time, event, risk, years=(5,)):
    """Approximate IPCW Brier scores by evaluating a logistic risk transform."""

    out = {f"brier_{year}y": np.nan for year in years}
    if Surv is None or brier_score is None:
        return out
    time_arr, event_arr, valid_surv = survival_arrays(time, event)
    risk_arr = pd.to_numeric(pd.Series(risk), errors="coerce").to_numpy(float)
    risk_arr = risk_arr[valid_surv]
    valid = np.isfinite(risk_arr)
    time_arr = time_arr[valid]
    event_arr = event_arr[valid]
    risk_arr = risk_arr[valid]
    if len(time_arr) < 20 or event_arr.sum() < 3:
        return out
    y = Surv.from_arrays(event_arr, time_arr)
    z = (risk_arr - np.nanmean(risk_arr)) / (np.nanstd(risk_arr) or 1.0)
    risk_prob = 1.0 / (1.0 + np.exp(-z))
    survival_prob = 1.0 - risk_prob
    eval_times = np.array([365.25 * year for year in years], dtype=float)
    allowed = eval_times[(eval_times > np.min(time_arr)) & (eval_times < np.max(time_arr))]
    if len(allowed) == 0:
        return out
    estimates = np.tile(survival_prob[:, None], (1, len(allowed)))
    try:
        _, scores = brier_score(y, y, estimates, allowed)
    except Exception:
        return out
    for year, eval_time in zip(years, eval_times):
        matches = np.where(np.isclose(allowed, eval_time))[0]
        if len(matches):
            out[f"brier_{year}y"] = float(scores[matches[0]])
    return out


def cox_hr_high_low(time, event, risk):
    """Fit a median-split Cox model and return HR, CI, and p-value."""

    time_arr, event_arr, valid_surv = survival_arrays(time, event)
    risk_arr = pd.to_numeric(pd.Series(risk), errors="coerce").to_numpy(float)
    risk_arr = risk_arr[valid_surv]
    valid = np.isfinite(risk_arr)
    if valid.sum() < 20 or event_arr[valid].sum() < 3:
        return {"hr": np.nan, "hr_ci_low": np.nan, "hr_ci_high": np.nan, "hr_p": np.nan}
    group = (risk_arr[valid] >= np.nanmedian(risk_arr[valid])).astype(int)
    if len(np.unique(group)) < 2:
        return {"hr": np.nan, "hr_ci_low": np.nan, "hr_ci_high": np.nan, "hr_p": np.nan}
    df = pd.DataFrame({"T": time_arr[valid], "E": event_arr[valid].astype(int), "high_risk": group})
    try:
        cph = CoxPHFitter(penalizer=0.01)
        cph.fit(df, "T", "E")
        row = cph.summary.loc["high_risk"]
        return {
            "hr": float(row["exp(coef)"]),
            "hr_ci_low": float(row["exp(coef) lower 95%"]),
            "hr_ci_high": float(row["exp(coef) upper 95%"]),
            "hr_p": float(row["p"]),
        }
    except Exception:
        return {"hr": np.nan, "hr_ci_low": np.nan, "hr_ci_high": np.nan, "hr_p": np.nan}


def paired_bootstrap_delta_c_index(time, event, risk_a, risk_b, n_bootstrap=1000, seed=42):
    """Bootstrap paired delta C-index (A minus B)."""

    rng = np.random.default_rng(seed)
    time_arr, event_arr, valid_surv = survival_arrays(time, event)
    a = pd.to_numeric(pd.Series(risk_a), errors="coerce").to_numpy(float)[valid_surv]
    b = pd.to_numeric(pd.Series(risk_b), errors="coerce").to_numpy(float)[valid_surv]
    valid = np.isfinite(a) & np.isfinite(b)
    time_arr = time_arr[valid]
    event_arr = event_arr[valid]
    a = a[valid]
    b = b[valid]
    if len(time_arr) < 20 or event_arr.sum() < 3:
        return {"delta": np.nan, "ci_low": np.nan, "ci_high": np.nan, "p": np.nan}
    delta = harrell_c_index(time_arr, event_arr, a) - harrell_c_index(time_arr, event_arr, b)
    boot = []
    n = len(time_arr)
    for _ in range(n_bootstrap):
        idx = rng.integers(0, n, n)
        if event_arr[idx].sum() < 3:
            continue
        ca = harrell_c_index(time_arr[idx], event_arr[idx], a[idx])
        cb = harrell_c_index(time_arr[idx], event_arr[idx], b[idx])
        if np.isfinite(ca) and np.isfinite(cb):
            boot.append(ca - cb)
    boot = np.array(boot)
    if len(boot) < 25:
        return {"delta": delta, "ci_low": np.nan, "ci_high": np.nan, "p": np.nan}
    ci_low, ci_high = np.quantile(boot, [0.025, 0.975])
    se = np.std(boot, ddof=1)
    p = 2 * stats.norm.sf(abs(delta / se)) if se > 0 else np.nan
    return {"delta": float(delta), "ci_low": float(ci_low), "ci_high": float(ci_high), "p": float(p)}
