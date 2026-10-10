#!/usr/bin/env python3
"""Evaluate every endpoint in the 2026-05-16 registration (Table S18)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from scripts.dca_analysis import nri_idi_stats  # noqa: E402
from scripts.phase3_common import (  # noqa: E402
    BOOTSTRAP_ITER,
    COX_MODEL,
    HEADLINE_MODEL,
    PAM50_MODEL,
    PROCESSED,
    RESULTS,
    SEED,
    ipcw_horizon,
    ci,
    load_modeling_frame,
    load_predictions,
    paired_delta,
    predict_headline,
    valid_survival_subset,
)
from src.meta import random_effects_from_ci  # noqa: E402
from src.survival import harrell_c_index  # noqa: E402

EXTERNAL = ["GSE96058", "METABRIC", "GSE20685"]
NI_MARGIN = -0.02
COX_MARGIN = 0.02
IDI_THRESHOLD = 0.05
HORIZON_DAYS = 5 * 365


def one_sided_lower(effect: float, se: float) -> float:
    return effect - stats.norm.ppf(0.95) * se


def delta_rows(pred: pd.DataFrame, comparator: str, subgroup: str, cohorts: list[str]) -> pd.DataFrame:
    rows = []
    for cohort in cohorts:
        sub = pred[pred["cohort"].eq(cohort)]
        if subgroup == "tnbc":
            sub = sub[sub["tnbc_flag"].fillna(False).astype(bool)]
        elif subgroup == "er_positive":
            sub = sub[sub["er_status"].eq("Positive")]
        sub = valid_survival_subset(sub, [HEADLINE_MODEL, comparator])
        if len(sub) < 10 or sub["os_event"].sum() < 2:
            continue
        d = paired_delta(sub, HEADLINE_MODEL, comparator)
        rows.append(
            {
                "cohort": cohort,
                "subgroup": subgroup,
                "comparator": comparator,
                "n": len(sub),
                "events": int(sub["os_event"].sum()),
                "delta": d["delta"],
                "ci_low": d["ci_low"],
                "ci_high": d["ci_high"],
                "p": d["p"],
            }
        )
    return pd.DataFrame(rows)


def pooled(rows: pd.DataFrame) -> dict:
    if rows.empty:
        return {"n_studies": 0, "effect": np.nan, "se": np.nan, "ci_low": np.nan, "ci_high": np.nan, "p": np.nan}
    return random_effects_from_ci(rows, "delta", "ci_low", "ci_high")


def with_oncotype(pred: pd.DataFrame) -> pd.DataFrame:
    onc = pd.read_parquet(PROCESSED / "baselines_oncotype.parquet")
    onc = onc[onc["status"].eq("ok")][["sample_id", "score", "er_status"]].rename(columns={"score": "Oncotype_surrogate"})
    return pred.merge(onc, on="sample_id", how="left")


def permutation_importance(df: pd.DataFrame, features: list[str], rng: np.random.Generator, repeats: int = 20) -> dict:
    base = harrell_c_index(df["os_days"], df["os_event"], predict_headline(df))
    out = {}
    for feature in features:
        drops = []
        for _ in range(repeats):
            shuffled = df.copy()
            shuffled[feature] = rng.permutation(shuffled[feature].to_numpy())
            drops.append(base - harrell_c_index(df["os_days"], df["os_event"], predict_headline(shuffled)))
        out[feature] = float(np.mean(drops))
    return out


def attribution_by_subgroup() -> pd.DataFrame:
    frame = load_modeling_frame()
    frame = frame[frame["cohort"].isin(EXTERNAL) & frame["os_days"].gt(0)]
    meta = json.loads((PROCESSED / "ml_model_zoo.metadata.json").read_text())
    features = meta["features"]
    groups = {
        "tnbc": frame["tnbc_flag"].astype("boolean").fillna(False).astype(bool),
        "her2_positive": frame["her2_status"].eq("Positive"),
    }
    rng = np.random.default_rng(SEED)
    rows = []
    for group, mask in groups.items():
        for cohort in EXTERNAL:
            sub = frame[mask & frame["cohort"].eq(cohort)]
            if len(sub) < 30 or sub["os_event"].sum() < 10:
                continue
            imp = permutation_importance(sub, features, rng)
            ranked = sorted(imp, key=imp.get, reverse=True)
            for feature in features:
                rows.append(
                    {
                        "subgroup": group,
                        "cohort": cohort,
                        "n": len(sub),
                        "events": int(sub["os_event"].sum()),
                        "feature": feature,
                        "importance_c_drop": imp[feature],
                        "rank": ranked.index(feature) + 1,
                    }
                )
    return pd.DataFrame(rows)


def five_year_prob(model: CoxPHFitter, df: pd.DataFrame) -> np.ndarray:
    surv = model.predict_survival_function(df, times=[HORIZON_DAYS])
    return 1.0 - surv.to_numpy().ravel()


def combined_idi(pred: pd.DataFrame) -> pd.DataFrame:
    cols = [HEADLINE_MODEL, PAM50_MODEL, "age_at_dx"]
    data = pred.dropna(subset=[*cols, "os_days"]).copy()
    data = data[data["os_days"] > 0]
    train = data[data["cohort"].eq("TCGA-BRCA")]
    full = CoxPHFitter(penalizer=0.01).fit(
        train[[HEADLINE_MODEL, PAM50_MODEL, "age_at_dx", "os_days", "os_event"]], "os_days", "os_event"
    )
    ref = CoxPHFitter(penalizer=0.01).fit(train[[PAM50_MODEL, "os_days", "os_event"]], "os_days", "os_event")
    rows = []
    for cohort in EXTERNAL:
        sub = data[data["cohort"].eq(cohort)]
        if sub.empty:
            continue
        y, w = ipcw_horizon(sub, HORIZON_DAYS)
        p_full = five_year_prob(full, sub[[HEADLINE_MODEL, PAM50_MODEL, "age_at_dx"]])
        p_ref = five_year_prob(ref, sub[[PAM50_MODEL]])
        point = nri_idi_stats(y, w, p_full, p_ref, 0.1)["idi"]
        rng = np.random.default_rng(SEED)
        boots = []
        for _ in range(BOOTSTRAP_ITER):
            idx = rng.integers(0, len(y), len(y))
            boots.append(nri_idi_stats(y[idx], w[idx], p_full[idx], p_ref[idx], 0.1)["idi"])
        low, high = ci(boots)
        rows.append(
            {
                "cohort": cohort,
                "n": len(sub),
                "events_5y": int(y.sum()),
                "idi_combined_vs_pam50": point,
                "idi_ci_low": low,
                "idi_ci_high": high,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    pred = with_oncotype(load_predictions())
    samples = load_modeling_frame()[["sample_id", "er_status", "her2_status", "age_at_dx"]]
    pred = pred.drop(columns=[c for c in ["er_status", "age_at_dx"] if c in pred]).merge(samples, on="sample_id", how="left")
    primary = json.loads((PROCESSED / "phase6_primary_endpoint.json").read_text())
    internal = pd.read_csv(RESULTS / "Table_S2_ml_internal_cv.csv")
    internal_delta = float(internal.loc[internal["model"].eq(HEADLINE_MODEL), "delta_vs_cox"].iloc[0])

    tnbc_all = delta_rows(pred, PAM50_MODEL, "tnbc", ["TCGA-BRCA", *EXTERNAL])
    tnbc_ext = delta_rows(pred, PAM50_MODEL, "tnbc", EXTERNAL)
    overall_pam = delta_rows(pred, PAM50_MODEL, "overall", EXTERNAL)
    er_onc = delta_rows(pred, "Oncotype_surrogate", "er_positive", EXTERNAL)
    overall_cox = delta_rows(pred, COX_MODEL, "overall", EXTERNAL)

    m_ext = pooled(tnbc_ext)
    m_pam = pooled(overall_pam)
    m_onc = pooled(er_onc)
    attr = attribution_by_subgroup()
    idi = combined_idi(pred)

    immune_top = {}
    for group, g in attr.groupby("subgroup"):
        mean_imp = g.groupby("feature")["importance_c_drop"].mean().sort_values(ascending=False)
        immune_top[group] = mean_imp.index[0]

    rows = [
        {
            "endpoint": "Primary: TNBC delta C vs PAM50 ROR-S (all cohorts with TNBC; TCGA out-of-fold)",
            "registered_criterion": "random-effects delta >= +0.03 and p < 0.05",
            "estimate": primary["meta_effect"],
            "ci_low": primary["meta_ci_low"],
            "ci_high": primary["meta_ci_high"],
            "p": primary["meta_p"],
            "met": bool(primary["met"]),
            "cohorts": ", ".join(primary["cohorts_in_meta"]),
        },
        {
            "endpoint": "Primary sensitivity: TNBC delta C vs PAM50 ROR-S (external cohorts only)",
            "registered_criterion": "random-effects delta >= +0.03 and p < 0.05",
            "estimate": m_ext["effect"],
            "ci_low": m_ext["ci_low"],
            "ci_high": m_ext["ci_high"],
            "p": m_ext["p"],
            "met": bool(m_ext["effect"] >= 0.03 and m_ext["p"] < 0.05),
            "cohorts": ", ".join(tnbc_ext["cohort"]),
        },
        {
            "endpoint": "Primary claim: TNBC superiority to Oncotype DX surrogate",
            "registered_criterion": "superior in TNBC",
            "estimate": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "p": np.nan,
            "met": None,
            "cohorts": "not evaluable: the 21-gene surrogate is defined for ER-positive disease only",
        },
        {
            "endpoint": "Primary claim: overall non-inferiority to Oncotype DX surrogate (ER-positive)",
            "registered_criterion": f"one-sided 95% lower bound of delta > {NI_MARGIN}",
            "estimate": m_onc["effect"],
            "ci_low": one_sided_lower(m_onc["effect"], m_onc["se"]),
            "ci_high": np.nan,
            "p": m_onc["p"],
            "met": bool(one_sided_lower(m_onc["effect"], m_onc["se"]) > NI_MARGIN),
            "cohorts": ", ".join(er_onc["cohort"]),
        },
        {
            "endpoint": "Secondary (a): ML minus Cox >= +0.02 in at least one cohort",
            "registered_criterion": f"delta >= {COX_MARGIN} in >= 1 cohort",
            "estimate": float(overall_cox["delta"].max()) if not overall_cox.empty else np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "p": np.nan,
            "met": bool((overall_cox["delta"] >= COX_MARGIN).any() or internal_delta >= COX_MARGIN),
            "cohorts": "external: "
            + "; ".join(f"{r.cohort} {r.delta:+.3f}" for r in overall_cox.itertuples())
            + f"; internal CV {internal_delta:+.3f}",
        },
        {
            "endpoint": "Secondary (b): immune pathway dominant attribution in TNBC and HER2-positive",
            "registered_criterion": "Pathway_Immune ranked first in both subgroups",
            "estimate": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "p": np.nan,
            "met": bool(all(v == "Pathway_Immune" for v in immune_top.values()) and len(immune_top) == 2),
            "cohorts": "; ".join(f"{k}: top={v}" for k, v in immune_top.items()),
        },
        {
            "endpoint": "Secondary (c): combined ML + PAM50 ROR-S + age IDI over PAM50 ROR-S alone (5-year)",
            "registered_criterion": f"IDI > {IDI_THRESHOLD}",
            "estimate": float(idi["idi_combined_vs_pam50"].mean()) if not idi.empty else np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "p": np.nan,
            "met": bool((idi["idi_combined_vs_pam50"] > IDI_THRESHOLD).all()) if not idi.empty else False,
            "cohorts": "; ".join(
                f"{r.cohort} {r.idi_combined_vs_pam50:+.3f} [{r.idi_ci_low:+.3f}, {r.idi_ci_high:+.3f}]" for r in idi.itertuples()
            ),
        },
        {
            "endpoint": "Secondary (d): overall non-inferiority to PAM50 ROR-S (external cohorts)",
            "registered_criterion": f"one-sided 95% lower bound of delta > {NI_MARGIN}",
            "estimate": m_pam["effect"],
            "ci_low": one_sided_lower(m_pam["effect"], m_pam["se"]),
            "ci_high": np.nan,
            "p": m_pam["p"],
            "met": bool(one_sided_lower(m_pam["effect"], m_pam["se"]) > NI_MARGIN),
            "cohorts": ", ".join(overall_pam["cohort"]),
        },
    ]
    table = pd.DataFrame(rows)
    table.to_csv(RESULTS / "Table_S18_registered_endpoints.csv", index=False)
    pd.concat([tnbc_all, tnbc_ext.assign(set="external_only"), overall_pam, er_onc, overall_cox], ignore_index=True).to_csv(
        RESULTS / "Table_S18_registered_endpoint_rows.csv", index=False
    )
    attr.to_csv(RESULTS / "Table_S18_subgroup_attribution.csv", index=False)
    idi.to_csv(RESULTS / "Table_S18_combined_idi.csv", index=False)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
