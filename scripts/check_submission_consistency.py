#!/usr/bin/env python3
"""Check repository values that support the submitted BMC manuscript."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PROCESSED = ROOT / "data" / "processed"


def close(actual: float, expected: float, tolerance: float = 5e-4) -> bool:
    return bool(np.isfinite(actual) and abs(actual - expected) <= tolerance)


def check(condition: bool, label: str, actual: object, expected: object) -> dict[str, object]:
    return {
        "check": label,
        "status": "PASS" if condition else "FAIL",
        "actual": actual,
        "expected": expected,
    }


def main() -> None:
    checks: list[dict[str, object]] = []

    primary = json.loads((PROCESSED / "phase6_primary_endpoint.json").read_text())
    checks.append(check(close(primary["meta_effect"], 0.0143978), "TNBC primary delta", primary["meta_effect"], 0.0143978))
    checks.append(check(close(primary["meta_p"], 0.6466270), "TNBC primary p-value", primary["meta_p"], 0.6466270))
    checks.append(check(primary["met"] is False, "TNBC primary endpoint status", primary["met"], False))

    internal = pd.read_csv(RESULTS / "Table_S2_ml_internal_cv.csv")
    headline = internal.loc[internal["model"].eq("Gradient_Boosted_Survival")].iloc[0]
    checks.append(
        check(
            close(float(headline["mean_cv_cindex"]), 0.641826, 1e-6),
            "Headline internal CV C-index",
            float(headline["mean_cv_cindex"]),
            0.641826,
        )
    )
    checks.append(
        check(
            close(float(headline["delta_vs_cox"]), 0.041964, 1e-6),
            "Headline internal delta versus Cox",
            float(headline["delta_vs_cox"]),
            0.041964,
        )
    )

    external = pd.read_csv(RESULTS / "Table_2_external_validation.csv")
    expected_external = {
        ("GSE20685", "Gradient_Boosted_Survival"): 0.565182,
        ("GSE20685", "PAM50_ROR_official"): 0.635178,
        ("GSE96058", "Gradient_Boosted_Survival"): 0.662624,
        ("GSE96058", "PAM50_ROR_official"): 0.632816,
        ("METABRIC", "Gradient_Boosted_Survival"): 0.554892,
        ("METABRIC", "PAM50_ROR_official"): 0.594897,
        ("TCGA-BRCA", "Gradient_Boosted_Survival"): 0.672968,
        ("TCGA-BRCA", "PAM50_ROR_official"): 0.504266,
    }
    for (cohort, model), expected in expected_external.items():
        row = external.loc[
            external["cohort"].eq(cohort)
            & external["model"].eq(model)
            & external["subgroup"].eq("overall")
            & external["status"].eq("ok")
        ].iloc[0]
        checks.append(
            check(
                close(float(row["harrell_c"]), expected, 1e-6),
                f"External Harrell C {cohort} {model}",
                float(row["harrell_c"]),
                expected,
            )
        )

    calibration = pd.read_csv(RESULTS / "Table_S10_calibration.csv")
    cal5 = calibration.loc[calibration["time_horizon"].eq("5y")]
    ml_mean_ici = float(cal5.loc[cal5["model"].eq("Gradient_Boosted_Survival"), "ici"].mean())
    pam_mean_ici = float(cal5.loc[cal5["model"].eq("PAM50_ROR_official"), "ici"].mean())
    checks.append(check(close(ml_mean_ici, 0.125013, 1e-6), "Mean 5-year ICI headline model", ml_mean_ici, 0.125013))
    checks.append(check(close(pam_mean_ici, 0.154658, 1e-6), "Mean 5-year ICI PAM50-ROR", pam_mean_ici, 0.154658))

    dca = pd.read_csv(RESULTS / "Table_S11_dca.csv")
    wide = dca.pivot_table(index=["cohort", "threshold"], columns="strategy", values="net_benefit")
    mean_net_benefit_delta = float((wide["Pathway ML"] - wide["PAM50-ROR"]).mean())
    checks.append(
        check(
            close(mean_net_benefit_delta, 0.015560, 1e-6),
            "Mean DCA net-benefit delta",
            mean_net_benefit_delta,
            0.015560,
        )
    )

    stability = pd.read_csv(RESULTS / "Table_S13_stability.csv").iloc[0]
    checks.append(
        check(
            close(float(stability["mean_pairwise_spearman"]), 0.395355, 1e-6),
            "Stability mean Spearman rho",
            float(stability["mean_pairwise_spearman"]),
            0.395355,
        )
    )

    table = pd.DataFrame(checks)
    print(table.to_string(index=False))
    failed = table.loc[table["status"].eq("FAIL")]
    if not failed.empty:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
