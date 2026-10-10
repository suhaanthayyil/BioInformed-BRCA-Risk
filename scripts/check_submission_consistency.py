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
    head = primary["headline_model"]
    checks.append(check(head == "DeepSurv", "Headline model", head, "DeepSurv"))
    checks.append(check(close(primary["meta_effect"], 0.0204261), "TNBC primary delta", primary["meta_effect"], 0.0204261))
    checks.append(check(close(primary["meta_p"], 0.5095809), "TNBC primary p-value", primary["meta_p"], 0.5095809))
    checks.append(check(primary["met"] is False, "TNBC primary endpoint status", primary["met"], False))

    internal = pd.read_csv(RESULTS / "Table_S2_ml_internal_cv.csv")
    row = internal.loc[internal["model"].eq(head)].iloc[0]
    checks.append(check(close(float(row["mean_cv_cindex"]), 0.592974, 1e-6), "Headline internal CV C-index", float(row["mean_cv_cindex"]), 0.592974))
    checks.append(check(close(float(row["delta_vs_cox"]), -0.021041, 1e-6), "Headline internal delta versus Cox", float(row["delta_vs_cox"]), -0.021041))

    external = pd.read_csv(RESULTS / "Table_2_external_validation.csv")
    expected_external = {
        ("GSE20685", head): 0.671836,
        ("GSE20685", "PAM50_RORS_genefu"): 0.635178,
        ("GSE96058", head): 0.568559,
        ("GSE96058", "PAM50_RORS_genefu"): 0.645972,
        ("METABRIC", head): 0.552385,
        ("METABRIC", "PAM50_RORS_genefu"): 0.594897,
        ("TCGA-BRCA", head): 0.591223,
        ("TCGA-BRCA", "PAM50_RORS_genefu"): 0.614782,
    }
    for (cohort, model), expected in expected_external.items():
        hit = external.loc[
            external["cohort"].eq(cohort)
            & external["model"].eq(model)
            & external["subgroup"].eq("overall")
            & external["status"].eq("ok")
        ].iloc[0]
        checks.append(check(close(float(hit["harrell_c"]), expected, 1e-6), f"Harrell C {cohort} {model}", float(hit["harrell_c"]), expected))

    calibration = pd.read_csv(RESULTS / "Table_S10_calibration.csv")
    cal5 = calibration.loc[calibration["time_horizon"].eq("5y")]
    ml_mean_ici = float(cal5.loc[cal5["model"].eq(head), "ici"].mean())
    pam_mean_ici = float(cal5.loc[cal5["model"].eq("PAM50_RORS_genefu"), "ici"].mean())
    checks.append(check(close(ml_mean_ici, 0.044659, 1e-6), "Mean 5-year ICI headline model", ml_mean_ici, 0.044659))
    checks.append(check(close(pam_mean_ici, 0.036063, 1e-6), "Mean 5-year ICI PAM50 ROR-S", pam_mean_ici, 0.036063))

    dca = pd.read_csv(RESULTS / "Table_S11_dca.csv")
    wide = dca.pivot_table(index=["cohort", "threshold"], columns="strategy", values="net_benefit")
    nb = float((wide["Pathway ML"] - wide["PAM50-ROR"]).mean())
    checks.append(check(close(nb, -0.002126, 1e-6), "Mean DCA net-benefit delta", nb, -0.002126))

    stability = pd.read_csv(RESULTS / "Table_S13_stability.csv").iloc[0]
    rho = float(stability["mean_pairwise_spearman"])
    checks.append(check(close(rho, 0.139458, 1e-6), "Stability mean Spearman rho", rho, 0.139458))

    table = pd.DataFrame(checks)
    print(table.to_string(index=False))
    if table["status"].eq("FAIL").any():
        raise SystemExit(1)


if __name__ == "__main__":
    main()
