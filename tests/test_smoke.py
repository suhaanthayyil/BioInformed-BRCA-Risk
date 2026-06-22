"""Recompute smoke test: headline discrimination from shipped weights + data.

Unlike check_submission_consistency.py (which re-reads committed result CSVs),
this test *recomputes* the per-cohort Harrell C-index for the headline model by
loading the committed Gradient-Boosted-Survival weight cold (no __main__ shim)
and scoring the committed pathway-feature / survival data, then asserts the
recomputed values match the locked Table_1 results. It is the end-to-end proof
that the shipped weights + intermediates reproduce the manuscript's headline
numbers on a fresh clone.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

HEADLINE = "Gradient_Boosted_Survival"
TABLE1 = REPO_ROOT / "results" / "Table_1_ml_external_validation.csv"
NEEDED = [
    REPO_ROOT / "models" / "gradient_boosted_survival.pkl",
    REPO_ROOT / "data" / "processed" / "04_pathway_features.parquet",
    REPO_ROOT / "data" / "processed" / "unified_cohorts.duckdb",
    TABLE1,
]

pytestmark = pytest.mark.skipif(
    not all(p.exists() for p in NEEDED),
    reason="committed weights/intermediates not present",
)


def _expected_harrell() -> dict[str, float]:
    table = pd.read_csv(TABLE1)
    rows = table[(table["model"] == HEADLINE) & (table["subgroup"] == "overall")]
    return dict(zip(rows["cohort"], rows["harrell_c"]))


def test_shipped_gbsa_reproduces_per_cohort_cindex() -> None:
    from scripts.external_validation import load_frame
    from scripts.predict import load_artifact, score
    from src.survival import harrell_c_index

    frame = load_frame()
    artifact = load_artifact("gradient_boosted_survival")  # cold load, no shim
    expected = _expected_harrell()
    assert expected, "no expected Table_1 rows found"

    recomputed: dict[str, float] = {}
    for cohort, cohort_df in frame.groupby("cohort", sort=True):
        risk = score(artifact, cohort_df)
        sub = cohort_df.assign(risk=risk)
        sub = sub[sub["risk"].notna() & sub["os_days"].notna() & sub["os_event"].notna() & (sub["os_days"] > 0)]
        if len(sub) < 10 or int(sub["os_event"].sum()) < 2:
            continue
        recomputed[cohort] = float(harrell_c_index(sub["os_days"], sub["os_event"], sub["risk"]))

    # Every locked cohort value must be reproduced from weights + data.
    for cohort, exp in expected.items():
        assert cohort in recomputed, f"could not recompute {cohort}"
        assert recomputed[cohort] == pytest.approx(exp, abs=1e-4), (
            f"{cohort}: recomputed {recomputed[cohort]:.6f} != locked {exp:.6f}"
        )


def test_headline_internal_cv_value_is_locked() -> None:
    """Guard the headline internal CV C-index recorded in the model metadata."""
    import json

    meta = json.loads((REPO_ROOT / "data" / "processed" / "ml_model_zoo.metadata.json").read_text())
    assert meta["headline_model"] == HEADLINE
    assert meta["headline_cv_cindex"] == pytest.approx(0.6418, abs=5e-4)
    assert meta["headline_delta_vs_cox"] == pytest.approx(0.0420, abs=5e-4)
