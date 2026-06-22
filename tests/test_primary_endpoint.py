"""Recompute the pre-registered TNBC primary endpoint from per-cohort deltas.

Feeds the committed per-cohort TNBC delta C-indexes (headline model vs official
PAM50-ROR) through the same random-effects meta-analysis used in the pipeline
and asserts the locked result: effect ~ +0.0144, p ~ 0.6466, and that the
pre-registered success criterion (delta >= +0.03 AND p < 0.05) is NOT met.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

HEAD2HEAD = REPO_ROOT / "results" / "Table_2_head_to_head.csv"
PAM50 = "PAM50_ROR_official"
THRESHOLD_DELTA = 0.03
THRESHOLD_P = 0.05

pytestmark = pytest.mark.skipif(not HEAD2HEAD.exists(), reason="Table_2 not present")


def test_tnbc_primary_endpoint_recomputed_not_met() -> None:
    from src.meta import random_effects_from_ci

    h2h = pd.read_csv(HEAD2HEAD)
    tnbc = h2h[(h2h["subgroup"] == "tnbc") & (h2h["comparator"] == PAM50)].dropna(
        subset=["delta_cindex", "ci_low", "ci_high"]
    )
    # Cohorts with evaluable TNBC subgroups (GSE20685 has no TNBC labels -> dropped).
    assert len(tnbc) == 3, f"expected 3 TNBC cohorts, got {len(tnbc)}"

    meta = random_effects_from_ci(tnbc, "delta_cindex", "ci_low", "ci_high")

    assert meta["n_studies"] == 3
    assert meta["effect"] == pytest.approx(0.0144, abs=1e-3)
    assert meta["p"] == pytest.approx(0.6466, abs=1e-3)

    met = bool(meta["effect"] >= THRESHOLD_DELTA and meta["p"] < THRESHOLD_P)
    assert met is False, "pre-registered TNBC superiority endpoint must read as NOT met"
