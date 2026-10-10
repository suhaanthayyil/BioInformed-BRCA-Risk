"""Recompute the pre-registered TNBC primary endpoint from per-cohort deltas.

Feeds the committed per-cohort TNBC delta C-indexes (headline model vs genefu
PAM50-ROR) through the same random-effects meta-analysis used in the pipeline
and asserts it matches the endpoint JSON written by external_validation.py.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

HEAD2HEAD = REPO_ROOT / "results" / "Table_2_head_to_head.csv"
PRIMARY = REPO_ROOT / "data" / "processed" / "phase6_primary_endpoint.json"
PAM50 = "PAM50_RORS_genefu"
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
    locked = json.loads(PRIMARY.read_text())

    assert meta["n_studies"] == 3
    assert meta["effect"] == pytest.approx(locked["meta_effect"], abs=1e-9)
    assert meta["p"] == pytest.approx(locked["meta_p"], abs=1e-9)

    met = bool(meta["effect"] >= THRESHOLD_DELTA and meta["p"] < THRESHOLD_P)
    assert met is locked["met"]
