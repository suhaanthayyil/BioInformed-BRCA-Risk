"""Tests for random-effects meta-analysis utilities."""

from __future__ import annotations

import pandas as pd

from src.meta import random_effects_from_ci


def test_random_effects_from_ci_recovers_common_positive_effect() -> None:
    rows = pd.DataFrame(
        {
            "effect": [0.04, 0.05, 0.03],
            "low": [0.01, 0.02, 0.00],
            "high": [0.07, 0.08, 0.06],
        }
    )

    result = random_effects_from_ci(rows, "effect", "low", "high")

    assert result["n_studies"] == 3
    assert 0.03 < result["effect"] < 0.05
    assert result["ci_low"] > 0
    assert result["p"] < 0.05
