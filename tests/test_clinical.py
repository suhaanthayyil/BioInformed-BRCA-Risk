from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.clinical import stage_to_ordinal  # noqa: E402


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("Stage I", 1.0),
        ("Stage IB", 1.0),
        ("Stage IIA", 2.0),
        ("Stage IIIC", 3.0),
        ("Stage IV", 4.0),
        ("2", 2.0),
        ("4", 4.0),
        ("T1.0N0M0", 1.0),
        ("T2.0N1M0", 2.0),
        ("T3.0N0M0", 2.0),
        ("T2.0N2M0", 3.0),
        ("T4.0N0M0", 3.0),
        ("T3.0N3M1", 4.0),
    ],
)
def test_stage_formats_from_every_cohort(value: str, expected: float) -> None:
    assert stage_to_ordinal(value) == expected


@pytest.mark.parametrize("value", ["No_Conversion", "Stage X", "Stage", "[Discrepancy]", "0", None, ""])
def test_unparseable_stage_is_missing(value) -> None:
    assert math.isnan(stage_to_ordinal(value))
