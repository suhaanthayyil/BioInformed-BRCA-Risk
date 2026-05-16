"""Regression tests for the official genefu PAM50/ROR repair."""

from __future__ import annotations

import numpy as np
import pandas as pd

from scripts.compute_clinical_baselines import load_pam50_map, official_pam50_score


def test_load_pam50_map_from_genefu() -> None:
    pam50_map = load_pam50_map()

    assert len(pam50_map) >= 50
    assert {"probe", "EntrezGene.ID"}.issubset(pam50_map.columns)
    assert pam50_map["probe"].str.len().gt(0).all()


def test_official_pam50_score_returns_genefu_rors_rows() -> None:
    pam50_map = load_pam50_map()
    probes = pam50_map["probe"].head(45).tolist()
    rng = np.random.default_rng(42)
    expression = pd.DataFrame(
        rng.normal(size=(8, len(probes))),
        index=[f"S{i}" for i in range(8)],
        columns=probes,
    )

    scores = official_pam50_score("TEST", expression, pam50_map)

    assert len(scores) == len(expression)
    assert scores["baseline"].eq("PAM50_ROR_official").all()
    assert scores["status"].eq("ok").all()
    assert scores["method_label"].eq("genefu_rorS").all()
    assert scores["score"].between(0, 100).all()
    assert scores["pam50_subtype_genefu"].notna().all()
