from __future__ import annotations

import matplotlib.figure
import pandas as pd

from scripts import calibration_analysis, dca_analysis, within_subtype_analysis


def _capture_legend_labels(monkeypatch):
    captured: list[list[str]] = []

    def capture(fig, *_args, **_kwargs):
        captured.append([text.get_text() for legend in fig.legends for text in legend.get_texts()])

    monkeypatch.setattr(matplotlib.figure.Figure, "savefig", capture)
    return captured


def test_calibration_legend_uses_publication_labels(monkeypatch) -> None:
    captured = _capture_legend_labels(monkeypatch)
    cohort = calibration_analysis.ALL_COHORTS[0]
    points = pd.DataFrame({"pred": [0.1, 0.2], "obs": [0.12, 0.19]})
    curves = {(cohort, 5, model): points for model in calibration_analysis.MODELS}

    calibration_analysis.plot_calibration(curves)

    assert captured
    assert captured[0] == ["Pathway ML", "Cox baseline", "PAM50 ROR-S"]


def test_dca_legend_uses_pam50_ror_s(monkeypatch) -> None:
    captured = _capture_legend_labels(monkeypatch)
    cohort = dca_analysis.ALL_COHORTS[0]
    dca = pd.DataFrame(
        {
            "cohort": [cohort] * 5,
            "strategy": ["Treat all", "Treat none", "Pathway ML", "Cox baseline", "PAM50-ROR"],
            "threshold": [0.1] * 5,
            "net_benefit": [0.1, 0.0, 0.08, 0.07, 0.09],
        }
    )

    dca_analysis.plot_dca(dca)

    assert captured
    assert "PAM50 ROR-S" in captured[0]
    assert "PAM50-ROR" not in captured[0]


def test_within_subtype_forest_keeps_labels_off_data_and_all_intervals_visible(monkeypatch) -> None:
    captured: dict[str, object] = {}

    def capture(fig, *_args, **_kwargs):
        ax = fig.axes[0]
        captured["xlabel"] = ax.get_xlabel()
        captured["xlim"] = ax.get_xlim()
        captured["free_text_count"] = len(ax.texts)
        captured["tick_labels"] = [label.get_text() for label in ax.get_yticklabels()]

    monkeypatch.setattr(matplotlib.figure.Figure, "savefig", capture)
    table = pd.DataFrame(
        [
            {
                "subtype": "Basal",
                "cohort": "GSE20685",
                "n": 48,
                "events": 10,
                "delta_vs_pam50": 0.36,
                "delta_ci_low": 0.20,
                "delta_ci_high": 0.45,
                "status": "ok",
            }
        ]
    )
    meta = pd.DataFrame(
        [
            {
                "subtype": "Basal",
                "meta_delta": 0.14,
                "meta_ci_low": -0.01,
                "meta_ci_high": 0.30,
                "meta_p": 0.064,
            }
        ]
    )

    within_subtype_analysis.plot_forest(table, meta)

    assert captured["xlabel"] == "Delta C-index: pathway ML minus PAM50 ROR-S"
    assert captured["xlim"][1] > 0.45
    assert captured["free_text_count"] == 0
    assert any("GSE20685" in label for label in captured["tick_labels"])
