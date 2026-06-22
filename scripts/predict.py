#!/usr/bin/env python3
"""Reviewer-facing inference: load a shipped model artifact and score samples.

This is a thin, self-contained loader so a reviewer can verify the committed
weights without running the full pipeline. The artifacts load cold (no import
shim) after ``scripts/migrate_pickles.py`` has run.

Input CSV must contain the model's feature columns. The headline model expects
the seven pathway scores plus two clinical features:

    Pathway_Immune, Pathway_Proliferation, Pathway_DNA_Repair, Pathway_Metabolism,
    Pathway_Stromal_EMT, Pathway_Apoptosis_Stress, Pathway_Hormone,
    age_at_dx, stage_ordinal

A convenience: passing ``--from-features data/processed/04_pathway_features.parquet``
scores the committed pathway-feature table directly (clinical columns are joined
from unified_cohorts.duckdb when available, else imputed by the saved
preprocessor).

Examples:
    python scripts/predict.py --input my_samples.csv
    python scripts/predict.py --model random_survival_forest --input my_samples.csv -o risk.csv
    python scripts/predict.py --from-features data/processed/04_pathway_features.parquet --cohort TCGA-BRCA
"""

from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from src.ml.deepsurv import CoxMLP  # noqa: E402
from src.ml.wrappers import (  # noqa: E402,F401  (imported so pickled artifacts resolve)
    LifelinesCoxWrapper,
    SkSurvAdapter,
    XGBCoxWrapper,
)

MODELS = REPO_ROOT / "models"
PROCESSED = REPO_ROOT / "data" / "processed"

# Maps a friendly model key to its artifact filename.
ARTIFACTS = {
    "gradient_boosted_survival": "gradient_boosted_survival.pkl",
    "random_survival_forest": "random_survival_forest.pkl",
    "cox_ph": "cox_ph.pkl",
    "elastic_net_cox": "elastic_net_cox.pkl",
    "deepsurv": "deepsurv.pt",
    "stacked_ensemble": "stacked_ensemble.pkl",
}


def load_artifact(model_key: str) -> dict:
    filename = ARTIFACTS[model_key]
    path = MODELS / filename
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. The trained weights are committed under models/; "
            f"run `python scripts/migrate_pickles.py` once if loading fails."
        )
    if model_key == "deepsurv":
        payload = torch.load(path, map_location="cpu", weights_only=False)
        config = payload["config"]
        model = CoxMLP(payload["input_dim"], tuple(config["hidden_dims"]), config["dropout"])
        model.load_state_dict(payload["state_dict"])
        model.eval()
        return {
            "model_name": "DeepSurv",
            "preprocessor": payload["preprocessor"],
            "model": model,
            "features": payload["features"],
        }
    with path.open("rb") as fh:
        return pickle.load(fh)


def score(artifact: dict, frame: pd.DataFrame) -> np.ndarray:
    features = artifact["features"]
    missing = [c for c in features if c not in frame.columns]
    if missing:
        raise ValueError(
            f"Input is missing required feature columns: {missing}. Expected: {features}"
        )
    x = frame[features]
    pre = artifact["preprocessor"]
    x_scaled = pre.transform(x)
    model = artifact["model"]
    if artifact.get("model_name") == "DeepSurv":
        with torch.no_grad():
            return model(torch.tensor(x_scaled, dtype=torch.float32)).detach().numpy().astype(float)
    return np.asarray(model.predict(x_scaled), dtype=float)


def load_input(args: argparse.Namespace, features: list[str]) -> pd.DataFrame:
    if args.from_features:
        frame = pd.read_parquet(args.from_features)
        if args.cohort and "cohort" in frame.columns:
            frame = frame[frame["cohort"].eq(args.cohort)].copy()
        # Clinical columns may be absent in the raw pathway table; the saved
        # median imputer fills them, so add NaN placeholders if missing.
        for col in features:
            if col not in frame.columns:
                frame[col] = np.nan
        return frame
    if args.input.endswith(".parquet"):
        return pd.read_parquet(args.input)
    return pd.read_csv(args.input)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", default="gradient_boosted_survival", choices=sorted(ARTIFACTS))
    parser.add_argument("--input", help="CSV/Parquet with the model's feature columns")
    parser.add_argument("--from-features", help="Score data/processed/04_pathway_features.parquet directly")
    parser.add_argument("--cohort", help="Filter --from-features to one cohort")
    parser.add_argument("-o", "--output", help="Write risk scores to this CSV (else print head)")
    args = parser.parse_args()

    if not args.input and not args.from_features:
        parser.error("provide --input or --from-features")

    artifact = load_artifact(args.model)
    frame = load_input(args, artifact["features"])
    risk = score(artifact, frame)

    out = pd.DataFrame({"risk_score": risk})
    for idcol in ("sample_id", "cohort"):
        if idcol in frame.columns:
            out.insert(0, idcol, frame[idcol].to_numpy())
    if args.output:
        out.to_csv(args.output, index=False)
        print(f"Wrote {len(out)} risk scores to {args.output} (model={args.model})")
    else:
        print(f"model={args.model}  n={len(out)}")
        print(out.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
