#!/usr/bin/env python3
"""Load a shipped model artifact and score samples.

Input must contain the model's feature columns (the seven pathway scores and
age_at_dx). --from-features scores data/processed/04_pathway_features.parquet
with age joined from unified_cohorts.duckdb.

Examples:
    python scripts/predict.py --input my_samples.csv
    python scripts/predict.py --model random_survival_forest --input my_samples.csv -o risk.csv
    python scripts/predict.py --from-features data/processed/04_pathway_features.parquet --cohort TCGA-BRCA
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import duckdb

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from src.clinical import stage_to_ordinal  # noqa: E402
from src.ml import registry  # noqa: E402

ARTIFACTS = {name.lower(): name for name in registry.FILES if name != "Stacked_Ensemble"}


def load_artifact(model_key: str) -> dict:
    return registry.load_artifact(ARTIFACTS[model_key])


def score(artifact: dict, frame: pd.DataFrame) -> np.ndarray:
    missing = [c for c in artifact["features"] if c not in frame.columns]
    if missing:
        raise ValueError(f"Input is missing required feature columns: {missing}. Expected: {artifact['features']}")
    return registry.predict(artifact, frame)


def load_input(args: argparse.Namespace, features: list[str]) -> pd.DataFrame:
    if args.from_features:
        frame = pd.read_parquet(args.from_features)
        if args.cohort and "cohort" in frame.columns:
            frame = frame[frame["cohort"].eq(args.cohort)].copy()
        con = duckdb.connect(str(REPO_ROOT / "data" / "processed" / "unified_cohorts.duckdb"), read_only=True)
        try:
            samples = con.execute("select sample_id, age_at_dx, stage_tnm from samples").fetchdf()
        finally:
            con.close()
        samples["stage_ordinal"] = samples["stage_tnm"].map(stage_to_ordinal)
        return frame.merge(samples, on="sample_id", how="left")
    if args.input.endswith(".parquet"):
        return pd.read_parquet(args.input)
    return pd.read_csv(args.input)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", default=registry.headline_name().lower(), choices=sorted(ARTIFACTS))
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
