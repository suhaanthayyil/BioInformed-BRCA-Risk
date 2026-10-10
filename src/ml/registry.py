"""Load, score and rebuild the trained survival models."""

from __future__ import annotations

import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import torch
from sksurv.ensemble import GradientBoostingSurvivalAnalysis, RandomSurvivalForest

from src.ml.deepsurv import CoxMLP, DeepSurvConfig, DeepSurvEstimator
from src.ml.wrappers import LifelinesCoxWrapper, SkSurvAdapter, XGBCoxWrapper  # noqa: F401

ROOT = Path(__file__).resolve().parents[2]
MODELS = ROOT / "models"
SEED = 42

FILES = {
    "Cox_PH": "cox_ph.pkl",
    "Elastic_Net_Cox": "elastic_net_cox.pkl",
    "Random_Survival_Forest": "random_survival_forest.pkl",
    "Gradient_Boosted_Survival": "gradient_boosted_survival.pkl",
    "DeepSurv": "deepsurv.pt",
    "Stacked_Ensemble": "stacked_ensemble.pkl",
}


def metadata() -> dict:
    return json.loads((ROOT / "data" / "processed" / "ml_model_zoo.metadata.json").read_text())


def headline_name() -> str:
    return metadata()["headline_model"]


def selected_config(name: str) -> dict:
    cv = pd.read_csv(ROOT / "results" / "Table_S2_ml_internal_cv.csv")
    return json.loads(cv.loc[cv["model"].eq(name), "config"].iloc[0])


def load_artifact(name: str, models_dir: Path = MODELS) -> dict:
    path = models_dir / FILES[name]
    if name == "DeepSurv":
        payload = torch.load(path, map_location="cpu", weights_only=False)
        config = payload["config"]
        model = CoxMLP(payload["input_dim"], tuple(config["hidden_dims"]), config["dropout"])
        model.load_state_dict(payload["state_dict"])
        model.eval()
        return {
            "model_name": "DeepSurv",
            "config": config,
            "preprocessor": payload["preprocessor"],
            "model": model,
            "features": payload["features"],
        }
    return joblib.load(path)


def predict(artifact: dict, frame: pd.DataFrame) -> np.ndarray:
    x = artifact["preprocessor"].transform(frame[artifact["features"]])
    if artifact.get("model_name") == "DeepSurv":
        with torch.no_grad():
            return artifact["model"](torch.tensor(x, dtype=torch.float32)).numpy().astype(float).ravel()
    inner = getattr(artifact["model"], "inner", None)
    if inner is not None and "n_jobs" in inner.get_params():
        inner.set_params(n_jobs=1)
    return np.asarray(artifact["model"].predict(x), dtype=float)


def make_estimator(name: str, config: dict, seed: int = SEED):
    if name == "Gradient_Boosted_Survival":
        return SkSurvAdapter(GradientBoostingSurvivalAnalysis(random_state=seed, **config))
    if name == "Random_Survival_Forest":
        return SkSurvAdapter(RandomSurvivalForest(random_state=seed, n_jobs=-1, **config))
    if name == "DeepSurv":
        cfg = dict(config)
        cfg["hidden_dims"] = tuple(cfg["hidden_dims"])
        return DeepSurvEstimator(DeepSurvConfig(**cfg))
    raise ValueError(f"No retraining recipe for {name}")
