#!/usr/bin/env python3
"""Pickleable survival-model wrappers used by the model zoo.

These adapters were historically defined inside ``scripts/train_ml_zoo.py`` and
were therefore pickled under the ``__main__`` namespace, which made the saved
model artifacts in ``models/`` impossible to ``pickle.load`` from a fresh
process without an import shim. They now live in this importable module so that
the shipped weights load cold via ``scripts/predict.py`` and
``scripts/external_validation.py`` without touching ``__main__``.

Each wrapper exposes the same ``fit(x, time, event)`` / ``predict(x)`` signature
so the training and validation code can treat heterogeneous survival estimators
(lifelines, XGBoost, scikit-survival) uniformly.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import xgboost as xgb
from lifelines import CoxPHFitter
from sksurv.util import Surv

# Single source of truth for the project random seed (mirrors train_ml_zoo.SEED).
SEED = 42


class LifelinesCoxWrapper:
    def __init__(self, penalizer: float):
        self.penalizer = penalizer
        self.model = CoxPHFitter(penalizer=penalizer)

    def fit(self, x, time, event):
        df = pd.DataFrame(x, columns=[f"x{i}" for i in range(x.shape[1])])
        df["T"] = time
        df["E"] = event.astype(int)
        self.model.fit(df, duration_col="T", event_col="E", show_progress=False)
        return self

    def predict(self, x):
        df = pd.DataFrame(x, columns=[f"x{i}" for i in range(x.shape[1])])
        return self.model.predict_partial_hazard(df).to_numpy(float)


class XGBCoxWrapper:
    def __init__(self, **params):
        self.params = params
        self.model = xgb.XGBRegressor(
            objective="survival:cox",
            eval_metric="cox-nloglik",
            random_state=SEED,
            tree_method="hist",
            n_jobs=1,
            verbosity=0,
            **params,
        )

    def fit(self, x, time, event):
        signed_time = np.where(event.astype(bool), time, -time)
        self.model.fit(x, signed_time, verbose=False)
        return self

    def predict(self, x):
        return self.model.predict(x)


class SkSurvAdapter:
    """Pickleable adapter giving sksurv estimators a common fit signature."""

    def __init__(self, estimator):
        self.inner = estimator

    def fit(self, x, time, event):
        self.inner.fit(x, Surv.from_arrays(event.astype(bool), time.astype(float)))
        return self

    def predict(self, x):
        return np.asarray(self.inner.predict(x), dtype=float)
