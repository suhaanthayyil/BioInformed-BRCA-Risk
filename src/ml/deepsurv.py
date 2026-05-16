"""Small PyTorch DeepSurv implementation for tabular Cox modeling."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch
from torch import nn


RANDOM_SEED = 42


class CoxMLP(nn.Module):
    """MLP that outputs a log-risk score for Cox partial likelihood."""

    def __init__(self, in_features: int, hidden_dims: tuple[int, ...] = (64, 32), dropout: float = 0.1):
        super().__init__()
        layers: list[nn.Module] = []
        prev = in_features
        for dim in hidden_dims:
            layers.extend([nn.Linear(prev, dim), nn.BatchNorm1d(dim), nn.ReLU(), nn.Dropout(dropout)])
            prev = dim
        layers.append(nn.Linear(prev, 1))
        self.network = nn.Sequential(*layers)

    def forward(self, x):
        return self.network(x).squeeze(-1)


def cox_partial_likelihood_loss(risk: torch.Tensor, time: torch.Tensor, event: torch.Tensor) -> torch.Tensor:
    """Negative Cox partial log likelihood for right-censored data."""

    order = torch.argsort(time, descending=True)
    risk = risk[order]
    event = event[order]
    log_cumsum_hazard = torch.logcumsumexp(risk, dim=0)
    observed = event > 0
    if observed.sum() == 0:
        return torch.tensor(0.0, dtype=risk.dtype, device=risk.device, requires_grad=True)
    return -torch.sum((risk - log_cumsum_hazard)[observed]) / observed.sum()


@dataclass
class DeepSurvConfig:
    hidden_dims: tuple[int, ...] = (64, 32)
    dropout: float = 0.1
    lr: float = 1e-3
    weight_decay: float = 1e-4
    epochs: int = 150
    patience: int = 20


class DeepSurvEstimator:
    """Minimal sklearn-like estimator around ``CoxMLP``."""

    def __init__(self, config: DeepSurvConfig | None = None, device: str | None = None):
        self.config = config or DeepSurvConfig()
        self.device = device or ("mps" if torch.backends.mps.is_available() else "cpu")
        self.model_: CoxMLP | None = None

    def fit(self, x, time, event, x_val=None, time_val=None, event_val=None):
        torch.manual_seed(RANDOM_SEED)
        np.random.seed(RANDOM_SEED)
        x_arr = np.asarray(x, dtype=np.float32)
        self.model_ = CoxMLP(x_arr.shape[1], self.config.hidden_dims, self.config.dropout).to(self.device)
        optimizer = torch.optim.AdamW(
            self.model_.parameters(),
            lr=self.config.lr,
            weight_decay=self.config.weight_decay,
        )
        x_t = torch.tensor(x_arr, dtype=torch.float32, device=self.device)
        time_t = torch.tensor(np.asarray(time, dtype=np.float32), dtype=torch.float32, device=self.device)
        event_t = torch.tensor(np.asarray(event, dtype=np.float32), dtype=torch.float32, device=self.device)

        if x_val is not None:
            xv_t = torch.tensor(np.asarray(x_val, dtype=np.float32), dtype=torch.float32, device=self.device)
            tv_t = torch.tensor(np.asarray(time_val, dtype=np.float32), dtype=torch.float32, device=self.device)
            ev_t = torch.tensor(np.asarray(event_val, dtype=np.float32), dtype=torch.float32, device=self.device)
        else:
            xv_t = tv_t = ev_t = None

        best_state = None
        best_loss = float("inf")
        stale = 0
        for _ in range(self.config.epochs):
            self.model_.train()
            optimizer.zero_grad(set_to_none=True)
            risk = self.model_(x_t)
            loss = cox_partial_likelihood_loss(risk, time_t, event_t)
            loss.backward()
            optimizer.step()

            self.model_.eval()
            with torch.no_grad():
                if xv_t is not None and tv_t is not None and ev_t is not None:
                    val_loss = cox_partial_likelihood_loss(self.model_(xv_t), tv_t, ev_t).item()
                else:
                    val_loss = loss.item()
            if val_loss < best_loss:
                best_loss = val_loss
                best_state = {k: v.detach().cpu().clone() for k, v in self.model_.state_dict().items()}
                stale = 0
            else:
                stale += 1
            if stale >= self.config.patience:
                break

        if best_state is not None:
            self.model_.load_state_dict(best_state)
        return self

    def predict(self, x) -> np.ndarray:
        if self.model_ is None:
            raise RuntimeError("DeepSurvEstimator is not fitted")
        self.model_.eval()
        with torch.no_grad():
            x_t = torch.tensor(np.asarray(x, dtype=np.float32), dtype=torch.float32, device=self.device)
            risk = self.model_(x_t).detach().cpu().numpy()
        return risk.astype(float)

    def state_dict(self) -> dict:
        if self.model_ is None:
            raise RuntimeError("DeepSurvEstimator is not fitted")
        return self.model_.state_dict()
