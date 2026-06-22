#!/usr/bin/env python3
"""One-time migration: re-pickle model artifacts so they load cold.

The original ``models/*.pkl`` artifacts were written while
``scripts/train_ml_zoo.py`` ran as ``__main__``, so the embedded wrapper
classes (``SkSurvAdapter`` / ``LifelinesCoxWrapper`` / ``XGBCoxWrapper``) were
recorded with module ``__main__`` and could not be ``pickle.load``-ed from a
fresh process without an import shim.

This script loads each artifact (installing a temporary ``__main__`` alias so
legacy pickles resolve), then re-dumps it. Because the loaded objects' classes
now live in ``src.ml.wrappers``, the re-dumped pickles reference
``src.ml.wrappers.*`` and load cold via ``scripts/predict.py`` /
``scripts/external_validation.py`` with no shim.

The script is idempotent: running it on already-migrated artifacts simply
re-dumps identical objects. Originals are backed up under
``models/_pre_migration/`` (gitignored) on first run.

Usage:
    python scripts/migrate_pickles.py
"""

from __future__ import annotations

import pickle
import shutil
import sys
import __main__
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from src.ml.wrappers import (  # noqa: E402
    LifelinesCoxWrapper,
    SkSurvAdapter,
    XGBCoxWrapper,
)

MODELS = REPO_ROOT / "models"
BACKUP = MODELS / "_pre_migration"

# Pickle artifacts that may embed the legacy __main__-scoped wrapper classes.
PICKLES = [
    "cox_ph.pkl",
    "elastic_net_cox.pkl",
    "random_survival_forest.pkl",
    "gradient_boosted_survival.pkl",
    "stacked_ensemble.pkl",
]


def install_legacy_aliases() -> None:
    """Make legacy ``__main__.<Wrapper>`` references resolve to src.ml.wrappers."""
    __main__.SkSurvAdapter = SkSurvAdapter
    __main__.LifelinesCoxWrapper = LifelinesCoxWrapper
    __main__.XGBCoxWrapper = XGBCoxWrapper


def migrate(name: str) -> None:
    path = MODELS / name
    if not path.exists():
        print(f"  skip (absent): {name}")
        return
    BACKUP.mkdir(exist_ok=True)
    backup_path = BACKUP / name
    if not backup_path.exists():
        shutil.copy2(path, backup_path)
    with path.open("rb") as fh:
        artifact = pickle.load(fh)
    with path.open("wb") as fh:
        pickle.dump(artifact, fh)
    print(f"  migrated: {name}")


def verify_cold(name: str) -> bool:
    """Load in a fresh subprocess with no __main__ shim to prove cold-loadability."""
    import subprocess

    code = (
        "import sys, pickle; sys.path.insert(0, r'%s'); "
        "pickle.load(open(r'%s','rb')); print('ok')"
        % (str(REPO_ROOT), str(MODELS / name))
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    return result.returncode == 0 and "ok" in result.stdout


def main() -> None:
    print("Migrating model pickles to src.ml.wrappers ...")
    install_legacy_aliases()
    for name in PICKLES:
        migrate(name)
    print("Verifying cold load (fresh process, no shim) ...")
    all_ok = True
    for name in PICKLES:
        if not (MODELS / name).exists():
            continue
        ok = verify_cold(name)
        print(f"  {'PASS' if ok else 'FAIL'}: {name}")
        all_ok = all_ok and ok
    if not all_ok:
        raise SystemExit("Cold-load verification failed for at least one artifact")
    print("All artifacts load cold without an import shim.")


if __name__ == "__main__":
    main()
