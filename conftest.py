"""Pytest configuration: deterministic imports + automatic R/genefu skip.

This file does two things:

1. Insert the repository root on ``sys.path`` so ``from src...`` and
   ``from scripts...`` imports resolve regardless of where pytest is invoked
   from (the test modules also do this themselves; doing it here makes
   collection itself robust).

2. Detect whether an R installation with the ``genefu`` Bioconductor package is
   available. The official PAM50-ROR test (``tests/test_pam50_official.py``)
   shells out to ``Rscript`` via ``scripts.compute_clinical_baselines``; when R
   or genefu is missing it is skipped cleanly instead of erroring at import.
"""

from __future__ import annotations

import functools
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# Tests that require R/genefu, identified by module filename.
_R_DEPENDENT_MODULES = {"test_pam50_official.py"}


@functools.lru_cache(maxsize=1)
def _r_with_genefu_available() -> bool:
    """Return True iff Rscript exists and `library(genefu)` loads."""
    if shutil.which("Rscript") is None:
        return False
    try:
        completed = subprocess.run(
            ["Rscript", "-e", "suppressMessages(library(genefu))"],
            capture_output=True,
            text=True,
            timeout=120,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return False
    return completed.returncode == 0


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line(
        "markers",
        "requires_r: test requires an R installation with the genefu package "
        "(PAM50-ROR); auto-skipped when absent",
    )


def pytest_collection_modifyitems(
    config: pytest.Config, items: list[pytest.Item]
) -> None:
    """Mark R-dependent tests and skip them when R/genefu is unavailable."""
    if _r_with_genefu_available():
        # Still tag the items so `-m requires_r` selection works, but do not skip.
        for item in items:
            if Path(str(item.fspath)).name in _R_DEPENDENT_MODULES:
                item.add_marker(pytest.mark.requires_r)
        return

    skip_no_r = pytest.mark.skip(
        reason="requires R with the genefu package (PAM50-ROR); not available"
    )
    for item in items:
        if Path(str(item.fspath)).name in _R_DEPENDENT_MODULES:
            item.add_marker(pytest.mark.requires_r)
            item.add_marker(skip_no_r)
