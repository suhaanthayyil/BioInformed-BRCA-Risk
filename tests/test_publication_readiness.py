from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "scripts" / "check_publication_readiness.py"
MANUSCRIPT = ROOT / "_submission" / "Manuscript_prepublication_revision_v1.2.1.docx"
SUPPLEMENT = ROOT / "_submission" / "BMC_Supplement_prepublication_revision_v1.2.1.docx"
RESPONSE = ROOT / "_submission" / "Response_to_prepublication_checks_v1.2.1.docx"


def run_audit(*arguments: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(AUDIT), *arguments],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )


def test_repository_audit_proves_publication_claims() -> None:
    result = run_audit()

    assert result.returncode == 0, result.stdout + result.stderr
    assert "PASS version: 1.2.1" in result.stdout
    assert "PASS cohorts: analysed=6644; PAM50-evaluable=6643" in result.stdout
    assert "PASS registered endpoints: 8 rows; Oncotype DX dispositions=2" in result.stdout
    assert "PASS registration deviations: 16 rows" in result.stdout
    assert "PASS stage availability:" in result.stdout
    assert "TCGA-BRCA=1042" in result.stdout
    assert "METABRIC=1764" in result.stdout
    assert "GSE20685=325" in result.stdout
    assert "GSE96058=0" in result.stdout


@pytest.mark.skipif(
    not all(path.exists() for path in (MANUSCRIPT, SUPPLEMENT, RESPONSE)),
    reason="private journal-submission DOCX files are not tracked in the code repository",
)
def test_document_audit_proves_submission_package_is_publication_ready() -> None:
    result = run_audit(
        "--manuscript",
        str(MANUSCRIPT),
        "--supplement",
        str(SUPPLEMENT),
        "--response",
        str(RESPONSE),
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert "PASS manuscript: chronology explicit; ChatGPT limited to language checks" in result.stdout
    assert "PASS supplement: Tables S18-S21 embedded; Oncotype DX dispositions=2" in result.stdout
    assert "PASS response: versioned files listed; black text on white background" in result.stdout
    assert "PASS clean DOCX: no tracked changes, comments, or prohibited AI attribution" in result.stdout
