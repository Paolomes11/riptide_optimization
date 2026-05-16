"""Shared fixtures and helpers for RIPTIDE regression tests."""

from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).parent.parent
BIN_DIR = REPO_ROOT / "build" / "Release"
ANALYSIS_BIN_DIR = REPO_ROOT / "build" / "analysis"


def read_source(rel_path: str) -> str:
    """Read a source file relative to the repo root."""
    return (REPO_ROOT / rel_path).read_text()


@pytest.fixture
def pareto_bin():
    p = ANALYSIS_BIN_DIR / "pareto_analysis" / "Release" / "test_pareto_selector"
    if not p.exists():
        pytest.skip(f"Binary not found: {p}")
    return p
