"""Regression test (5): lens_scan and optimizer must use the same gap margin.

Bug: the two components used different hardcoded margin values (e.g. 1.0 mm vs
3.0 mm), causing the lens pair grids to be inconsistent: the optimizer would
try lens positions that the scanner had not simulated.
Fix: both read from the same config key "lens_gap_margin" with the same default.
"""

import json
import re

from conftest import REPO_ROOT, read_source


def _extract_margin_default(src: str) -> str | None:
    m = re.search(r'lens_gap_margin",\s*([\d.]+)', src)
    return m.group(1) if m else None


def test_lens_scan_and_optimizer_use_same_margin_default():
    ls_src = read_source("src/lens_simulation/lens_scan.cpp")
    op_src = read_source("src/optimization/optimizer.cpp")
    ls_default = _extract_margin_default(ls_src)
    op_default = _extract_margin_default(op_src)
    assert ls_default is not None, "lens_scan.cpp: lens_gap_margin default not found"
    assert op_default is not None, "optimizer.cpp: lens_gap_margin default not found"
    assert ls_default == op_default, (
        f"Margin defaults differ: lens_scan={ls_default}mm, optimizer={op_default}mm. "
        "They must be equal to keep the pair grids consistent."
    )


def test_config_json_has_lens_gap_margin():
    cfg = json.loads((REPO_ROOT / "config" / "config.json").read_text())
    assert "lens_gap_margin" in cfg, (
        "config/config.json must declare 'lens_gap_margin' so both components "
        "read from the same authoritative source"
    )


def test_both_use_config_value_not_hardcoded():
    for path in ["src/lens_simulation/lens_scan.cpp", "src/optimization/optimizer.cpp"]:
        src = read_source(path)
        assert 'config.value("lens_gap_margin"' in src, (
            f"{path}: margin must come from config.value(\"lens_gap_margin\"), not a hardcoded literal"
        )
