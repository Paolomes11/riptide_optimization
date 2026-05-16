"""Regression test (10): x_virtual shifting must not be applied as scan_min.

Bug: resolution_map.cpp clamped cfg_scan_min = max(scan_min, x_virtual).
Because x_virtual ≈ 180 mm and the optical focus is often below that, the
entire focus region was excluded from the scan, producing empty results.
Also, the default dof_x_scan_min in config.json was 180 mm, compounding the issue.

Fix:
  - cfg_scan_min = max(scan_min, c.x2 + lens_det_gap)  (physical clearance, not x_virtual)
  - dof_x_scan_min default changed from 180 mm → 100 mm in config.json
"""

import json

from conftest import REPO_ROOT, read_source


def test_cfg_scan_min_not_clamped_to_x_virtual():
    src = read_source("analysis/resolution_analysis/resolution_map.cpp")
    idx = src.index("cfg_scan_min")
    line_end = src.index("\n", idx)
    line = src[idx:line_end]
    assert "x_virtual" not in line, (
        f"resolution_map.cpp: cfg_scan_min line contains 'x_virtual' — "
        f"this is the scan-exclusion bug. Line: {line.strip()}"
    )


def test_cfg_scan_min_uses_lens_det_gap():
    src = read_source("analysis/resolution_analysis/resolution_map.cpp")
    idx = src.index("cfg_scan_min")
    line_end = src.index("\n", idx)
    line = src[idx:line_end]
    assert "lens_det_gap" in line, (
        f"resolution_map.cpp: cfg_scan_min should use 'lens_det_gap' for physical clearance. "
        f"Line: {line.strip()}"
    )


def test_x_virtual_still_used_for_ray_propagation():
    src = read_source("analysis/resolution_analysis/resolution_map.cpp")
    # x_virtual must still appear in the ray-propagation arithmetic (dx = x_focus - x_virtual)
    assert "x_virtual" in src, (
        "resolution_map.cpp: x_virtual disappeared entirely — "
        "it must still be used for ray propagation arithmetic"
    )


def test_dof_x_scan_min_default_below_150mm():
    cfg = json.loads((REPO_ROOT / "config" / "config.json").read_text())
    scan_min = cfg.get("dof_x_scan_min", 180.0)
    assert scan_min < 150.0, (
        f"config/config.json: dof_x_scan_min={scan_min} mm is too large "
        "(old bug used 180 mm which excluded the focal region). Expected ~100 mm."
    )
