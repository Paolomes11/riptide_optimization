"""Regression test (8): detector position must not be hardcoded at 180 mm.

Bug: the simulation's detector (photocathode) position was fixed at 180 mm
regardless of the config, making all efficiency/PSF scans use the wrong
detector plane.
Fix: the simulation reads x_det from config.json (default 579.8 mm) and
SetDetectorPosition() applies it each scan step.
"""

import json

from conftest import REPO_ROOT, read_source


def test_config_json_xdet_not_180():
    cfg = json.loads((REPO_ROOT / "config" / "config.json").read_text())
    x_det = cfg.get("x_det")
    assert x_det is not None, "config/config.json must contain 'x_det'"
    assert x_det != 180.0, (
        f"config/config.json: x_det={x_det} — detector position is 180 mm (the old hardcoded bug value)"
    )


def test_optimizer_reads_xdet_from_config():
    src = read_source("src/optimization/optimizer.cpp")
    assert 'config.value("x_det"' in src or '"x_det"' in src, (
        "optimizer.cpp: x_det must be read from config, not hardcoded"
    )


def test_optimizer_does_not_hardcode_xdet_180():
    src = read_source("src/optimization/optimizer.cpp")
    # A literal assignment `x_det = 180` in simulation code is the bug
    assert "x_det = 180" not in src, (
        "optimizer.cpp: found 'x_det = 180' — hardcoded detector position bug"
    )


def test_set_detector_position_called_in_optimizer():
    src = read_source("src/optimization/optimizer.cpp")
    assert "SetDetectorPosition" in src, (
        "optimizer.cpp: SetDetectorPosition() must be called to apply the configured x_det"
    )
