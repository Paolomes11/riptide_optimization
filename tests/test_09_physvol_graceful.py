"""Regression test (9): missing photocathode physvol must be handled gracefully.

Bug: detector_construction.cpp called methods on m_photocathode_phys without
checking for null, causing a segfault when the GDML lacked the physvol.
Fix: an explicit null-check logs a spdlog::warn and returns early, so the
simulation continues without a sensitive detector rather than crashing.
"""

from conftest import read_source

DETECTOR_SRCS = [
    "src/optimization/detector_construction.cpp",
    "src/lens_simulation/detector_construction.cpp",
]


def test_null_check_before_physvol_use():
    for path in DETECTOR_SRCS:
        src = read_source(path)
        assert "if (!m_photocathode_phys)" in src, (
            f"{path}: null-check 'if (!m_photocathode_phys)' missing — "
            "missing physvol would cause a segfault"
        )


def test_warning_logged_when_physvol_missing():
    for path in DETECTOR_SRCS:
        src = read_source(path)
        assert "Photocathode physical volume not found" in src, (
            f"{path}: spdlog warning for missing physvol not found — "
            "silent failure makes debugging impossible"
        )


def test_set_detector_position_also_guards_null():
    for path in DETECTOR_SRCS:
        src = read_source(path)
        # SetDetectorPosition must also guard against null before using the pointer
        fn_idx = src.find("SetDetectorPosition")
        assert fn_idx != -1, f"{path}: SetDetectorPosition not found"
        fn_body = src[fn_idx : fn_idx + 300]
        assert "m_photocathode_phys" in fn_body, (
            f"{path}: SetDetectorPosition does not reference m_photocathode_phys"
        )
        assert "if (!m_photocathode_phys)" in fn_body or "if (!m_photocathode_phys)" in src[fn_idx:fn_idx+500], (
            f"{path}: SetDetectorPosition lacks null-guard for m_photocathode_phys"
        )
