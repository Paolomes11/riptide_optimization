"""Regression test (2): EE80 must be computed as diameter, not radius.

Bug: EE80 was returned as the raw 80th-percentile radius R80.
Fix: dof_map.cpp multiplies by 2.0 (EE80 = 2*R80); resolution_map.cpp
uses k_EE80 = 2.0 * 1.7941 (diameter factor for Gaussian σ → EE80 diameter).
"""

from conftest import read_source


def test_dof_map_ee80_multiplied_by_two():
    src = read_source("analysis/dof_analysis/dof_map.cpp")
    assert "2.0 * R80" in src, (
        "dof_map.cpp: EE80 must be 2.0 * R80 (diameter), not the raw radius R80"
    )


def test_dof_map_ee80_not_stored_as_raw_r80():
    src = read_source("analysis/dof_analysis/dof_map.cpp")
    # Ensure there is no assignment of raw R80 to out.EE80 without the 2.0 factor
    for line in src.splitlines():
        if "out.EE80" in line and "R80" in line:
            assert "2.0 * R80" in line, (
                f"dof_map.cpp: suspicious EE80 assignment without diameter factor: {line.strip()}"
            )


def test_resolution_map_k_ee80_is_diameter_factor():
    src = read_source("analysis/resolution_analysis/resolution_map.cpp")
    assert "k_EE80 = 2.0 *" in src, (
        "resolution_map.cpp: k_EE80 must include the 2.0 diameter factor"
    )
