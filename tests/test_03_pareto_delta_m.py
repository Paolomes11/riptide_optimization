"""Regression test (3): Pareto selector must include ΔM as a dominance dimension.

Bug: pareto_dominates() only compared η and Q, ignoring magnification error ΔM.
A config with good η/Q but large ΔM would incorrectly appear on the Pareto front.
Fix: dominance requires A.M_abs_err <= B.M_abs_err (with at least one strict inequality).
"""

import subprocess

from conftest import read_source


def test_pareto_dominates_includes_M_abs_err_in_signature():
    src = read_source("analysis/pareto_analysis/pareto_core.hpp")
    # Find the body of pareto_dominates and verify M_abs_err appears there
    fn_start = src.index("pareto_dominates")
    fn_body = src[fn_start : fn_start + 400]
    assert "M_abs_err" in fn_body, (
        "pareto_core.hpp: pareto_dominates() must reference M_abs_err "
        "so that ΔM is a dominance dimension"
    )


def test_pareto_dominates_all_three_conditions_strict():
    src = read_source("analysis/pareto_analysis/pareto_core.hpp")
    fn_start = src.index("pareto_dominates")
    fn_body = src[fn_start : fn_start + 400]
    # The strict inequality on M_abs_err must be present
    assert "M_abs_err < B.M_abs_err" in fn_body or "A.M_abs_err < B" in fn_body, (
        "pareto_core.hpp: strict M_abs_err inequality missing from pareto_dominates()"
    )


def test_pareto_selector_binary_passes(pareto_bin):
    """The compiled C++ test suite must exit 0 (covers T1–T7 including ΔM cases)."""
    result = subprocess.run([str(pareto_bin)], capture_output=True, text=True)
    assert result.returncode == 0, (
        f"test_pareto_selector exited {result.returncode}.\n"
        f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    )
