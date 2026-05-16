"""Regression test (4): weighted-sum objective must be properly normalised.

Bug: old formula used direct subtraction (Mtot = w_eta*η - w_Q*Q + …),
which could yield negative Mtot values and was not bounded in [0,1].
Fix: pareto_core.hpp uses (1 − Q/Q_max) and (1 − M_abs_err/M_max) so that
every term is in [0,1] and Mtot ∈ [0,1]. Weights are also re-normalised to
sum to 1.0 before use.

The C++ test T4 expected values were stale (0.025/0.140/0.110, matching the
old subtraction formula). They must be updated to 0.525/0.640/0.610.
"""

import subprocess

from conftest import read_source


def test_mtot_formula_uses_one_minus_Q_n():
    src = read_source("analysis/pareto_analysis/pareto_core.hpp")
    assert "1.0 - Q_n" in src, (
        "pareto_core.hpp: Mtot must use (1.0 - Q_n) to keep Q contribution in [0,1]; "
        "direct subtraction was the bug"
    )


def test_mtot_formula_uses_one_minus_M_n():
    src = read_source("analysis/pareto_analysis/pareto_core.hpp")
    assert "1.0 - M_n" in src, (
        "pareto_core.hpp: Mtot must use (1.0 - M_n) for magnification error term"
    )


def test_weights_are_renormalised_by_w_sum():
    src = read_source("analysis/pareto_analysis/pareto_core.hpp")
    assert "w_sum" in src, (
        "pareto_core.hpp: weights must be divided by w_sum to guarantee they sum to 1"
    )


def test_pareto_selector_binary_all_pass(pareto_bin):
    """After T4 fix, all C++ tests must pass — including the Mtot assertions."""
    result = subprocess.run([str(pareto_bin)], capture_output=True, text=True)
    assert result.returncode == 0, (
        f"test_pareto_selector failed.\n"
        f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    )
    # "[FAIL]" marks individual sub-test failures; "0 FAIL" in the summary is expected
    assert "[FAIL]" not in result.stdout, (
        f"One or more C++ sub-tests failed:\n{result.stdout}"
    )
