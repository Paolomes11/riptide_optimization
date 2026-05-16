"""Regression test (7): --focus-tsv must generate pairs by reading the TSV file.

Bug: the --focus-tsv flag was accepted but the pair generation still fell
through to the geometric loop, ignoring the TSV contents and potentially
simulating positions where no focus data existed.
Fix: when use_focus_map is true, pairs come exclusively from
get_pairs_from_focus_map(); the geometric loop is the else-branch fallback.
"""

from conftest import read_source


def test_focus_tsv_branch_reads_from_get_pairs_from_focus_map():
    src = read_source("src/lens_simulation/lens_scan.cpp")
    assert "use_focus_map" in src, (
        "lens_scan.cpp: use_focus_map flag missing — TSV-based pair generation not implemented"
    )
    assert "get_pairs_from_focus_map" in src, (
        "lens_scan.cpp: get_pairs_from_focus_map() call missing in TSV branch"
    )


def test_geometric_loop_is_else_fallback():
    src = read_source("src/lens_simulation/lens_scan.cpp")
    idx_tsv = src.index("get_pairs_from_focus_map")
    # The geometric loop (identified by lens_gap_margin) must appear AFTER an
    # else-branch that follows the TSV block, not before it.
    remaining = src[idx_tsv:]
    assert "} else {" in remaining, (
        "lens_scan.cpp: geometric loop is not an else-fallback after the TSV branch"
    )
    idx_else = remaining.index("} else {")
    idx_margin = remaining.index("lens_gap_margin", idx_else)
    assert idx_margin > idx_else, (
        "lens_scan.cpp: lens_gap_margin loop appears before the else-branch — "
        "geometric loop is not a proper fallback"
    )


def test_focus_map_hpp_extracts_pairs_from_map():
    src = read_source("include/common/focus_map.hpp")
    assert "get_pairs_from_focus_map" in src, (
        "focus_map.hpp: get_pairs_from_focus_map() helper not defined"
    )
