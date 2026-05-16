"""Regression test (6): n_hits must be stored and read with a consistent type.

Bug: mismatch between writer type (Int_t) and reader type (double) for the
"n_hits" ROOT branch caused silent data corruption or read errors.
Fix: all writers use CreateNtupleDColumn("n_hits") and all readers bind to
a `double` variable via SetBranchAddress.
"""

from conftest import read_source

WRITER_FILES = [
    "src/lens_simulation/lens_scan.cpp",
    "src/optimization/optimizer.cpp",
    "src/psf_dof_scan/psf_dof_scan.cpp",
]

READER_FILES = [
    "analysis/plot2D.cpp",
    "analysis/resolution_analysis/resolution_map.cpp",
    "analysis/psf_extractor.cpp",
]


def test_nhits_created_as_double_not_int():
    for path in WRITER_FILES:
        src = read_source(path)
        assert 'CreateNtupleDColumn("n_hits"' in src, (
            f"{path}: n_hits must be created with CreateNtupleDColumn (double), not I/F"
        )
        assert 'CreateNtupleIColumn("n_hits"' not in src, (
            f"{path}: found CreateNtupleIColumn(\"n_hits\") — this is the type-mismatch bug"
        )


def test_nhits_read_with_setbranchaddress():
    for path in READER_FILES:
        src = read_source(path)
        assert 'SetBranchAddress("n_hits"' in src, (
            f"{path}: expected SetBranchAddress(\"n_hits\", ...) binding"
        )
