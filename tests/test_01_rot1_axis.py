"""Regression test (1): rot1 axis correctness in DOF source geometry.

Bug: wrong rot1 axis in GPS commands rotated the source plane off-axis,
misaligning the photon beam with the optical axis.
Fix: /gps/pos/rot1 must be X=(1,0,0); /gps/ang/rot1 must be -Y=(0,-1,0).
"""

from conftest import read_source

# dof_scan.cpp uses a rectangular Plane source → needs pos/rot axes.
# psf_dof_scan.cpp uses a Point source → no pos/rot1 needed (different geometry).
DOF_PLANE_SRC = "src/dof_simulation/dof_scan.cpp"

# Both plane-source and point-source sims share the same angular distribution axes.
DOF_ANG_SRCS = [
    "src/dof_simulation/dof_scan.cpp",
    "src/psf_dof_scan/psf_dof_scan.cpp",
]


def test_pos_rot1_is_x_axis():
    src = read_source(DOF_PLANE_SRC)
    assert '"/gps/pos/rot1 1 0 0"' in src, (
        f"{DOF_PLANE_SRC}: /gps/pos/rot1 must be (1 0 0) — X axis keeps source plane horizontal"
    )


def test_pos_rot2_is_y_axis():
    src = read_source(DOF_PLANE_SRC)
    assert '"/gps/pos/rot2 0 1 0"' in src, (
        f"{DOF_PLANE_SRC}: /gps/pos/rot2 must be (0 1 0)"
    )


def test_ang_rot1_is_neg_y_axis():
    for f in DOF_ANG_SRCS:
        src = read_source(f)
        assert '"/gps/ang/rot1 0 -1 0"' in src, (
            f"{f}: /gps/ang/rot1 must be (0 -1 0) — neg-Y points angular distribution toward detector"
        )


def test_ang_rot2_is_z_axis():
    for f in DOF_ANG_SRCS:
        src = read_source(f)
        assert '"/gps/ang/rot2 0 0 1"' in src, (
            f"{f}: /gps/ang/rot2 must be (0 0 1)"
        )
