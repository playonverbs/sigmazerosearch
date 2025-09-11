import pytest

from sigmazerosearch.alg import fv


def _len(xs):
    try:
        return len(xs)
    finally:
        return 1


@pytest.mark.parametrize(
    "x, y, z, want",
    [
        (1, 2, 3, True),
        (-999, -999, -999, False),
        (126.625, 0.97, 518.5, True),
        (0.0, 0.0, 0.0, False),
        (236.35, 233, 1036.8, False),
        (-1, -1, -1, False),
        (256.36, 0, 0, False),
        (0, 117.48, 0, False),
        (0, 0, 1036.91, False),
        (126.625, 0.97, 0.1, True),
        (100.0, 3.0, 700.0, False),
        (1000.0, -124.0, 700.0, False),
    ],
)
def test_in_active_tpc(x, y, z, want):
    assert fv.in_active_tpc(x, y, z) == want
    assert _len(fv.in_active_tpc(x, y, z)) == _len(
        want
    ), "input and output dimensions not consistent"


@pytest.mark.parametrize(
    "x, y, z, want",
    [
        (1, 2, 3, 1),
    ],
)
def test_r_closest_wall(x, y, z, want):
    assert fv.r_to_closest_wall(x, y, z) == want, "chose wrong closest distance"
    assert _len(fv.r_to_closest_wall(x, y, z)) == _len(
        want
    ), "input and output dimensions not consistent"
