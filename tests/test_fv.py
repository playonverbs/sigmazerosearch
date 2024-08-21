import pytest

from sigmazerosearch.alg.fv import in_active_tpc


@pytest.mark.parametrize(
    "x, y, z, want",
    [
        (1, 2, 3, True),
        (-999, -999, -999, False),
        # from Claude
        (126.625, 0.97, 518.5, True),
        (0.0, 0.0, 0.0, False),
        (236.35, 233, 1036.8, False),
        (-1, -1, -1, False),
        (256.36, 0, 0, False),
        (0, 117.48, 0, False),
        (0, 0, 1036.91, False),
        (126.625, 0.97, 0.1, True),
    ],
)
def test_in_active_tpc(x, y, z, want):
    assert in_active_tpc(x, y, z) == want
