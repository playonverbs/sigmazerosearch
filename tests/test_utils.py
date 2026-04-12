import awkward as ak
import numpy as np
import pytest

from sigmazerosearch.alg import fv
from sigmazerosearch.utils import WireGeometry, WireView, filter_by_rse


@pytest.fixture
def example_array():
    return ak.Array(
        {
            "run": [1, 2, 3],
            "subrun": [3213, 321, 4932],
            "event": [12310, 33842, 24352],
        }
    )


def test_filter_by_rse(example_array):
    out = filter_by_rse(example_array, 3, 4932, 24352)

    assert "run" in out.fields
    assert "subrun" in out.fields
    assert "event" in out.fields


def test_wireview():
    for i in (0, 1, 2):
        wv = WireView(i)
        assert WireView.from_suffix(wv.to_suffix()) == wv, "failed roundtripping"

    assert WireView.from_suffix("ct_islands_plane0") == WireView.U
    assert WireView.from_suffix("ct_islands_plane1") == WireView.V
    assert WireView.from_suffix("ct_islands_plane2") == WireView.Y

    for val in ("hello_there_uhh", "hellothere", "hi_plane_"):
        with pytest.raises(ValueError):
            WireView.from_suffix(val)


@pytest.mark.skip(reason="Currently unused conversion")
def test_wiregeometry_conversions():
    """
    Generate a random bunch of x,y,z detector coordinates and convert them to
    plane coordinates
    """
    n_points = 1_000_000
    rng = np.random.default_rng()

    xs = rng.uniform(*fv.FV_x, size=n_points)
    ys = rng.uniform(*fv.FV_y, size=n_points)
    zs = rng.uniform(*fv.FV_z, size=n_points)

    wire_us = WireGeometry.pos_to_u(xs, ys, zs)
    wire_vs = WireGeometry.pos_to_v(xs, ys, zs)
    wire_ys = WireGeometry.pos_to_y(xs, ys, zs)
    time_ticks = WireGeometry.pos_to_time(xs, ys, zs)

    assert np.all(wire_us <= np.min(wire_vs)) and np.all(
        wire_vs <= np.min(wire_ys)
    ), "wire numbers are not sequential in order of U V Y"

    assert np.all(time_ticks >= 0), "negative time ticks are present"
