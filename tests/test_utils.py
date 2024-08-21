import awkward as ak
import pytest

from sigmazerosearch.utils import filter_by_rse


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
