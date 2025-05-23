"""
Functions dedicated to the island-finding algorithm, also known as the
Connectedness (CT) Test.

This aims to build a map of adjacent ionisation hits above an ionisation
threshold value in an event. Different sets of connected hits will be labelled
differently, and their positions recorded. Reconstructed particles can be
associated to islands by converting their positions from 3D detector
coordinates to wire-time coordinates for each wire plane.

This algorithm creates a _window_ of wires to collect hits from, this is
centered on the wire-time position of the reconstructed neutrino vertex.
Currently the data structure that represents the hits is a nested list of
floats that each represent the peak time tick value of the hit. An approximate
datashape expression for this is:
```
<window width> * var * float
```

The size of the wire window (number of wires) and the minimum ADC threshold
value are used and applied upstream in the
[LArSoft module](https://github.com/playonverbs/HyperonProduction).

> Unresponsive wire regions have yet to be handled by this code.
"""

from typing import Literal

import awkward as ak
import hist
import numpy as np
import scipy.ndimage as ndi

import sigmazerosearch.utils as utils
from sigmazerosearch.selection import ParameterSet

LABEL_STRUCTURE = [
    [1, 1, 1],
    [1, 1, 1],
    [1, 1, 1],
]
"""Structure passed to <inv:#scipy.ndimage.label> to decide adjacent hits"""


def compute_island_sizes(
    islands, metric: Literal["hits", "size"] = "hits"
) -> list[float] | list[int]:
    """
    Given the output of the <inv:#scipy.ndimage.label>, an option of which
    metric for size to use and a <project:#ParameterSet>.

    This either computes the length of the longest axes of each island or the
    number of hits encompassed by an island.

    Note that when passing 'size', there is an axis mismatch, one axis
    corresponds to drift time ticks while the other corresponds to anode plane
    wires. These are not equal and should be transformed into detector
    coordinates.
    """
    lbl, _ = islands

    if metric == "hits":
        return [lbl[obj].size for obj in ndi.find_objects(lbl)]
    elif metric == "size":
        # Computes the euclidean distance
        # FIXME: Not quite correct, assumes equal time-wire unit lengths.
        return [
            np.sqrt(np.add.reduce(np.square(lbl[obj].shape)))
            for obj in ndi.find_objects(lbl)
        ]


def filter_window_sizes(islands, pset: ParameterSet):
    """
    Given the output of <project:#_window_to_map> or
    <inv:#scipy.ndimage.label>, filter the resulting islands and return a new
    NDArray representing the filtered array and the number of surviving
    islands.
    """
    assert pset.ct_island_size is not None

    lbls, _ = islands

    lengths = ak.Array(compute_island_sizes(islands))
    indices = ak.local_index(lengths)[lengths >= pset.ct_island_size] + 1

    return np.where(np.isin(lbls, ak.to_numpy(indices)), lbls, 0), len(indices)


def convert_seeds_to_window(seeds: dict):
    """
    Takes a map of names to x,y,z positions and converts them to wire,time
    coordinates
    """
    for value in seeds.values():
        utils.WireGeometry.pos_to_u(*value)
    raise NotImplementedError


# TODO: change from using hist.Hist to np.histogram2d
def _window_to_map(window: ak.Array, time_bins: int, wires_max: int):
    """
    Converts a single window object to a labelled map. Returns the labelled map
    and the number of islands found.
    """
    h = (
        hist.Hist.new.Regular(time_bins, 0, 7500, name="time", label="Time [tick]")
        .Integer(0, wires_max, name="wire", label="Relative Wire")
        .Double()
    )

    wires, times = ak.broadcast_arrays(ak.local_index(window, axis=0), window)

    h.fill_flattened(time=times, wire=wires)

    # XXX: think about this unsafe casting
    clipped = np.clip(h.values(), 0, 1, dtype=np.int32, casting="unsafe")

    # NOTE: also consider a function to return the intermediate 'event display'
    return ndi.label(clipped, structure=LABEL_STRUCTURE)


def count_event_islands(arr: ak.Array, pset: ParameterSet) -> ak.Array:
    """
    Takes a sample <inv:#ak.Array> and returns the number of islands found in
    each plane in each event
    """
    filter_sizes = pset.ct_island_size is not None

    islands = {
        "ct_test_islands_plane0": [],
        "ct_test_islands_plane1": [],
        "ct_test_islands_plane2": [],
    }

    # TODO: add option for selecting wire view

    for window in arr.ct_test_window_plane0:
        labelled = _window_to_map(window, pset.ct_time_bins, pset.ct_wire_window)
        if filter_sizes:
            labelled = filter_window_sizes(labelled, pset)
        islands["ct_test_islands_plane0"].append(labelled[1])

    for window in arr.ct_test_window_plane1:
        labelled = _window_to_map(window, pset.ct_time_bins, pset.ct_wire_window)
        if filter_sizes:
            labelled = filter_window_sizes(labelled, pset)
        islands["ct_test_islands_plane1"].append(labelled[1])

    for window in arr.ct_test_window_plane2:
        labelled = _window_to_map(window, pset.ct_time_bins, pset.ct_wire_window)
        if filter_sizes:
            labelled = filter_window_sizes(labelled, pset)
        islands["ct_test_islands_plane2"].append(labelled[1])

    return ak.from_iter(islands)


def make_event_map(arr: ak.Array, pset: ParameterSet) -> ak.Array:
    """
    Creates event maps per event when given an awkward array that
    """

    raise NotImplementedError
