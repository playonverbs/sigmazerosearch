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


def convert_seeds_to_window(seeds: dict):
    """
    Takes a map of names to x,y,z positions and converts them to wire,time
    coordinates
    """
    for key, value in seeds.items():
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


def count_event_islands(
    arr: ak.Array, pset: ParameterSet, view: utils.WireView | None = None
) -> ak.Array:
    """
    Takes a sample <inv:#ak.Array> and returns the number of islands found in
    each plane in each event
    """

    islands = {
        "ct_test_islands_plane0": [],
        "ct_test_islands_plane1": [],
        "ct_test_islands_plane2": [],
    }

    # TODO: add option for selecting wire view

    for window in arr.ct_test_window_plane0:
        labelled = _window_to_map(window, pset.ct_time_bins, pset.ct_wire_window)
        islands["ct_test_islands_plane0"].append(labelled[1])

    for window in arr.ct_test_window_plane1:
        labelled = _window_to_map(window, pset.ct_time_bins, pset.ct_wire_window)
        islands["ct_test_islands_plane1"].append(labelled[1])

    for window in arr.ct_test_window_plane2:
        labelled = _window_to_map(window, pset.ct_time_bins, pset.ct_wire_window)
        islands["ct_test_islands_plane2"].append(labelled[1])

    return ak.from_iter(islands)


def make_event_map(arr: ak.Array, pset: ParameterSet) -> ak.Array:
    """
    Creates event maps per event when given an awkward array that
    """

    raise NotImplementedError
