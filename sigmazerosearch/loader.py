"""
Abstracted data-file loading procedures.

Currently ROOT NTuple files are accepted (as per the output of
HyperonProduction). However other data files may be added in the future.
"""

from os.path import isabs
from pathlib import Path

import awkward as ak
import numpy as np
import uproot as up
from numpy.typing import ArrayLike
from uproot.behaviors.TBranch import HasBranches

from sigmazerosearch.general import Config
from sigmazerosearch.utils import WireGeometry


def _yield_array_from_ttree(tree: HasBranches, config: Config):
    """
    Iterate through a given TTree object in chunks, yielding the current
    chunked array.

    The size of the chunks and the branches read from the TTree is determined
    by the selection config.
    """
    # if config.iterate:
    for arr in tree.iterate(
        expressions=config.branch_list, step_size=config.iterate_step, report=None
    ):  # type: ignore
        yield arr
    # else:
    #     return tree.arrays(config.branch_list)


def load_ntuple(filename: str) -> HasBranches:
    """
    Wraps the uproot.open method, taking a filename and outputting some ROOT
    object that has branches
    """
    data = up.open(filename)
    if not isabs(filename):
        raise OSError("Please provide an absolute file path")
    if isinstance(data, HasBranches):
        return data
    else:
        raise TypeError("please read a TTree using the ':' separator.")


def get_POT(filename: str) -> float:
    """Sums the POT value of each subrun from a given ROOT file"""
    if not isabs(filename):
        raise OSError("Please provide an absolute file path")
    data = up.open(filename)
    pots_per_subrun: ak.Array = data.get("ana/MetaTree").get("POT").array()  # type: ignore
    return ak.sum(pots_per_subrun)


def load_dead_wire_map(filename: str | Path, transpose: bool = True) -> ArrayLike:
    """
    Given a path to a dead-wire channel map, produce an array of all valid
    wires. Optionally produced for each plane.

    The dead-wire channel map is expected to be a text file with the format:
    ```
    {wire-number} {status}
    ...
    ```

    - `{wire-number}` is an integer that ranges from 0 to the total number of
    wires - 1 in the detector across all planes (`0..8255` for MicroBooNE).
    - `{status}` is an integer with 1 representing a dead channel and 4
    representing a live channel.

    :param str filename: The filename containing the dead-wire map.
    :param bool transpose:
        Whether to transpose the wire numbers to represent active wires using a
        per-plane ID.

    :return:
        Either a numpy array with all dead wire channel numbers; or a tuple of
        arrays with dead wire channel numbers for each plane
    """
    data = np.loadtxt(filename, dtype=np.int16, unpack=True)

    if not transpose:
        wires, status = data
        return wires[status != 4]
    else:
        bounds = WireGeometry.wire_boundaries()
        return (
            data[..., bounds[0] : bounds[1]][0][
                data[..., bounds[0] : bounds[1]][1] != 4
            ]
            - bounds[0],
            data[..., bounds[1] : bounds[2]][0][
                data[..., bounds[1] : bounds[2]][1] != 4
            ]
            - bounds[1],
            data[..., bounds[2] : bounds[3]][0][
                data[..., bounds[2] : bounds[3]][1] != 4
            ]
            - bounds[2],
        )
