import awkward as ak

from sigmazerosearch.selection import ParameterSet
from sigmazerosearch.utils import displacement


def select_mu_candidate(arr: ak.Array, pset: ParameterSet) -> ak.Array:
    """
    Apply conditions for a muon ID cut using parameters.

    Parameters used:

    1. Muon PID cut: <project:#ParameterSet.pid_cut>
    2. Minimum length: <project:#ParameterSet.min_length>
    3. Max separation: <project:#ParameterSet.max_separation>
    """
    mask = (
        (arr["trk_llrpid"] < pset.pid_cut)
        & (arr["trk_length"] < pset.min_length)
        & ak.fill_none(
            ak.any(
                displacement(arr, "trk_start_x", "trk_start_y", "trk_start_z")
                < pset.max_separation,
                axis=1,
            ),
            False,
        )
    )
    return mask
