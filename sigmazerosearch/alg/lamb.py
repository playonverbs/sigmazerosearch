"""
Functions focussing on selecting and capturing information regarding the lambda
decay subsystem.

The targeted decay channel is that of a neutral lambda decaying to a proton and
negative pion: {math}`\\Lambda \\rightarrow p + \\pi^-`.

> This sub-package must be named `lamb` as `lambda` is a reserved Python keyword
"""

import awkward as ak

import sigmazerosearch.utils as utils
from sigmazerosearch.selection import ParameterSet


def select_p_pi_candidates(arr: ak.Array, pset: ParameterSet) -> ak.Array:
    """
    This aims to recreate the proton-pion candidate selection algorithm used in
    the Lambda analysis using a Boosted-Decision Tree method.

    :return:
        An <inv:#ak.Array> with shape equal to [`number of events in the
        analysis`][`2`]. The inner values mark the PFP indices of the proton
        and pion candidates, respectively.
    """

    raise NotImplementedError


def select_p_pi_candidates_box(arr: ak.Array, pset: ParameterSet) -> ak.Array:
    """
    A stop-gap method of selecting the decay products from Lambda decay.

    This uses simpler box cuts to select the appropriate tracks in the event
    using the main three criteria that would be fed into the nominal BDT-based
    selection cut:

    1. Proton hypothesis PID score: <project:#ParameterSet.proton_pid_cut>.
    2. Pion hypothesis PID score: <project:#ParameterSet.pion_pid_cut>.
    3. Separation between two tracks: <project:#ParameterSet.separation_cut>.
    """

    var_names = ["trk_start_x", "trk_start_y", "trk_start_z", "trk_llrpid"]

    # zip is a bit of array restructuring needed for combinations
    # XXX: this only constructs combinations in one way: AxB and not BxA
    combos = ak.combinations(ak.zip({k: arr[k] for k in var_names}), 2)
    protons, pions = ak.unzip(combos)

    # apply PID cuts to the un-zipped per-particle arrays
    # use array.mask to keep the dimensions the same for later calculations
    protons = protons.mask[ak.any(protons.trk_llrpid >= pset.proton_pid_cut, axis=1)]
    pions = pions.mask[ak.any(pions.trk_llrpid >= pset.pion_pid_cut, axis=1)]

    # zip together the list of protons and pions and apply the separation cut
    mask = utils.separation(protons, pions, "trk_start_") >= pset.separation_cut

    return ~ak.is_none(
        # XXX: this reduces any multiple combos to one combo
        ak.firsts(
            ak.fill_none(
                ak.drop_none(
                    ak.zip({"p": protons, "pi": pions})[mask],
                    axis=1,
                ),
                [],
                axis=0,
            )
        )
    )  # type: ignore


def invariant_mass_cut(arr: ak.Array, pset: ParameterSet) -> ak.Array:
    """
    Select events that have a reconstructed {math}`p-\\pi^-` invariant mass
    within the provided parameter range.

    The lower and upper bounds for this cut are defined by
    <project:#ParameterSet.w_lambda_min> and
    <project:#ParameterSet.w_lambda_max> in GeV, respectively.
    """
    raise NotImplementedError
