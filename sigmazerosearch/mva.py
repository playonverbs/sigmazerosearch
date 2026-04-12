"""
Provdies handlers that wrap common mutli-variate analysis preparation steps.

This class does not operate algorithms such as gradient-descent BDTs in this
package, but prepares data to be passed to the ROOT TMVA package.
"""

import dataclasses
import logging
from pathlib import Path
from typing import Callable

import awkward as ak
import uproot as up

logger = logging.getLogger(__name__)


@dataclasses.dataclass(kw_only=True, eq=False, frozen=True)
class MVAHandler:
    """
    Class encapsulating data needed for reading, filtering and writing MVA
    TTrees.

    Given an array of per-event information output by
    <project:#Selection.apply_cut> and signal/background criteria
    """

    data: ak.Array
    """The array produced by <project:#Selection.apply_cut>, acting as the input."""
    target_branches: set[str]
    """The subset of input branches to be written to the BDT trees."""
    signal_condition: Callable
    """A function returning a boolean <inv:#ak.Array> indicating which entries to keep."""
    background_condition: Callable
    precondition: Callable | None = None
    """A filter applied to all events before the signal/background split"""
    extra_fields: dict[str, Callable[[ak.Array], ak.Array]] = dataclasses.field(
        default_factory=dict
    )
    """Additional fields to be computed using and added to `data`."""
    mva_method: str = "bdt"
    """MVA method string, used as a TDirectory name inside the output file"""
    flatten: bool = False
    """
    Flag to flatten arrays if picking variables at an axis deeper than 0. e.g.
    per-PFP variables
    """

    def __post_init__(self):
        for k, func in self.extra_fields.items():
            self.data[k] = func(self.data)

    def __restructure_array(self) -> ak.Array:
        d = dict(
            zip(
                ak.fields(self.data[self.target_branches]),
                ak.unzip(
                    self.data[self.target_branches][self.precondition(self.data)]
                    if self.precondition is not None
                    else self.data[self.target_branches]
                ),
            )
        )

        if self.flatten:
            for k, v in d.items():
                d[k] = ak.flatten(v)

        return ak.zip(d)

    def save_bdt_trees(self, output_path: Path) -> tuple[int, int]:
        """
        Restructure the input data and output a TFile with a TTree for signal
        entries and one for background entries.

        :param Path output_path: The file to output the tree to.

        :return: A tuple counting `(signal_entries, background_entries)`.
        """
        signal_key = f"{self.mva_method}/SignalTree"
        background_key = f"{self.mva_method}/BackgroundTree"

        output = self.__restructure_array()
        bdt_dir = up.recreate(output_path)

        logger.info(f"writing trees to {output_path}")

        # apply signal and background to preconditioned original data array to
        # allow the use of all included fields, and keeping dimensions correct.
        bdt_dir[signal_key] = output[
            self.signal_condition(self.data[self.precondition(self.data)])
            if self.precondition is not None
            else self.signal_condition(self.data)
        ]
        bdt_dir[background_key] = output[
            self.background_condition(self.data[self.precondition(self.data)])
            if self.precondition is not None
            else self.background_condition(self.data)
        ]

        logger.info(
            f"Saved {bdt_dir[signal_key].num_entries} Signal entries and {bdt_dir[background_key].num_entries} Background entries"  # type: ignore
        )

        return bdt_dir[signal_key].num_entries, bdt_dir[background_key].num_entries


def precond_has_pfps(array: ak.Array):
    """Check if PFPs exist in the event via the `pfp_true_pdg` field"""
    return ak.num(array.pfp_true_pdg) > 0
