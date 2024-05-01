from enum import Enum
import matplotlib.pyplot as plt

import awkward as ak
import numpy as np
from uproot.behaviors.TBranch import HasBranches
import sigmazerosearch.utils as utils
import sigmazerosearch.alg.fv as fv
from sigmazerosearch.types.truth import GenType
from sigmazerosearch.types.general import PDG
from sigmazerosearch.loader.loader import get_POT, load_ntuple

from os.path import isabs
from typing import Callable

"""
ValueUnc represents a central value with either a symmetric or (upper,
lower) associated error.
"""
ValueUnc = tuple[float, float] | tuple[float, float, float]

BRANCH_LIST = None
# BRANCH_LIST = [
#     "reco_primary_vtx_inFV" "mc_hyperon_pdg" "pfp_true_pdg" "NPrimaryTrackDaughters",
#     "NPrimaryShowerDaughters",
#     "RecoPrimaryVertex",
# ]


def signal_def(arr: ak.Array) -> ak.Array:
    """takes an awkward.Array with fields corresponding to ntuple branches,
    applies a mask and returns a boolean array"""
    return np.logical_and.reduce(
        (
            arr["mc_nu_pdg"] == PDG.NuMu.anti,
            arr["mc_hyperon_pdg"] == PDG.Sigma0.value,
            fv.in_active_tpc(
                arr["mc_nu_pos_x"], arr["mc_nu_pos_y"], arr["mc_nu_pos_z"]
            ),
        )
    )  # type: ignore


class Cut:
    """Cut represents a single selection cut and the selection state for it."""

    total_signal = 0.0
    n_signal: ValueUnc
    n_passing: ValueUnc
    n_background: ValueUnc

    def __init__(self, name: str, cutfunc: Callable):
        self.name: str = name
        self.cutfunc: Callable[[ak.Array], ak.Array] = cutfunc
        self.n_passing: ValueUnc = (0.0, 0.0, 0.0)
        self.n_signal: ValueUnc = (0.0, 0.0, 0.0)
        self.n_background: ValueUnc = (0.0, 0.0, 0.0)
        self.applied: bool = False

    def eff(self) -> float:
        """Calculate the selection efficiency at the current Cut"""
        return self.n_signal[0] / self.total_signal

    def pur(self) -> float:
        """Calculate the selection purity at the current Cut"""
        return self.n_signal[0] / self.n_passing[0]

    def _validate_(self) -> bool:
        """validate the cut data"""
        if self.n_signal[0] < 0 or self.n_passing[0] < 0 or self.n_background[0] < 0:
            return False

        return True


# TODO: assign to real type representing a set of parameters, gets passed into
# functions that require a tweakable parameters
ParameterSet = dict[str, bool | float | int]  # STUB


class SampleType(Enum):
    """represents different types of samples being fed into the selection"""

    Data = 0
    Background = 1
    Hyperon = 2
    Dirt = 3
    Ext = 4


class Sample:
    """represents samples and their associated data types"""

    name: str
    file_name: str
    type: SampleType | None = None
    POT: float
    gen_type: GenType
    is_data: bool
    df: HasBranches | None = None

    def __init__(
        self,
        name: str,
        file_name: str,
        type: SampleType,
        POT: float | None,
        is_data: bool = False,
    ):
        self.name = name
        self.file_name = file_name
        self.type = type
        self.POT = POT if POT else get_POT(file_name)  # type: ignore
        self.is_data = is_data

    @classmethod
    def from_dict(cls, kv: dict):
        return cls(kv["name"], kv["file_name"], kv["type"], kv["POT"])

    def load_df(self):
        """read file_name into a pandas dataframe"""
        if not isabs(self.file_name):
            raise OSError
        self.df = load_ntuple(self.file_name + ":ana/OutputTree")

    def _validate_(self) -> bool:
        if self.POT < 0:
            return False

        if not utils.file_ok(self.file_name):
            return False

        return True


class SampleSet(list[Sample]):
    def __init__(self, *samples, **kwargs):
        super().__init__(samples)
        # self.base_dir: str = kwargs["base_dir"] if kwargs["base_dir"] else "."


class Selection:
    # parameters: ParameterSet
    # samples: SampleSet
    # cuts: list[Cut]

    def __init__(self, **kwargs):
        self.parameters: ParameterSet = ParameterSet(kwargs["params"])
        self.samples: SampleSet = kwargs["samples"]
        self.cuts: list[Cut] = kwargs["cuts"]

    @classmethod
    def empty(cls):
        return cls

    @classmethod
    def with_cuts(cls, cuts: list[Cut]):
        out = cls()
        out.cuts = cuts
        return out

    def apply_cut(self, cutname: str):
        """
        Apply a given selection cut's cut function to the sample arrays and
        accumulates the resulting number of signal, background and total
        passing particles per cut.
        """
        for i, cut in enumerate(self.cuts):
            if cut.name != cutname:
                continue
            for s in self.samples:
                if isinstance(s.df, HasBranches):
                    arr = s.df.arrays(BRANCH_LIST)
                    if Cut.total_signal == 0.0:
                        Cut.total_signal = ak.sum(signal_def(arr), axis=None)  # type: ignore
                    cond = np.logical_and.reduce(
                        [c.cutfunc(arr) for c in self.cuts[: i + 1]]
                    )
                    arr = arr[cond]
                    cut.n_signal = (ak.sum(signal_def(arr), axis=None), 0, 0)  # type: ignore
                    cut.n_background = (ak.sum(~signal_def(arr), axis=None), 0, 0)  # type: ignore
                    cut.n_passing = (ak.sum(cond, axis=None), 0, 0)  # type: ignore
                    cut.applied = True
                else:
                    raise TypeError(f"sample {s.file_name} has not been loaded")

    def plot_reco_effs(self, signal=True) -> None:
        pdgs = [PDG.Photon.value, PDG.Proton.value, PDG.Pi.anti, PDG.Muon.anti]
        lost = []
        counted = []
        # TODO: fix for multiple samples

        for s in self.samples:
            if isinstance(s.df, HasBranches):
                arr = s.df.arrays(BRANCH_LIST)
                for pdg in pdgs:
                    cond = ak.sum(arr["pfp_true_pdg"] == pdg, axis=1) >= 1  # type: ignore
                    if signal:
                        counted.append(len(arr[signal_def(arr[cond])]))
                        lost.append(len(arr[signal_def(arr[~cond])]))
                    else:
                        counted.append(len(arr[cond]))
                        lost.append(len(arr[~cond]))

        labels = [r"$\gamma$", r"$p$", r"$\pi^-$", r"$\mu^+$"]

        fig, ax = plt.subplots()
        ax.set_title("Reco. Efficiency", loc="right", color="grey", weight="bold")
        ax.set_ylabel("# Events with Particles")
        ax.bar(labels, counted, label="found", bottom=[0, 0, 0, 0])
        ax.bar(labels, lost, label="lost", bottom=counted)
        ax.legend()
        fig.tight_layout()
        plt.show()

    def plot_eff_pur(self, exp: bool = False) -> None:
        """plot progressive change in selection purity and efficiency"""
        names: list[str] = [c.name for c in self.cuts]
        effs: list[float] = [c.eff() for c in self.cuts]
        purs: list[float] = [c.pur() for c in self.cuts]

        fig, ax = plt.subplots()
        ax.set_title("Selection Performance", loc="right", color="grey", weight="bold")
        if exp:
            ax.set_ylabel(r"$\epsilon p$")
            ax.plot(
                names,
                [e * p for e, p in zip(effs, purs)],
                label="eff * pur",
                marker="o",
            )
            fig.tight_layout()
            plt.show()
        else:
            (e,) = ax.plot(
                names, effs, label="efficiency", color="tab:blue", marker="o"
            )
            ax2 = ax.twinx()
            (p,) = ax2.plot(names, purs, label="purity", color="tab:orange", marker="o")
            ax.set_ylabel("Efficiency")
            ax2.set_ylabel("Purity")
            ax.legend([e, p], ["Efficiency", "Purity"])
            fig.tight_layout()
            plt.show()

    def plot_slice_info(self, type="both", signal=True) -> None:
        if type not in ["both", "purity", "completeness"]:
            raise TypeError(
                '''type must be one of "purity", "completeness" or "both"'''
            )

        fig, ax = plt.subplots(tight_layout=True)

        arr = self.samples[0].df.arrays(BRANCH_LIST)  # type: ignore

        cond = signal_def(arr) if signal else True

        if type == "both":
            ax.set_xlabel(r"True $\nu$ slice completeness")
            ax.set_ylabel(r"True $\nu$ slice purity")
            h2 = ax.hist2d(
                arr["true_nu_slice_completeness"][
                    (arr["true_nu_slice_completeness"] != -999) & cond
                ].to_numpy(),
                arr["true_nu_slice_purity"][
                    (arr["true_nu_slice_purity"] != -999) & cond
                ].to_numpy(),
            )
            fig.colorbar(h2[3], ax=ax, label="Number of Slices")
        elif type == "purity":
            ax.set_xlabel(r"True $\nu$ slice purity")
            ax.hist(
                arr["true_nu_slice_completeness"][
                    (arr["true_nu_slice_completeness"] != -999) & cond
                ].to_numpy(),
                histtype="step",
            )
        elif type == "completeness":
            ax.set_xlabel(r"True $\nu$ slice completeness")
            ax.hist(
                arr["true_nu_slice_completeness"][
                    (arr["true_nu_slice_completeness"] != -999) & cond
                ].to_numpy(),
                histtype="step",
            )

        plt.show()

    def sample_types(self) -> list[SampleType | None]:
        """list types of all assoc samples"""
        return [sample.type for sample in self.samples]

    def _validate_cuts_(self) -> bool:
        """run validation for all assoc cuts"""
        for cut in self.cuts:
            if not cut._validate_():
                return False

        return True

    def open_files(self) -> None:
        """load all samples into dataframes asynchronously (for now)"""
        for sample in self.samples:
            sample.load_df()

    def close_files(self) -> None:
        """
        delete dataframe objects from memory, to be run after IO operations
        have been run on samples
        """
        for s in self.samples:
            del s.df  # NOTE: maybe naive; refactor when final DataFrame chosen

    # def to_latex(self) -> str:
    #     """output selection details to a latex table"""
    #     out = f""
    #     for c in self.cuts:
    #         out += f"{c.name} & {c.n_passing[0]} &"
    #     return out

    def cut_summary(self, header: bool = False, format: str = "text"):
        if format not in ["csv", "text", "latex"]:
            raise TypeError('''format must be one of "csv", "text", "latex"''')

        if format == "text":
            header_row = "{:<20} {:>10} {:>10} {:>10} {:>10}"
            row = "{:<20} {:>10} {:>10} {:>10.5f} {:>10.5f}"
            if header:
                print(
                    header_row.format(
                        "Cut name", "Signal", "Background", "Eff.", "Pur."
                    )
                )
            for cut in self.cuts:
                print(
                    row.format(
                        cut.name,
                        cut.n_signal[0],
                        cut.n_background[0],
                        cut.eff(),
                        cut.pur(),
                    )
                )
        elif format == "csv":
            header_row = "{:<},{:>},{:>},{:>},{:>}"
            row = "{:<},{:>},{:>},{:>.5f},{:>.5f}"
            if header:
                print(
                    header_row.format(
                        "Cut name", "Signal", "Background", "Eff.", "Pur."
                    )
                )
            for cut in self.cuts:
                print(
                    row.format(
                        cut.name,
                        cut.n_signal[0],
                        cut.n_background[0],
                        cut.eff(),
                        cut.pur(),
                    )
                )
        elif format == "latex":
            header_row = "{:<20} & {:>10} & {:>10} & {:>10} & {:>10} \\"
            row = "{:<20} &  ${:>10}$ &  ${:>10}$ & ${:>10.5f}$ & ${:>10.5f}$ \\"
            if header:
                print(
                    header_row.format(
                        "Cut name", "Signal", "Background", "Eff.", "Pur."
                    )
                )
            for cut in self.cuts:
                print(
                    row.format(
                        cut.name,
                        cut.n_signal[0],
                        cut.n_background[0],
                        cut.eff(),
                        cut.pur(),
                    )
                )
