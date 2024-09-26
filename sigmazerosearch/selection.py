"""
Selection contains the main objects for handling the physics selection.
"""

from enum import Enum, IntEnum
from os.path import isabs
from typing import Callable

import awkward as ak
import matplotlib.pyplot as plt
import numpy as np
from tabulate import tabulate
from uproot.behaviors.TBranch import HasBranches

import sigmazerosearch.alg.fv as fv
import sigmazerosearch.utils as utils
from sigmazerosearch.general import PDG, Config
from sigmazerosearch.loader import _yield_array_from_ttree, get_POT, load_ntuple
from sigmazerosearch.truth import GenType

# ValueUnc = tuple[float, float] | tuple[float, float, float]
ValueUnc = list[float]
"""
ValueUnc represents a central value with either a symmetric or (upper,
lower) associated error.
"""


class EventCategory(IntEnum):
    Signal = 0
    Lambda = 1
    NuMuCC = 2
    NC = 3
    Other = 4

    @staticmethod
    def _signal(arr):
        return signal_def(arr) * EventCategory.Signal.value

    @staticmethod
    def _lambda(arr):
        return (arr["mc_hyperon_pdg"] == PDG.Lambda.value) * EventCategory.Lambda.value

    @staticmethod
    def _numucc(arr):
        return (
            np.logical_and.reduce(
                [
                    np.abs(arr["mc_nu_pdg"]) == PDG.NuMu.value,
                    np.abs(arr["mc_lepton_pdg"]) == PDG.Muon.value,
                ]
            )
            * EventCategory.NuMuCC.value
        )

    @staticmethod
    def _nc(arr):
        return (
            np.logical_and.reduce(
                [
                    np.abs(arr["mc_nu_pdg"]) == PDG.NuMu.value,
                    np.abs(arr["mc_lepton_pdg"]) == PDG.NuMu.value,
                ]
            )
            * EventCategory.NC.value
        )

    @classmethod
    def from_arr(cls, arr):
        return cls._signal(arr) | cls._lambda(arr) | cls._numucc(arr) | cls._nc(arr)


def signal_def(arr: ak.Array) -> ak.Array:
    """Takes an <inv:#ak.Array> with fields corresponding to ntuple branches,
    applies a mask and returns a boolean array"""
    return np.logical_and.reduce(
        (
            arr["mc_nu_pdg"] == PDG.NuMu.anti,
            arr["mc_hyperon_pdg"] == PDG.Sigma0.value,
            fv.in_active_tpc(
                arr["mc_nu_pos_x"], arr["mc_nu_pos_y"], arr["mc_nu_pos_z"]
            ),
            ak.sum(arr["mc_decay_pdg"] == PDG.Proton.value, axis=1) >= 1,
            ak.sum(arr["mc_decay_pdg"] == PDG.Pi.anti, axis=1) >= 1,
        )
    )  # type: ignore


class Cut:
    """Cut represents a single selection cut and the selection state for it."""

    total_signal = 0.0

    def __init__(self, name: str, cutfunc: Callable):
        self.name: str = name
        self.cutfunc: Callable[[ak.Array], ak.Array] = cutfunc
        self.n_passing: ValueUnc = [0.0, 0.0, 0.0]
        self.n_signal: ValueUnc = [0.0, 0.0, 0.0]
        self.n_background: ValueUnc = [0.0, 0.0, 0.0]
        self.applied: bool = False

    def eff(self) -> float:
        """Calculate the selection efficiency at the current Cut"""
        return self.n_signal[0] / self.total_signal

    def pur(self) -> float:
        """Calculate the selection purity at the current Cut"""
        return self.n_signal[0] / self.n_passing[0]

    def update(self, arr, cond):
        self.n_signal[0] += ak.sum(signal_def(arr[cond]))
        self.n_background[0] += ak.sum(~signal_def(arr[cond]), axis=None)
        self.n_passing[0] += ak.sum(cond, axis=None)

    def __call__(self, *args):
        """Allow an instance of Cut to be used like its cutfunc"""
        return self.cutfunc(*args)

    def _validate_(self) -> bool:
        """validate the cut data"""
        if self.n_signal[0] < 0 or self.n_passing[0] < 0 or self.n_background[0] < 0:
            return False

        return True


# TODO: assign to real type representing a set of parameters, gets passed into
# functions that require a tweakable parameters
ParameterSet = dict[str, bool | float | int]  # STUB


class SampleType(Enum):
    """Represents different types of samples being fed into the selection"""

    Data = 0
    Background = 1
    Hyperon = 2
    Dirt = 3
    Ext = 4


class Sample:
    """Represents samples and their associated data types"""

    def __init__(
        self,
        name: str,
        file_name: str,
        type: SampleType,
        POT: float | None,
        is_data: bool = False,
        gen_type: GenType = GenType.GENIE,
    ):
        self.name: str = name
        self.file_name: str = file_name
        self.type: SampleType = type
        self.gen_type: GenType = gen_type
        self.POT: float = POT if POT else get_POT(file_name)  # type: ignore
        self.is_data: bool = is_data
        self.df: HasBranches | None = None

    @classmethod
    def from_dict(cls, kv: dict):
        return cls(kv["name"], kv["file_name"], kv["type"], kv["POT"])

    def load_df(self):
        """Read file_name into an awkward.Array"""
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
    """An ordered collection of <project:#Sample> objects."""

    def __init__(self, *samples, **kwargs):
        super().__init__(samples)
        self.target_POT: float | None = kwargs.get("target_POT")
        # self.base_dir: str = kwargs["base_dir"] if kwargs["base_dir"] else "."


class Selection:
    def __init__(self, **kwargs):
        self.parameters: ParameterSet = ParameterSet(kwargs["params"])
        self.samples: SampleSet = kwargs["samples"]
        self.cuts: list[Cut] = kwargs["cuts"]
        self.label: str = "_" + kwargs["label"] if kwargs.get("label") else ""
        self.config: Config = kwargs.get("config", Config.default())
        self.config.validate()

    def apply_cut(self, cutname: str):
        """
        Apply a given selection cut's cut function to the sample arrays and
        accumulates the resulting number of signal, background and total
        passing particles per cut.

        Returns the destructively masked array after cut conditions have been
        applied.
        """
        for i, cut in enumerate(self.cuts):
            if cut.name != cutname:
                continue
            for s in self.samples:
                if isinstance(s.df, HasBranches):
                    # arr = s.df.arrays(self.config.branch_list)
                    for arr in _yield_array_from_ttree(s.df, self.config):
                        # only count total signal for the initial cut
                        if i == 0:
                            Cut.total_signal += ak.sum(signal_def(arr), axis=None)  # type: ignore
                        cond = np.logical_and.reduce(
                            [c.cutfunc(arr) for c in self.cuts[: i + 1]]
                        )
                        cut.update(arr, cond)
                    # return arr
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
                arr = s.df.arrays(self.config.branch_list)
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
        if self.config.plot_save:
            utils._save_plot(self.config, fig, f"particle_reco_efficiency{self.label}")
        plt.show()

    def plot_eff_pur(self, exp: bool = False) -> None:
        """
        Plot progressive change in selection purity and efficiency as a
        function of `Cut`
        """
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
        if self.config.plot_save:
            utils._save_plot(self.config, fig, f"selection_performance{self.label}")
        plt.show()

    def plot_slice_info(self, type="both", signal=True) -> None:
        if type not in ["both", "purity", "completeness"]:
            raise TypeError(
                '''type must be one of "purity", "completeness" or "both"'''
            )

        fig, ax = plt.subplots(tight_layout=True)

        title = "Slice Info (Signal)" if signal else "Slice Info (All)"
        ax.set_title(title, loc="right", color="grey", weight="bold")

        arr = self.samples[0].df.arrays(self.config.branch_list)  # type: ignore

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

        fig.tight_layout()
        if self.config.plot_save:
            utils._save_plot(self.config, fig, f"slice_info{self.label}")
        plt.show()

    def sample_types(self) -> list[SampleType | None]:
        """List types of all samples associated with this Selection"""
        return [sample.type for sample in self.samples]

    def _validate_cuts_(self) -> bool:
        """run validation for all assoc cuts"""
        for cut in self.cuts:
            if not cut._validate_():
                return False

        return True

    def open_files(self) -> None:
        """load all samples into dataframes synchronously (for now)"""
        for sample in self.samples:
            sample.load_df()

    def close_files(self) -> None:
        """
        delete dataframe objects from memory, to be run after IO operations
        have been run on samples
        """
        for s in self.samples:
            del s.df  # NOTE: maybe naive; refactor when final DataFrame chosen

    def cut_summary(self, header: bool = False, format: str = "text"):
        def print_table(format: str = "simple"):
            headers = (
                ["Cut name", "Signal", "Background", "Eff.", "Pur."] if header else []
            )
            print(
                tabulate(
                    [
                        [
                            cut.name,
                            cut.n_signal[0],
                            cut.n_background[0],
                            cut.eff(),
                            cut.pur(),
                        ]
                        for cut in self.cuts
                    ],
                    headers=headers,
                    tablefmt=format,
                    floatfmt=("", "", "", ".5f", ".5f"),
                )
            )

        if format == "text":
            print_table()
        elif format in ["markdown", "md"]:
            print_table("github")
        elif format == "latex":
            print_table("latex")
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
        else:
            print_table(format)
