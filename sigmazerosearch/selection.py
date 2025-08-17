"""
Selection contains the main objects for handling the physics selection.
"""

import logging
from dataclasses import asdict, dataclass, field
from enum import Enum, IntEnum
from os.path import isabs
from typing import Callable, Literal, Optional

import awkward as ak
import matplotlib.pyplot as plt
import numpy as np
from tabulate import tabulate
from uproot.behaviors.TBranch import HasBranches

import sigmazerosearch.alg.fv as fv
import sigmazerosearch.utils as utils
from sigmazerosearch import loader
from sigmazerosearch.general import PDG, Config
from sigmazerosearch.truth import GenType

# ValueUnc = tuple[float, float] | tuple[float, float, float]
ValueUnc = list[float]
"""
ValueUnc represents a central value with either a symmetric or (upper,
lower) associated error.
"""

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class EventCategoryMixin:
    index: int
    func: Callable | None = field(repr=False)


class EventCategory(EventCategoryMixin, Enum):
    """
    Represents categories that MC neutrino interaction events can fit into.
    Members have an enum index and a function for applying to an array of
    events.

    If members are used as values (for comparison, etc...) they should be
    accessed via `.index` or `.name`:
    ```{code-block} python
    EventCategory.from_arr(array) == EventCategory.Signal.name

    EventCategory.from_arr(array, "code") == EventCategory.Signal.index
    ```
    """

    Other = -1, None
    Signal = 0, lambda arr: signal_def(arr)
    Lambda = 1, lambda arr: arr["mc_hyperon_pdg"] == PDG.Lambda.value
    NuMuCC = (
        2,
        lambda arr: np.logical_and.reduce(
            [np.abs(arr["mc_nu_pdg"]) == PDG.NuMu.value, arr["mc_ccnc"] == "CC"]
        ),
    )
    NC = (3, lambda arr: arr["mc_ccnc"] == "NC")
    NuE = (4, lambda arr: np.abs(arr["mc_nu_pdg"]) == PDG.NuE.value)
    NuMuDIS = (
        5,
        lambda arr: np.logical_and.reduce(
            [np.abs(arr["mc_nu_pdg"]) == PDG.NuMu.value, arr["mc_mode"] == "DIS"]
        ),
    )
    NuMuRES = (
        6,
        lambda arr: np.logical_and.reduce(
            [np.abs(arr["mc_nu_pdg"]) == PDG.NuMu.value, arr["mc_mode"] == "RES"]
        ),
    )

    def __call__(self, values):
        if self.func is None:
            raise TypeError(f"{self}.func is None")

        return self.func(values)

    @staticmethod
    def from_arr(arr, output: Literal["string", "code"] = "string"):
        """
        Returns the enum value per event given which category each event
        satisfies.
        """
        if output == "string":
            return np.select(
                [ec(arr) for ec in list(EventCategory) if ec.func is not None],
                [
                    [ec.name] * len(arr)
                    for ec in list(EventCategory)
                    if ec.func is not None
                ],
                default=EventCategory.Other.index,
            )
        else:
            return np.select(
                [ec(arr) for ec in list(EventCategory) if ec.func is not None],
                [
                    [ec.index] * len(arr)
                    for ec in list(EventCategory)
                    if ec.func is not None
                ],
                default=EventCategory.Other.index,
            )


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

    def __init__(self, name: str, cutfunc: Callable):
        self.name: str = name
        self.cutfunc: Callable[[ak.Array], ak.Array] = cutfunc
        self.n_passing: ValueUnc = [0.0, 0.0, 0.0]
        self.n_signal: ValueUnc = [0.0, 0.0, 0.0]
        self.n_background: ValueUnc = [0.0, 0.0, 0.0]
        self.applied: bool = False
        self.total_signal: float = 0.0

    def eff(self) -> float:
        """Calculate the selection efficiency at the current Cut"""
        return self.n_signal[0] / self.total_signal

    def pur(self) -> float:
        """Calculate the selection purity at the current Cut"""
        return self.n_signal[0] / self.n_passing[0]

    def update(self, arr, cond, scale: float = 1.0, sample=None):
        """
        Update the cut with the numbers of signal, background and selected
        events.

        These are scaled to the given weighting (eg. POT-based), and signal
        events are only counted if they are from a Hyperon sample.
        """
        if sample and sample.type == SampleType.Hyperon:
            self.n_signal[0] += scale * ak.sum(signal_def(arr[cond]))
        self.n_background[0] += scale * ak.sum(~signal_def(arr[cond]), axis=None)
        self.n_passing[0] += scale * ak.sum(cond, axis=None)

    def __call__(self, *args):
        """Allow an instance of Cut to be used like its cutfunc"""
        return self.cutfunc(*args)

    def _validate_(self) -> bool:
        """validate the cut data"""
        if self.n_signal[0] < 0 or self.n_passing[0] < 0 or self.n_background[0] < 0:
            return False

        return True

    @classmethod
    def null(cls):
        """Returns the 'null' or identity Cut that passes all events"""
        return Cut("null", lambda arr: len(arr) * [True])

    def __repr__(self) -> str:
        return f"<Cut name={self.name} passing={self.n_passing} signal={self.n_signal} background={self.n_background}>"


@dataclass(frozen=True)
class ParameterSet:
    """
    Wraps all selection parameter values.

    All parameters are required even if only a subset of selection cuts are
    chosen.
    """

    max_separation: float
    """Units: cm"""
    min_length: float
    """Units: cm"""
    pid_cut: float
    proton_pid_cut: float
    pion_pid_cut: float
    separation_cut: float
    """Units: cm"""
    w_lambda_min: float
    """Units: GeV"""
    w_lambda_max: float
    """Units: GeV"""
    ct_time_bins: int = 250  # FIXME: remove from defaults
    ct_wire_window: int = 100  # FIXME: remove from defaults
    ct_island_size: int | None = 6  # FIXME: remove from defaults

    @staticmethod
    def from_dict(kwargs):
        return ParameterSet(**kwargs)

    def as_table(self, file=None):
        """
        Print this ParameterSet in a tabular form
        """
        print(
            tabulate(
                asdict(self),
                tablefmt="rounded_outline",
                showindex="always",
            ),
            file=file,
        )

    def with_(self, kvs: dict):
        return ParameterSet.from_dict(asdict(self) | kvs)


class SampleType(IntEnum):
    """
    Represents different types of samples being fed into the selection, both
    Monte-Carlo (simulated) files and directly recorded files.
    """

    Data = 0
    """Produced from files recorded from the detector."""
    Background = 1
    """Combined Background MC."""
    Hyperon = 2
    """All enriched hyperon interactions."""
    Dirt = 3
    """Out of TPC interactions."""
    EXT = 4
    """Interactions originating from cosmic origins."""


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
        self.POT: float = POT if POT else loader.get_POT(file_name)  # type: ignore
        self.is_data: bool = is_data
        self.df: HasBranches | None = None

    @classmethod
    def from_dict(cls, kv: dict):
        return cls(kv["name"], kv["file_name"], kv["type"], kv["POT"])

    def load_df(self):
        """Read file_name into an awkward.Array"""
        if not isabs(self.file_name):
            raise OSError
        self.df = loader.load_ntuple(self.file_name + ":ana/OutputTree")

    def _validate_(self) -> bool:
        if self.POT < 0:
            return False

        if not utils.file_ok(self.file_name):
            return False

        return True

    def __repr__(self) -> str:
        return f"<Sample name={self.name} file={self.file_name} type={self.type.name} gen={self.gen_type.name} POT={self.POT}>"


class SampleSet(list[Sample]):
    """An ordered collection of <project:#Sample> objects."""

    def __init__(self, *samples, **kwargs):
        super().__init__(samples)
        self.target_POT: float | None = kwargs.get("target_POT")
        # self.base_dir: str = kwargs["base_dir"] if kwargs["base_dir"] else "."

    def as_table(self, file=None):
        """
        Print samples and the target POT in a tabular form.
        """
        print(f"target_POT = {self.target_POT} POT")
        print(
            tabulate(
                [
                    [sam.file_name.split("/")[-1], sam.type.name, sam.POT]
                    for sam in self
                ],
                ["filename", "type", "POT"],
                tablefmt="rounded_outline",
                showindex="always",
            ),
            file=file,
        )


class Selection:
    def __init__(self, **kwargs):
        self.parameters: ParameterSet = kwargs["params"]
        self.samples: SampleSet = kwargs["samples"]
        self.cuts: list[Cut] = kwargs["cuts"]
        self.label: str = "_" + kwargs["label"] if kwargs.get("label") else ""
        self.config: Config = kwargs.get("config", Config.default())
        self.dead_wire_map = loader.load_dead_wire_map(
            self.config.data_dir / "MCC9_channel_list.txt"
        )
        self.config.validate()

    def apply_cut(
        self, cuts: list[Cut], accumulate: bool = False
    ) -> Optional[ak.Array]:
        """
        Apply a given selection cut's cut function to the sample arrays and
        accumulates the resulting number of signal, background and total
        passing particles per cut.

        When `accumulate = True`, this also adds the following fields to the
        output array:

        :sample: The string representation of the sample type the entry is from
        :cut: The name of the latest cut the event has passed.
        :weight: A per-event weighting (currently based only on the sample POT)
        """

        accum = ak.Array([]) if accumulate else None
        for s in self.samples:
            scale = self.samples.target_POT / s.POT
            logger.info(f"Opening {s.file_name.split('/')[-1]}:")
            if isinstance(s.df, HasBranches):
                for i, cut in enumerate(cuts):
                    logger.debug(f"Applying cut {cut.name} to {s.name}")
                    for arr in loader._yield_array_from_ttree(s.df, self.config):
                        if s.type == SampleType.Hyperon:
                            cut.total_signal += scale * ak.sum(
                                signal_def(arr), axis=None
                            )
                        cond = np.logical_and.reduce([c(arr) for c in cuts[: i + 1]])
                        cut.update(arr, cond, scale=scale, sample=s)

                        filter_arr = arr[cond]

                        # build array only for the last cut (avoid double-counting)
                        if accumulate and (i == len(cuts) - 1):
                            filter_arr["sample"] = s.type.name
                            filter_arr["cut"] = cut.name

                            # POT weight: null weight if signal found in non-hyperon file
                            filter_arr["weight"] = (
                                ak.where(signal_def(filter_arr), 0.0, scale)
                                if s.type != SampleType.Hyperon
                                else scale
                            )
                            accum = ak.concatenate((accum, filter_arr), axis=0)
            else:
                raise TypeError(f"sample {s.file_name} has not been loaded")

        logger.info(f"Applied cuts: {', '.join([cut.name for cut in self.cuts])}")
        if accumulate:
            return accum

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

        fig, ax = plt.subplots(figsize=(6, 7), layout="constrained")
        ax.set_title("Selection Performance", loc="right", weight="bold")
        if exp:
            ax.set_ylabel(r"$\epsilon p$")
            ax.plot(
                names,
                [e * p for e, p in zip(effs, purs)],
                label="eff * pur",
                marker="o",
            )
        else:
            # change to percentage if flag set
            effs = [eff * 100 for eff in effs] if self.config.perf_percent else effs
            purs = [pur * 100 for pur in purs] if self.config.perf_percent else purs

            (e,) = ax.plot(
                names, effs, label="efficiency", color="tab:blue", marker="o"
            )
            ax2 = ax.twinx()
            (p,) = ax2.plot(names, purs, label="purity", color="tab:orange", marker="^")
            ax.set_ylabel(
                "Efficiency{percent}".format(
                    percent=r" [\%]" if self.config.perf_percent else ""
                )
            )
            ax2.set_ylabel(
                "Purity{percent}".format(
                    percent=r" [\%]" if self.config.perf_percent else ""
                )
            )
            ax.legend([e, p], ["Efficiency", "Purity"], loc="upper center")

        ax.grid(True, "major", "y")
        ax.set_ylim(0, 100)
        ax.tick_params("x", rotation=90)

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
        """Run validation for all associated cuts"""
        for cut in self.cuts:
            if not cut._validate_():
                return False

        return True

    def open_files(self) -> None:
        """Load all samples into dataframes synchronously (for now)"""
        for sample in self.samples:
            sample.load_df()

    def close_files(self) -> None:
        """
        Delete dataframe objects from memory, to be run after IO operations
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
                            cut.eff() * 100 if self.config.perf_percent else cut.eff(),
                            cut.pur() * 100 if self.config.perf_percent else cut.pur(),
                        ]
                        for cut in self.cuts
                    ],
                    headers=headers,
                    tablefmt=format,
                    floatfmt=("", ".2f", ".2f", ".5f", ".5f"),
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
            row = "{:<},{:>.2f},{:>.2f},{:>.5f},{:>.5f}"
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
