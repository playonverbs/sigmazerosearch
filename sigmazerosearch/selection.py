"""
Selection contains the main objects for handling the physics selection.
"""

import inspect
import logging
import pathlib
from collections.abc import Iterable
from dataclasses import asdict, dataclass, field
from enum import Enum, IntEnum, auto
from os.path import isabs
from typing import Callable, Literal, Optional

import awkward as ak
import matplotlib.pyplot as plt
import numpy as np
import uproot as up
from tabulate import tabulate
from uproot.behaviors.TBranch import HasBranches

import sigmazerosearch.alg.fv as fv
import sigmazerosearch.utils as utils
from sigmazerosearch import loader
from sigmazerosearch.general import PDG, Config
from sigmazerosearch.truth import GenEventType, GenType

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

    Added enum members should be added in increasing specificity, due to the
    characteristics of <inv:#numpy.select>.
    """

    @staticmethod
    def _generate_next_value_(name, start, count, last_values):
        if isinstance(last_values[-1], tuple):
            return last_values[-1][0] + 1
        else:
            return last_values[-1] + 1

    Other = -1, None
    Signal = auto(), lambda arr: signal_def(arr)
    Lambda = auto(), lambda arr: arr["mc_hyperon_pdg"] == PDG.Lambda.value
    OOFV = (auto(), lambda arr: arr["sample"] == SampleType.Dirt.name)
    NuMuCC = (
        auto(),
        lambda arr: np.logical_and.reduce(
            [np.abs(arr["mc_nu_pdg"]) == PDG.NuMu.value, arr["mc_ccnc"] == "CC"]
        ),
    )
    NC = (auto(), lambda arr: arr["mc_ccnc"] == "NC")
    NuE = (auto(), lambda arr: np.abs(arr["mc_nu_pdg"]) == PDG.NuE.value)
    NuMuDIS = (
        auto(),
        lambda arr: np.logical_and.reduce(
            [np.abs(arr["mc_nu_pdg"]) == PDG.NuMu.value, arr["mc_mode"] == "DIS"]
        ),
    )
    NuMuRES = (
        auto(),
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
            np.logical_or(
                arr["mc_mode"] == GenEventType.QEL.name,
                arr["mc_mode"] == GenEventType.HYP.name,
            ),
            arr["mc_nu_pdg"] == PDG.NuMu.anti,
            arr["mc_hyperon_pdg"] == PDG.Sigma0.value,
            fv.in_active_tpc(
                arr["mc_nu_pos_x"], arr["mc_nu_pos_y"], arr["mc_nu_pos_z"]
            ),
            ak.sum(arr["mc_decay_pdg"] == PDG.Proton.value, axis=1) >= 1,
            ak.sum(arr["mc_decay_pdg"] == PDG.Pi.anti, axis=1) >= 1,
        )
    )  # type: ignore


def has_hyperon(arr: ak.Array) -> ak.Array:
    """Checks if an event has a S=1 baryon (hyperon) present"""
    return np.logical_or.reduce(
        (
            arr["mc_hyperon_pdg"] == PDG.Lambda.value,
            arr["mc_hyperon_pdg"] == PDG.Sigma0.value,
            arr["mc_hyperon_pdg"] == PDG.SigmaM.value,
            arr["mc_hyperon_pdg"] == PDG.SigmaP.value,
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

    def update(self, arr, cond, scale: float = 1.0, sample=None, signal_def=signal_def):
        """
        Update the cut with the numbers of signal, background and selected
        events.

        These are scaled to the given weighting (eg. POT-based), and signal
        events are only counted if they are from a Hyperon sample.
        Events from dirt samples are scaled to 65% of their original value as
        per NuMI dirt handling procedures.
        """

        if sample:
            match sample.type:
                case SampleType.Hyperon:
                    self.n_signal[0] += scale * ak.sum(signal_def(arr[cond]))
                case SampleType.Dirt:
                    scale *= 0.65

        self.n_background[0] += scale * ak.sum(~signal_def(arr[cond]), axis=None)
        self.n_passing[0] += scale * ak.sum(cond, axis=None)

    def save_state(self) -> dict:
        """
        :::{seealso}
        The documentation for <project:#Selection.save_state>.
        :::
        """
        try:
            cutfunc_source = inspect.getsource(self.cutfunc)
        except Exception:
            cutfunc_source = None

        return {
            "name": self.name,
            "cutfunc": cutfunc_source,
            "n_passing": self.n_passing,
            "n_signal": self.n_signal,
            "n_background": self.n_background,
            "total_signal": self.total_signal,
        }

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


@dataclass(frozen=True, slots=True)
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
    ct_dead_wire_removal: bool = True  # FIXME: remove from defaults

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

    def save_state(self) -> dict:
        """
        :::{seealso}
        The documentation for <project:#Selection.save_state>.
        :::
        """
        return {
            "POT": self.POT,
            "is_data": self.is_data,
            "gen_type": self.gen_type.name,
        }

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
    """
    Represents a hyperon selection: the selection cuts and event samples the
    cuts are applied to. This class also keeps config settings and selection
    cut parameters.
    """

    def __init__(self, **kwargs):
        """
        Creates a selection object through passed kwargs. The keyword arguments
        are given here.

        :param ParameterSet params:
            Parameters that the selection cuts depend upon. Must be defined
            before the Selection is bound.
        :param SampleSet samples:
            The Samples that make up the selection and also the target POT to
            scale events to.
        :param list[Cut] cuts:
            A list of Cuts that will be applied consecutively to the samples.
        :param str label:
            A basic tag for this selection.
        :param Config config:
            Configuration options that change how the selection operates:
            sample iteration step size, plotting directory, etc...
        """

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
                                ak.where(has_hyperon(filter_arr), 0.0, scale)
                                if s.type != SampleType.Hyperon
                                else scale
                            )
                            match s.type:
                                case SampleType.Dirt:
                                    filter_arr["weight"] = (
                                        filter_arr["weight"] * 0.65
                                    )  # NuMI dirt weighting

                            accum = ak.concatenate((accum, filter_arr), axis=0)
            else:
                raise TypeError(f"sample {s.file_name} has not been loaded")

        logger.info(f"Applied cuts: {', '.join([cut.name for cut in self.cuts])}")
        if accumulate:
            return accum

    def optimise_cut(
        self,
        cutfunc: Callable[[ak.Array, dict], ak.Array],
        values: dict[str, Iterable],  # | Iterable[ParameterSet],
        consts: dict = {},
        signal_def: Callable[[ak.Array], ak.Array] = signal_def,
        # eff_func: Callable[[ak.Array], ak.Array] | None = None,
        # pur_func: Callable[[ak.Array], ak.Array] | None = None,
        prev_cuts: list[Cut] | None = None,
    ):
        """
        Repeatedly applies a cut with a varying parameter. Returns the
        efficiency and purity at each variation.
        """

        params_set: list[dict] = []
        for key, vals in values.items():
            for val in vals:
                params_set.append({key: val})

        proxy_cuts = [Cut("proxy_cut", cutfunc) for _ in range(len(params_set))]

        for i, iteration in enumerate(params_set):
            logger.info(f"Iteration {i + 1}")
            logger.info(iteration | consts)
            for s in self.samples:
                scale = self.samples.target_POT / s.POT
                if isinstance(s.df, HasBranches):
                    for arr in loader._yield_array_from_ttree(s.df, self.config):
                        if s.type == SampleType.Hyperon:
                            proxy_cuts[i].total_signal += scale * ak.sum(
                                signal_def(arr), axis=None
                            )
                        cond = (
                            np.logical_and.reduce(
                                [c(arr) for c in prev_cuts]
                                + [proxy_cuts[i](arr, iteration | consts)]
                            )
                            if prev_cuts
                            else proxy_cuts[i](arr, iteration | consts)
                        )
                        proxy_cuts[i].update(
                            arr, cond, scale=scale, sample=s, signal_def=signal_def
                        )
                else:
                    raise TypeError(f"sample {s.file_name} has not been loaded")

        return np.array([[cut.eff(), cut.pur()] for cut in proxy_cuts])

    def save_state(self, arr: ak.Array, output_path: str | pathlib.Path):
        """
        Save the current state to an intermediate ROOT file as a set of
        RDataFrames.

        Includes all defined <project:#Cut>, <project:#Sample> objects and by
        default the filtered entries at the end point of the selection.

        :::{seealso}
        The ROOT documentation for `RDataFrame` can be seen at
        https://root.cern/doc/master/classROOT_1_1RDataFrame.html.
        :::
        """
        samples = ak.Array([sample.save_state() for sample in self.samples])
        cuts = ak.Array([cut.save_state() for cut in self.cuts])
        pset = ak.Array([asdict(self.parameters)])
        # config = ak.Array([asdict(self.config)])
        # branch_list = self.config.branch_list

        with up.create(output_path) as fd:
            fd.mkrntuple("sigmazerosearch/cuts", cuts)
            # fd.mkrntuple("sigmazerosearch/config", config)
            fd.mkrntuple("sigmazerosearch/parameters", pset)
            fd.mkrntuple("sigmazerosearch/samples", samples)
            fd.mkrntuple("sigmazerosearch/events", arr)

        logger.info(f"successfully wrote state to {output_path}")

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

    @staticmethod
    def load_state(filename: str | pathlib.Path):
        """Creates a Selection from a path to a state file"""
        # FIXME: import other exported data: samples, cuts, parameters

        with up.open(filename) as fd:
            root_folder = "sigmazerosearch"

            parameters = ParameterSet(
                **fd.get(f"{root_folder}/parameters").arrays().to_list().pop()  # type: ignore
            )
            events = fd.get(f"{root_folder}/events").arrays()  # type: ignore

        return events, parameters
