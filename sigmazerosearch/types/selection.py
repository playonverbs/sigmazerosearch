from enum import Enum
import matplotlib.pyplot as plt

# import pandas as pd
import polars as pl
import sigmazerosearch.utils as utils
from sigmazerosearch.types.truth import GenType
from sigmazerosearch.loader.loader import get_POT

"""
ValueUnc represents a central value with either a symmetric or (upper,
lower) associated error
"""
ValueUnc = tuple[float, float] | tuple[float, float, float]


class Cut:
    """Cut represents a single selection cut and the selection state for it."""

    n_signal: ValueUnc
    n_passing: ValueUnc
    n_background: ValueUnc

    def __init__(self, name: str, cutexpr: str):
        self.name: str = name
        self.cutexpr: str = cutexpr
        self.n_passing: ValueUnc = (0.0, 0.0, 0.0)
        self.n_signal: ValueUnc = (0.0, 0.0, 0.0)
        self.n_background: ValueUnc = (0.0, 0.0, 0.0)
        self.applied: bool = False

    def eff(self) -> float:
        """Calculate the selection efficiency at the current Cut"""
        return -1.0

    def pur(self) -> float:
        """Calculate the selection purity at the current Cut"""
        return -1.0

    def _validate_(self) -> bool:
        """validate the cut data"""
        if self.n_signal[0] < 0 or self.n_passing[0] < 0 or self.n_background[0] < 0:
            return False

        return True


# TODO: assign to real type representing a set of parameters
ParameterSet = list[float]  # STUB


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
    # df: pd.DataFrame | None = None
    df: pl.DataFrame | pl.LazyFrame | None = None

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

    async def load_df(self):
        """read file_name into a pandas dataframe"""
        # self.df = pd.read_csv(self.file_name)
        self.df = pl.scan_csv(self.file_name)

    def _validate_(self) -> bool:
        if self.POT < 0:
            return False

        if not utils.file_ok(self.file_name):
            return False

        return True


SampleSet = list[Sample]


class Selection:
    parameters: ParameterSet
    samples: SampleSet
    cuts: list[Cut]

    def __init__(self, **kwargs):
        self.parameters = ParameterSet()
        self.samples = SampleSet()
        self.cuts = []

    @classmethod
    def empty(cls):
        return cls

    @classmethod
    def with_cuts(cls, cuts: list[Cut]):
        out = cls()
        out.cuts = cuts
        return out

    def apply_cut(self, cutname: str):
        for cut in self.cuts:
            if cut.name != cutname:
                continue
            for s in self.samples:
                if isinstance(s.df, pl.DataFrame):
                    s.df[cut.cutexpr]
                else:
                    raise TypeError(f"sample {s.file_name} has not been loaded")

    def plot_eff_pur(self) -> None:
        """plot progressive change in selection purity and efficiency"""
        names: list[str] = [c.name for c in self.cuts]
        effs: list[float] = [c.eff() for c in self.cuts]
        purs: list[float] = [c.pur() for c in self.cuts]

        plt.scatter(names, effs, label="efficiency")
        plt.scatter(names, purs, label="purity")
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

    # def to_latex(self) -> str:
    #     """output selection details to a latex table"""
    #     out = f""
    #     for c in self.cuts:
    #         out += f"{c.name} & {c.n_passing[0]} &"
    #     return out

    # def cut_summary(self):
    #     # for cut in self.cuts:
    #     #     for item in zip(cut.n_signal, cut.n_passing, cut.n_background):
    #     #         item[0]
