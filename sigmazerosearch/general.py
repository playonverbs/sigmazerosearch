"""
General types to aid in configuring the framework.
"""

from dataclasses import dataclass
from enum import IntEnum
from pathlib import Path
from typing import Iterable


@dataclass
class Config:
    """
    Config represents configuration options that will affect the operation of
    the selection framework but should not effect the physics choices made by
    the framework.
    """

    plot_save: bool = False
    plot_dir: Path | None = None
    plot_format: str | Iterable[str] = "png"
    branch_list: Iterable[str] | None = None
    iterate: bool = False
    iterate_step: int | str | None = None

    def __post_init__(self):
        self.validate()

    @classmethod
    def default(cls):
        return cls()

    def validate(self):
        """
        Checks if certain config value combinations have been set together.

        For example: plotting requires an output format(s) and an output
        directory
        """
        if self.plot_save and self.plot_dir is None:
            raise Exception("ensure all required plot options are set")

        if self.plot_dir is not None and not self.plot_dir.is_dir():
            raise ValueError("plot_dir should be a path to a directory")

        if self.iterate and self.iterate_step is None:
            raise ValueError("iterate should be set with iterate_step")

        if not self.iterate and self.iterate_step is not None:
            raise ValueError("iterate should be set with iterate_step")


class PDG(IntEnum):
    """
    Enum representing the PDG codes relevant to the analysis.
    """

    E = 11
    NuE = 12
    Muon = 13
    NuMu = 14
    Photon = 22

    Pi0 = 111
    Pi = 211
    Kaon0 = 311
    Kaon = 321
    Neutron = 2112
    Proton = 2212
    Lambda = 3122
    Sigma0 = 3212

    @property
    def anti(self) -> int:
        """Gets the anti-particle pdg code value excluding neutral particles"""
        neutrals = [PDG.Lambda, PDG.Neutron, PDG.Photon] + [
            x for x in PDG if x.name.endswith("0")
        ]
        if self in neutrals:
            return self.value
        return self.value.__neg__()

    def __invert__(self):
        return self.anti

    def __neg__(self):
        return self.anti
