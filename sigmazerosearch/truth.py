"""
Types for handling truth-level information.

This includes event generator information and that produced by GEANT4 that are
propagated to the ntuples.
"""

from enum import Enum


class GenType(Enum):
    """Represents the neutrino event generator used for an MC sample"""

    GENIE = 0  #: Produced with the [GENIE neutrino event generator](https://github.com/GENIE-MC/Generator)
    NuWro = 1  #: Produced with the [Wroclaw neutrino event generator](https://github.com/nuwro/nuwro)


class GenEventType(Enum):
    """Enum representing the Generator codes for MC events"""

    QEL = 0  #: Quasi-elastic interaction
    RES = 1  #: Resonance
    DIS = 2  #: Deep Inelastic Scattering
    COH = 3  #: Coherent
    ElectronScatter = 5
    MEC = 10  #: Meson-Exchange Current
    Diffractive = 11
    HYP = 1095  #: Hyperon (QEL-event with a hyperon in the final-state)
    Other = -1
