from enum import Enum


class GenType(Enum):
    """Represents the neutrino event generator used for an MC sample"""

    GENIE = 0
    NuWro = 1


class GenEventType(Enum):
    """Enum representing the Generator codes for MC events"""

    QEL = 0
    RES = 1
    DIS = 2
    COH = 3
    ElectronScatter = 5
    MEC = 10
    Diffractive = 11
    HYP = 1095
    Other = -1
