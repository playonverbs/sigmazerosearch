from enum import Enum


class PDG(Enum):
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
        """gets the anti-particle pdg code value excluding neutral particles"""
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
