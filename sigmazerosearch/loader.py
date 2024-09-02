"""
Abstracted data-file loading procedures.

Currently ROOT NTuple files are accepted (as per the output of
HyperonProduction). However other data files may be added in the future.
"""

from os.path import isabs

import awkward as ak
import uproot as up
from uproot.behaviors.TBranch import HasBranches


def load_ntuple(filename: str) -> HasBranches:
    """
    Wraps the uproot.open method, taking a filename and outputting some ROOT
    object that has branches
    """
    data = up.open(filename)
    if not isabs(filename):
        raise OSError("Please provide an absolute file path")
    if isinstance(data, HasBranches):
        return data
    else:
        raise TypeError("please read a TTree using the ':' separator.")


def get_POT(filename: str) -> float:
    """Sums the POT value of each subrun from a given ROOT file"""
    if not isabs(filename):
        raise OSError("Please provide an absolute file path")
    data = up.open(filename)
    pots_per_subrun: ak.Array = data.get("ana/MetaTree").get("POT").array()  # type: ignore
    return ak.sum(pots_per_subrun)
