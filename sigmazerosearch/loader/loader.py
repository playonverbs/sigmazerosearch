import uproot as up
import awkward as aw
from uproot.behaviors.TBranch import HasBranches


def load_ntuple(filename: str) -> HasBranches:
    """wraps the uproot.open method, taking a filename and outputting some ROOT object that has branches"""
    data = up.open(filename)
    if isinstance(data, HasBranches):
        return data
    else:
        raise TypeError("please read a TTree using the ':' separator.")


def get_POT(filename: str):
    """sums the POT value of each subrun from a"""
    data = up.open(filename)
    pots_per_subrun: aw.Array = data.get("ana/MetaTree").get("POT").array()  # type: ignore
    try:
        return aw.sum(pots_per_subrun)
    except:
        return -1.0
