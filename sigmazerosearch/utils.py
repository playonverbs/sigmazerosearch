import logging
import sys

import awkward as ak
import vector


def displacement(x_i, y_i, z_i, arr) -> ak.Array:
    """
    compute displacement array from given x,y,z array indices to reco_primary_vtx
    """
    v = vector.zip(
        {
            "x": arr["reco_primary_vtx_x"],
            "y": arr["reco_primary_vtx_y"],
            "z": arr["reco_primary_vtx_z"],
        }
    )
    u = vector.zip(
        {
            "x": arr[x_i],
            "y": arr[y_i],
            "z": arr[z_i],
        }
    )

    res = (v - u).mag

    return ak.mask(res, ak.num(res) != 0)


def print_rse(arr: ak.Array, file=sys.stdout):
    """
    print the event details in the format "run subrun event",
    the file parameter can be used to pipe this output to a file
    (useful for filtering upstream).
    """
    for elem in arr:
        print(f"{elem.run} {elem.subrun} {elem.event}", file=file)  # type: ignore
    logging.info(f"output rse numbers for {len(arr)} events")


def file_ok(filename: str, mode: str = "read") -> bool:
    """
    Test file for read/write availability and convert Exceptions to boolean
    Defaults to checking for read availability
    """
    try:
        fp = open(filename) if mode == "read" else open(filename, "w")
    except FileNotFoundError:
        return False
    else:
        fp.close()

    return True
