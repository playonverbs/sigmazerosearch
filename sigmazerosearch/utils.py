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
