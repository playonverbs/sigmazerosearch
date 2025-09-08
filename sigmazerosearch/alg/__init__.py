"""
Methods that operate on physics information to aid in selections.

Many of the methods here accept
{py:obj}`sigmazerosearch.selection.ParameterSet` as an argument which contains
values that affect the physics decisions of this framework.

Submodules are organised by what problems they aim to solve:

- `fv`: provides a fiducial volume cut and some other definitions regarding the
  detector geometry.
- `muon`: provides methods of selecting a muon-like object.
- `lamb`: methods for selecting the $\\Lambda \\rightarrow p + \\pi^-$ subsystem.
"""

import awkward as ak
import vector as vec


def p_to_p3D(arr: ak.Array, prefix: str):
    """
    Given an input <inv:#ak.Array> and a prefix string, computes the 3D
    momentum vector for each `{prefix}` sub-object.

    This assumes for a prefix string `{prefix}` the fields `{prefix}_mom`,
    `{prefix}_start_{x,y,z}` and `{prefix}_end_{x,y,z}` exist -- corresponding
    to the scalar momentum, start position and end position components. Uses
    the range $\\vec{s} = \\vec{x}_\\mathrm{end} - \\vec{x}_\\mathrm{start}$ as
    a proxy for direction.

    :param ak.Array arr:
        An array matching the format of an input ntuple.
    :param str prefix:
        A string matching the beginning of the position and momentum fields of
        the object. To use a set of fields like `mc_decay_start_x`, provide a
        prefix of `mc_decay`.
    """
    mom = arr[prefix + "_mom"]

    output = vec.zip(
        {
            "px": arr[f"{prefix}_end_x"] - arr[f"{prefix}_start_x"],
            "py": arr[f"{prefix}_end_y"] - arr[f"{prefix}_start_y"],
            "pz": arr[f"{prefix}_end_z"] - arr[f"{prefix}_start_z"],
        }
    )

    return output.unit().scale3D(mom)


def p_to_p4D(arr: ak.Array, prefix: str, mass=None, energy=None):
    """
    Given an input <inv:#ak.Array>, a prefix string and either a mass _or_
    energy parameter, computes the 4D momentum vector for each `{prefix}`
    sub-object.

    See <project:#p_to_p3D> for information on the `arr` and `prefix` variables.

    :param mass:
        An array-like object with the same shape as the `{prefix}_` fields in
        `arr`. Represents the mass of each object in units of GeV.

        Mutually exclusive with the `energy` parameter.
    :param energy:
        An array-like object with the same shape as the `{prefix}_` fields in
        `arr`. Represents the energy of each object in units of GeV.

        Mutually exclusive with the `mass` parameter.
    """
    assert (mass is None) != (
        energy is None
    ), "Please provide either a mass or energy parameter"

    output = p_to_p3D(arr, prefix)

    if mass is not None:
        return output.to_pxpypzmass(mass=mass)
    elif energy is not None:
        return output.to_pxpypzenergy(energy=energy)
    else:
        raise ValueError
