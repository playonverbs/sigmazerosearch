"""
Functions and definitions for a Fiducial Volume (FV) within the detector.

An FV aims to remove areas of suboptimal detector response and areas where
large amounts of 'Dirt' interactions occur.
"""

import awkward as ak
import numpy as np

__all__ = ["in_active_tpc"]

TPC_CENTER = (126.625, 0.97, 518.5)
TPC_SIDE_LENGTHS = (236.35, 233.0, 1036.8)

Sides = tuple[float, float]

FV_x: Sides = (0.0, 256.35)
FV_y: Sides = (-115.53, 117.47)
FV_z: Sides = (0.1, 1036.9)


def in_active_tpc(x, y, z):
    """
    Given a set of scalar or vector values corresponding to x, y, z
    coordinates, returns if the point(s) are inside of the fiducial volume as
    currently defined by this file.
    """
    return np.logical_and.reduce(
        (
            x >= FV_x[0],
            x <= FV_x[1],
            y >= FV_y[0],
            y <= FV_y[1],
            z >= FV_z[0],
            z <= FV_z[1],
        )
    )


def r_to_closest_wall(x, y, z):
    """
    Computes the distance from the given point to the closest TPC Fiducial
    Volume boundary as defined by the `FV_{x,y,z}` variables.

    The point is assumed to be in the detector coordinate frame and given in
    units of centimeters.
    """
    return ak.min(
        (
            ak.min((x - FV_x[0], FV_x[1] - x), axis=0),
            ak.min((y - FV_y[0], FV_y[1] - y), axis=0),
            ak.min((z - FV_z[0], FV_z[1] - z), axis=0),
        ),
        axis=0,
    )
