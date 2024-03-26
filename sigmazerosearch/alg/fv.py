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
    return np.logical_or.reduce(
        (x > FV_x[1], x < FV_x[0], y > FV_y[1], y < FV_y[0], z > FV_z[1], z < FV_z[0])
    )
