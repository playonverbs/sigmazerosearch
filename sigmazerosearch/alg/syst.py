"""
Contains functions relevant for computing systematics, including those
beamline-related weights specific to NuMI analyses.
"""

from typing import Literal

import numpy as np
import vector as vec
from numpy.typing import ArrayLike
from scipy.spatial.transform import Rotation as Rot

ROTATION_MATRIX: Rot = Rot.from_matrix(
    np.array(
        [
            [0.92103853804025681562, 0.022713504803924120662, 0.38880857519374290021],
            [
                4.6254001262154668408e-05,
                0.99829162468141474651,
                -0.058427989452906302359,
            ],
            [-0.38947144863934973769, 0.053832413938664107345, 0.91946400794392302291],
        ],
    ).T
)
"""
Describes the rotation from the detector coordinate system to that of the NuMI
beamline.
"""


def get_numi_angle(
    px: ArrayLike,
    py: ArrayLike,
    pz: ArrayLike,
    direction: Literal["beam", "target"] = "beam",
):
    """
    Takes the direction vector components of a neutrino in the detector
    frame and produces the angle between it and the NuMI beamline, or the
    target (if specified).

    A co-ordinate transform is carried out from the detector coordinate system
    to the NuMI beamline frame.
    """
    det_xyz: vec.MomentumObject3D = vec.zip({"px": px, "py": py, "pz": pz})

    beam_coords = det_xyz.rotate_quaternion(
        *ROTATION_MATRIX.as_quat(canonical=True, scalar_first=True)
    )

    beamdir = vec.obj(x=0, y=0, z=1).unit()

    if direction == "target":
        beamdir = vec.obj(x=5502, y=7259, z=67270).unit()

    angle = beam_coords.deltaangle(beamdir) * (180.0 / np.pi)

    return angle
