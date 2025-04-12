"""
Mixed utility functions.
"""

import logging
import pathlib
import sys

import awkward as ak
import numpy as np
import vector
from matplotlib.figure import Figure

from sigmazerosearch.general import Config

logger = logging.getLogger(__name__)


class WireGeometry:
    """
    A suite of methods to convert from MicroBooNE detector coordinates to
    wire,time coordinates for each anode plane.
    """

    A_W = 3.33328
    C_U = 338.140
    C_V = 2732.53
    C_Y = 4799.19
    A_T = 18.2148
    C_T = 818.351

    cos60 = 0.5
    sin60 = np.sqrt(3) / 2.0

    @classmethod
    def pos_to_u(cls, x, y, z):
        return cls.A_W * (-cls.sin60 * y + cls.cos60 * z) + cls.C_U

    @classmethod
    def pos_to_v(cls, x, y, z):
        return cls.A_W * (cls.sin60 * y + cls.cos60 * z) + cls.C_V

    @classmethod
    def pos_to_y(cls, x, y, z):
        return cls.A_W * z + cls.C_Y

    @classmethod
    def pos_to_time(cls, x, y, z):
        return cls.A_T * x + cls.C_T

    @classmethod
    def wire_time_to_window(cls, wire, nu_wire, window_width):
        window_origin = nu_wire - (window_width // 2)
        return wire - window_origin


def _save_plot(config: Config, fig: Figure, title: str):
    if isinstance(config.plot_format, list):
        for format in config.plot_format:
            fig.savefig(
                pathlib.Path(config.plot_dir / f"{title}.{format}"),
                dpi=300,
                bbox_inches="tight",
            )
            logger.info(
                "saved plot to %s", pathlib.Path(config.plot_dir / f"{title}.{format}")
            )
    else:
        fig.savefig(
            pathlib.Path(config.plot_dir / f"{title}.{config.plot_format}"),
            dpi=300,
            bbox_inches="tight",
        )
        logger.info(
            "saved plot to %s",
            pathlib.Path(config.plot_dir / f"{title}.{config.plot_format}"),
        )


def npfp(arr, opt: str | None = None) -> ak.Array:
    """
    Returns the number of pfps in the given array that are either track-like,
    shower-like, or both.
    """
    trk = arr["pfp_trk_shr_score"] > 0.5
    shr = arr["pfp_trk_shr_score"] < 0.5

    if opt == "both" or opt is None:
        return ak.sum((trk | shr), axis=1)
    elif "track" in opt:
        return ak.sum(trk, axis=1)
    elif "shower" in opt:
        return ak.sum(shr, axis=1)
    else:
        raise TypeError


def displacement(arr, x_i, y_i, z_i) -> ak.Array:
    """
    Compute displacement array from given `x,y,z` array indices to
    `reco_primary_vtx` equivalents.
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

    return ak.mask(res, ak.num(res) != 0)  # type: ignore


def separation(
    arr_i: ak.Array,
    arr_j: ak.Array,
    index: str,
    suffixes: list[str] = ["x", "y", "z"],
) -> ak.Array:
    """
    Compute the separation between two <inv:#ak.Array> subsets using the
    partial field name `index` + `suffixes`.
    """

    v = vector.zip(
        {
            "x": arr_i[index + suffixes[0]],
            "y": arr_i[index + suffixes[1]],
            "z": arr_i[index + suffixes[2]],
        }
    )

    u = vector.zip(
        {
            "x": arr_j[index + suffixes[0]],
            "y": arr_j[index + suffixes[1]],
            "z": arr_j[index + suffixes[2]],
        }
    )

    res = (v - u).mag

    return res.mask[ak.num(res) != 0]


def filter_by_rse(arr: ak.Array, run: int, subrun: int, event: int) -> ak.Array:
    """
    Takes an array and either three numbers corresponding to run, subrun, event
    numbers or a list of such numbers
    """
    cond = np.logical_and.reduce(
        (
            arr["run"] == run,
            arr["subrun"] == subrun,
            arr["event"] == event,
        )
    )
    return arr[cond]  # type: ignore


def print_rse(arr: ak.Array, file=sys.stdout):
    """
    Print the event details in the format `run subrun event`, the `file`
    parameter can be used to pipe this output to a file (useful for filtering
    upstream).
    """
    for elem in arr:
        print(f"{elem.run} {elem.subrun} {elem.event}", file=file)  # type: ignore
    logger.info(f"output rse numbers for {len(arr)} events to {file.name}")


def file_ok(filename: str, mode: str = "read") -> bool:
    """
    Test file for read/write availability and convert Exceptions to boolean.

    Defaults to checking for read availability.
    """
    try:
        fp = open(filename) if mode == "read" else open(filename, "w")
    except FileNotFoundError:
        logger.warning(f"file {filename} not found")
        return False
    else:
        fp.close()

    return True
