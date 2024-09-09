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


def _save_plot(config: Config, fig: Figure, title: str):
    if isinstance(config.plot_format, list):
        for format in config.plot_format:
            fig.savefig(
                pathlib.Path(config.plot_dir / f"{title}.{format}"),
                dpi=300,
                bbox_inches="tight",
            )
            logging.info(
                "saved plot to %s", pathlib.Path(config.plot_dir / f"{title}.{format}")
            )
    else:
        fig.savefig(
            pathlib.Path(config.plot_dir / f"{title}.{config.plot_format}"),
            dpi=300,
            bbox_inches="tight",
        )
        logging.info(
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
    Compute displacement array from given x,y,z array indices to
    reco_primary_vtx
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
    Print the event details in the format "run subrun event",
    the file parameter can be used to pipe this output to a file
    (useful for filtering upstream).
    """
    for elem in arr:
        print(f"{elem.run} {elem.subrun} {elem.event}", file=file)  # type: ignore
    logging.info(f"output rse numbers for {len(arr)} events to {file.name}")


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
