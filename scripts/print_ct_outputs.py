#!/usr/bin/env python3
import sys

import awkward as ak
import hist
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndi
import uproot as up

import sigmazerosearch.alg.island as island
from sigmazerosearch.defaults import cthorpe_RHC_Tune_397
from sigmazerosearch.selection import signal_def
from sigmazerosearch.utils import WireGeometry

plt.style.use("plots/sigmazerosearch.tex.mplstyle")


def funky(ntuple: str, nevent: int = 0, nwires: int = 100, timebins: int = 250):
    tree = up.open(ntuple + ":ana/OutputTree")

    arr = tree.arrays()

    arr = arr[signal_def(arr)]
    event: int = nevent

    print(f"RSE({arr.run[event]}, {arr.subrun[event]}, {arr.event[event]})")
    print(
        f"neutrino_pdg = {arr.mc_nu_pdg[event]}\nhyperon_pdg = {arr.mc_hyperon_pdg[event]}"
    )
    print(f"decay_pdgs = {arr.mc_decay_pdg[event]}")

    labelled = island._window_to_map(arr.ct_test_window_plane0[event], timebins, nwires)

    print(labelled[0])
    print("Islands:", labelled[1])

    _ = island.count_event_islands(arr, cthorpe_RHC_Tune_397)


def main(ntuple: str, nevent: int = 0, nwires: int = 100, timebins: int = 250):
    tree = up.open(ntuple + ":ana/OutputTree")
    meta_tree = up.open(ntuple + ":ana/MetaTree")

    meta_arr = meta_tree.arrays(["ConnectednessThreshold"])

    arr = tree.arrays()

    arr = arr[signal_def(arr)]

    h0 = (
        hist.Hist.new.Regular(
            timebins, 0, 7500, name="plane0", label="Plane 0 Time [tick]"
        )
        .Integer(0, nwires, name="wire", label="Relative Wire")
        .Double()
    )
    # h1 = hist.Hist.new.Regular(
    #     100, 1, 7500, name="plane1", label="Plane 1 [tick]"
    # ).Double()
    # h2 = hist.Hist.new.Regular(
    #     100, 1, 7500, name="plane2", label="Plane 2 [tick]"
    # ).Double()

    event: int = nevent

    print(f"RSE({arr.run[event]}, {arr.subrun[event]}, {arr.event[event]})")
    print(
        f"neutrino_pdg = {arr.mc_nu_pdg[event]}\nhyperon_pdg = {arr.mc_hyperon_pdg[event]}"
    )
    print(f"decay_pdgs = {arr.mc_decay_pdg[event]}")

    wires, times = ak.broadcast_arrays(
        ak.local_index(arr.ct_test_window_plane0[event], axis=0),
        arr.ct_test_window_plane0[event],
    )

    h0.fill_flattened(plane0=times, wire=wires)
    # h1.fill_flattened(plane1=ak.ravel(arr.ct_test_window_plane1))
    # h2.fill_flattened(plane2=ak.ravel(arr.ct_test_window_plane2))

    fig, ax = plt.subplots(figsize=(14, 10))

    h0.plot(ax=ax)
    # h1.plot(ax=ax)
    # h2.plot(ax=ax)
    ax.axhline(
        50, color="lightgrey", linestyle="dashed", alpha=0.4
    )  # draw central line for reco vertex position.

    ax.set_title("Collection Plane Heatmap", loc="left")
    ax.set_title(f"RSE({arr.run[event]}, {arr.subrun[event]}, {arr.event[event]})")
    ax.yaxis.set_ticks(range(0, 101, 5), [str(x) for x in range(0, 101, 5)])

    ax.set_xlabel(f"{h0.axes[0].label} / {h0.axes[0].widths[0]:.1f} ticks")

    plt.savefig("plots/ct_test_plane0.png", bbox_inches="tight")
    lepton = [arr.mc_lepton_start_x, arr.mc_lepton_start_y, arr.mc_lepton_start_z]
    lamb = [
        arr.mc_sigmazero_lambda_end_x,
        arr.mc_sigmazero_lambda_end_y,
        arr.mc_sigmazero_lambda_end_z,
    ]
    nu = [arr.reco_primary_vtx_x, arr.reco_primary_vtx_y, arr.reco_primary_vtx_z]
    lepton_y_t = ak.zip(
        {"y": WireGeometry.pos_to_y(*lepton), "t": WireGeometry.pos_to_time(*lepton)}
    )
    lambda_y_t = ak.zip(
        {"y": WireGeometry.pos_to_y(*lamb), "t": WireGeometry.pos_to_time(*lamb)}
    )

    print(
        lepton_y_t[event]["t"],
        WireGeometry.wire_time_to_window(
            lepton_y_t[event]["y"], WireGeometry.pos_to_y(*nu)[event], 100
        ),
    )
    ax.scatter(
        lepton_y_t[event]["t"],
        WireGeometry.wire_time_to_window(
            lepton_y_t[event]["y"], WireGeometry.pos_to_y(*nu)[event], 100
        ),
        100,
        color="red",
        marker="+",
        label=r"true $\mu$ start",
    )
    ax.scatter(
        lambda_y_t[event]["t"],
        WireGeometry.wire_time_to_window(
            lambda_y_t[event]["y"], WireGeometry.pos_to_y(*nu)[event], 100
        ),
        100,
        color="red",
        marker="x",
        label=r"true $\Lambda$ end",
    )
    ax.legend()

    plt.show()

    clipped = np.clip(h0.values(), 0, 1, dtype=np.int32, casting="unsafe")

    print(clipped)

    fig, ax = plt.subplots(figsize=(14, 10))

    clipped_labelled = ndi.label(clipped, island.LABEL_STRUCTURE)

    print(clipped_labelled[0])
    print("found", clipped_labelled[1], "features")

    im = ax.imshow(
        np.flip(clipped_labelled[0].T, 0),
        norm="linear",
        interpolation="none",
        extent=(0, 250, 0, 100),
    )

    ax.figure.colorbar(im, ax=ax, label="Island Number")
    ax.axhline(50, color="lightgrey", linestyle="dashed", alpha=0.4)
    ax.set_title(
        f"Collection Plane Islands ($> {meta_arr.ConnectednessThreshold[0]}$ ADC)",
        loc="left",
    )
    ax.set_title(
        f"{clipped_labelled[1]} features, RSE({arr.run[event]}, {arr.subrun[event]}, {arr.event[event]})"
    )
    ax.set_ylabel("Relative Wire")
    # ax.set_xlabel(f"{h0.axes[0].label} / {h0.axes[0].widths[0]:.1f} ticks")
    plt.show()


if __name__ == "__main__":
    # parser = argparse.ArgumentParser(__file__)
    # parser.add_argument("wires")
    # parser.parse_args()
    main(sys.argv[1], int(sys.argv[2]))
    # funky(sys.argv[1], int(sys.argv[2]))
