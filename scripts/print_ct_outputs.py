#!/usr/bin/env python3
"""
Plot heat-maps of activity centered on the reconstructed neutrino position,
applying dead-wire-removal, island-finding and blip removal.

Usage: [uv run] python print_ct_outputs.py <path-to-ntuple> <event-index>
"""

import sys

import awkward as ak
import hist
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndi
import uproot as up

import sigmazerosearch.alg.island as island
from sigmazerosearch import loader
from sigmazerosearch.defaults import cthorpe_RHC_Tune_397
from sigmazerosearch.selection import signal_def
from sigmazerosearch.utils import WireGeometry, npfp

plt.style.use("plots/sigmazerosearch.tex.mplstyle")

PARTICLE_POS = False

figsizes = ((6, 6), (10, 4))

dead_wires = loader.load_dead_wire_map("data/MCC9_channel_list.txt")


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

    labelled = island._window_to_map(arr.ct_test_window_plane2[event], timebins, nwires)

    print(labelled[0])
    print("Islands:", labelled[1])

    _ = island.count_event_islands(arr, cthorpe_RHC_Tune_397)


def main(ntuple: str, nevent: int = 0, nwires: int = 100, timebins: int = 500):
    tree = up.open(ntuple + ":ana/OutputTree")

    arr = tree.arrays()

    # arr = arr[signal_def(arr)]

    h0 = (
        hist.Hist.new.Regular(
            timebins, 0, 7500, name="plane0", label="Plane 0 Time [tick]"
        )
        .Integer(0, nwires, name="wire", label="Relative Wire")
        .Double()
    )
    h1 = (
        hist.Hist.new.Regular(timebins, 1, 7500, name="plane1", label="Plane 1 [tick]")
        .Integer(0, nwires, name="wire", label="Relative Wire")
        .Double()
    )
    h2 = (
        hist.Hist.new.Regular(timebins, 1, 7500, name="plane2", label="Plane 2 [tick]")
        .Integer(0, nwires, name="wire", label="Relative Wire")
        .Double()
    )

    # arr = arr[arr.event == 361]

    event: int = nevent

    primary_vtx_wires = arr.ct_test_primary_vtx_wires[event]

    print(f"RSE({arr.run[event]}, {arr.subrun[event]}, {arr.event[event]})")
    print(
        f"neutrino_pdg = {arr.mc_nu_pdg[event]}\nhyperon_pdg = {arr.mc_hyperon_pdg[event]}"
    )
    print(f"decay_pdgs = {arr.mc_decay_pdg[event]}")
    print(f"n_showers = {npfp(arr, 'shower')[event]}")

    wires, times = ak.broadcast_arrays(
        ak.local_index(arr.ct_test_window_plane0[event], axis=0),
        arr.ct_test_window_plane0[event],
    )

    h0.fill_flattened(plane0=times, wire=wires)
    h1.fill_flattened(
        *ak.broadcast_arrays(
            arr.ct_test_window_plane1[event],
            ak.local_index(arr.ct_test_window_plane1[event], axis=0),
        )
    )
    h2.fill_flattened(
        *ak.broadcast_arrays(
            arr.ct_test_window_plane2[event],
            ak.local_index(arr.ct_test_window_plane2[event], axis=0),
        )
    )

    fig, axs = plt.subplots(
        ncols=1,
        nrows=3,
        figsize=(10, 12.5),
        layout="constrained",
        squeeze=True,
        sharex=True,
    )

    for hi, ax in zip([h0, h1, h2], axs):
        hi.plot(ax=ax, cbar=False)
        ax.axhline(
            50, color="lightgrey", linestyle="dashed", alpha=0.4
        )  # draw central line for reco vertex position.

        ax.yaxis.set_ticks(range(0, 101, 5), [str(x) for x in range(0, 101, 5)])
        ax.set_xlabel(f"{hi.axes[0].label} / {hi.axes[0].widths[0]:.1f} ticks")

    axs[0].set_title(f"RSE({arr.run[event]}, {arr.subrun[event]}, {arr.event[event]})")

    # plt.savefig(
    #     f"plots/ct/ct_heatmap_plane2_{arr.run[event]}_{arr.subrun[event]}_{arr.event[event]}.png",
    #     bbox_inches="tight",
    # )
    if PARTICLE_POS:
        lepton = [arr.mc_lepton_start_x, arr.mc_lepton_start_y, arr.mc_lepton_start_z]
        lamb = [
            arr.mc_sigmazero_lambda_end_x,
            arr.mc_sigmazero_lambda_end_y,
            arr.mc_sigmazero_lambda_end_z,
        ]
        nu = [arr.reco_primary_vtx_x, arr.reco_primary_vtx_y, arr.reco_primary_vtx_z]
        lepton_y_t = ak.zip(
            {
                "y": WireGeometry.pos_to_y(*lepton),
                "t": WireGeometry.pos_to_time(*lepton),
            }
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
        axs[0].scatter(
            lepton_y_t[event]["t"],
            WireGeometry.wire_time_to_window(
                lepton_y_t[event]["y"], WireGeometry.pos_to_y(*nu)[event], 100
            ),
            100,
            color="red",
            marker="+",
            label=r"true $\mu$ start",
        )
        axs[0].scatter(
            lambda_y_t[event]["t"],
            WireGeometry.wire_time_to_window(
                lambda_y_t[event]["y"], WireGeometry.pos_to_y(*nu)[event], 100
            ),
            100,
            color="red",
            marker="x",
            label=r"true $\Lambda$ end",
        )
        fig.legend()

    plt.savefig(
        f"plots/ct/ct_heatmap_all_{arr.run[event]}_{arr.subrun[event]}_{arr.event[event]}.png",
        bbox_inches="tight",
    )
    plt.show()

    # print(clipped)

    fig, axs = plt.subplots(
        ncols=1, nrows=3, figsize=(10, 12.5), layout="constrained", squeeze=True
    )

    for i, (hi, ax) in enumerate(zip([h0, h1, h2], axs)):
        print(f"Plane {i}: {np.ptp(primary_vtx_wires[i])}")
        clipped = np.clip(
            np.flip(hi.values().T, 0), 0, 1, dtype=np.int32, casting="unsafe"
        )
        clipped = island.remove_dead_window_channels(
            clipped, dead_wires[i], primary_vtx_wires[i]
        )
        clipped_labelled = ndi.label(clipped, island.LABEL_STRUCTURE)

        lengths = island.compute_island_sizes(clipped_labelled)

        # print(clipped_labelled[0])
        print("found", clipped_labelled[1], "features")

        im = ax.imshow(
            clipped_labelled[0],
            norm="linear",
            interpolation="none",
            extent=(0, 250, 0, 100),
        )

        fig.colorbar(im, ax=ax, label="Island Number")
        ax.axhline(50, color="lightgrey", linestyle="dashed", alpha=0.4)
        ax.set_ylabel("Relative Wire")
        ax.yaxis.set_ticks(range(0, 101, 25), [str(x) for x in range(0, 101, 25)])
        # ax.xaxis.set_ticks(range(0, 7501, 100), [str(x) for x in range(0, 7501, 100)])

        length_str = (
            ["Hits:"]
            + [f"{i + 1}: {length:.0f}" for i, length in enumerate(lengths)]  # if hits
            if isinstance(lengths[0], int)
            else ["Lengths:"]
            + [
                f"{i + 1}: {length:.2f}" for i, length in enumerate(lengths)
            ]  # if lengths
        )

        ax.text(
            0.95,
            0.95,
            "\n".join(length_str),
            ha="right",
            va="top",
            c="white",
            size="small",
            transform=ax.transAxes,
        )

    plt.savefig(
        f"plots/ct/ct_islands_{arr.run[event]}_{arr.subrun[event]}_{arr.event[event]}.png",
        bbox_inches="tight",
    )
    plt.show()

    min_length = 4.0

    lengths = ak.Array(lengths)
    ind = ak.local_index(lengths)[lengths > min_length] + 1

    fig, ax = plt.subplots(figsize=figsizes[1], layout="constrained")

    im = ax.imshow(
        np.where(
            np.isin(clipped_labelled[0], ak.to_numpy(ind)), clipped_labelled[0], 0
        ),
        norm="linear",
        interpolation="none",
    )

    ax.figure.colorbar(im, ax=ax, label="Island Number")
    ax.axhline(50, color="lightgrey", linestyle="dashed", alpha=0.4)
    # ax.set_title(
    #     f"Collection Plane Islands ($> {meta_arr.ConnectednessThreshold[0]}$ ADC)",
    #     loc="left",
    # )
    # ax.set_title(
    #     f"{len(ind)} features, RSE({arr.run[event]}, {arr.subrun[event]}, {arr.event[event]})"
    # )
    ax.set_ylabel("Relative Wire")
    ax.yaxis.set_ticks(range(0, 101, 25), [str(x) for x in range(0, 101, 25)])
    plt.savefig(
        f"plots/ct/ct_islands_trimmed_plane2_{arr.run[event]}_{arr.subrun[event]}_{arr.event[event]}.png",
        bbox_inches="tight",
    )

    plt.show()


def show_removal():
    fig, axs = plt.subplots(
        nrows=2, ncols=2, figsize=(10, 10), layout="constrained", squeeze=True
    )

    m = np.eye(100)
    axs[0].matshow(m)
    plt.show()


if __name__ == "__main__":
    # show_removal()
    if len(sys.argv) != 3:
        print("Please give two arguments")
        sys.exit(1)
    main(sys.argv[1], int(sys.argv[2]))
    # funky(sys.argv[1], int(sys.argv[2]))
