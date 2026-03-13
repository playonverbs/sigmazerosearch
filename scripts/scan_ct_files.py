import re
from pathlib import Path

import awkward as ak
import hist
import matplotlib.pyplot as plt
import numpy as np
import uproot as up

import sigmazerosearch.alg.island as island

plt.style.use("plots/sigmazerosearch.tex.mplstyle")
pat = re.compile(r".*_(\d+.\d+)_(\d+)w.*.root")
inverted_cmap = plt.colormaps.get("viridis")

nevent: int = 0


def scan_both():
    files = sorted(
        list(
            Path("/home/niam/phd/data/hyperons/multiSlice/run3b_RHC/cttest/scan").glob(
                "*.root"
            )
        )
    )

    isles = {"threshold": [], "window": [], "islands": []}

    fig, axs = plt.subplots(
        layout="constrained",
        nrows=5,
        ncols=5,
        sharex="col",
        sharey="row",
        figsize=(14, 10),
        gridspec_kw=dict(hspace=0.1, wspace=0.075),
    )

    axs = axs.T

    rse = (0, 0, 0)

    for file, ax in zip(files, axs.reshape(-1)):
        threshold, window = pat.findall(file.name)[0]
        threshold, window = float(threshold), int(window)

        print(f"Minimum threshold: {threshold} ADC, Wire window: {window} wires")

        nt = up.open(str(file) + ":ana/OutputTree")
        mt = up.open(str(file) + ":ana/MetaTree")
        arr = nt.arrays()
        marr = mt.arrays()
        assert threshold == marr.ConnectednessThreshold[0]
        assert window == marr.ConnectednessWindowW[0]

        islands = island._window_to_map(arr[nevent].ct_test_window_plane1, 500, window)[
            1
        ]
        isles["threshold"].append(threshold)
        isles["window"].append(window)
        isles["islands"].append(islands)

        timebins = 250
        nwires = window
        h0 = (
            hist.Hist.new.Regular(
                timebins, 0, 7500, name="plane1", label="Plane 1 Time [tick]"
            )
            .Integer(0, nwires, name="wire", label="Relative Wire")
            .Double()
        )

        wires, times = ak.broadcast_arrays(
            ak.local_index(arr.ct_test_window_plane1[nevent], axis=0),
            arr.ct_test_window_plane1[nevent],
        )

        h0.fill_flattened(plane1=times, wire=wires)

        # fig, ax = plt.subplots(figsize=(14, 10))

        h0.plot(ax=ax, cbar=False)
        ax.axhline(
            window / 2, color="lightgrey", linestyle="dashed", alpha=0.4
        )  # draw central line for reco vertex position.
        ax.set_xlabel(None)
        ax.set_ylabel(None)
        ax.yaxis.set_ticks([])

        ax.text(
            0.95,
            0.95,
            f"{threshold} ADC\n{window} wires\n{islands} islands",
            transform=ax.transAxes,
            ha="right",
            va="top",
            size="small",
            c="white",
        )

        # ax.set_title(f"Collection Plane ({threshold} ADC, {window} wires)", loc="left")
        ax.yaxis.set_ticks(
            range(0, window + 1, 25), [str(x) for x in range(0, window + 1, 25)]
        )
        rse = (arr.run[nevent], arr.subrun[nevent], arr.event[nevent])

    fig.supxlabel("Time ticks")
    fig.supylabel("Relative Wire Numbers")
    # fig.suptitle(f"Collection Plane {rse}", x=0.125, y=0.905, ha="left")
    fig.suptitle(f"Collection Plane {rse}")
    plt.savefig(
        f"plots/ct_scanned_heatmaps_{rse[0]}-{rse[1]}-{rse[2]}.png", bbox_inches="tight"
    )
    plt.show()


def scan_adc():
    files = sorted(
        list(
            Path(
                "/home/niam/phd/data/hyperons/multiSlice/run3b_RHC/cttest/scan_adc_only"
            ).glob("*.root")
        )
    )

    # files = sorted(files, key=lambda x: int(re.search(r'\d+', str(x)).group()))

    isles = {"threshold": [], "window": [], "islands": []}

    fig, axs = plt.subplots(
        layout="constrained",
        nrows=len(files),
        sharex="col",
        sharey="row",
        figsize=(7, 10),
        # gridspec_kw=dict(hspace=0.1, wspace=0.075),
    )

    rse = (0, 0, 0)

    for file, ax in zip(files, axs.reshape(-1)):
        nt = up.open(str(file) + ":ana/OutputTree")
        mt = up.open(str(file) + ":ana/MetaTree")
        arr = nt.arrays()
        marr = mt.arrays()

        threshold = marr.ConnectednessThreshold[0]
        window = marr.ConnectednessWindowW[0]
        print(f"Minimum threshold: {threshold} ADC, Wire window: {window} wires")

        pic, islands = island._window_to_map(
            arr[nevent].ct_test_window_plane1, 500, window
        )
        isles["threshold"].append(threshold)
        isles["window"].append(window)
        isles["islands"].append(islands)

        timebins = 250
        nwires = window
        h0 = (
            hist.Hist.new.Regular(
                timebins, 0, 7500, name="plane1", label="Plane 0 Time [tick]"
            )
            .Integer(0, nwires, name="wire", label="Relative Wire")
            .Double()
        )

        wires, times = ak.broadcast_arrays(
            ak.local_index(arr.ct_test_window_plane1[nevent], axis=0),
            arr.ct_test_window_plane1[nevent],
        )

        h0.fill_flattened(plane1=times, wire=wires)

        # fig, ax = plt.subplots(figsize=(14, 10))

        ax.matshow(np.flip(pic.T, 0))
        # ax.imshow(
        #     np.flip(pic.T, 0),
        #     norm="linear",
        #     interpolation="none",
        #     # extent=(0, 250, 0, 100),
        # )
        # h0.plot(ax=ax, cbar=False)
        ax.axhline(
            window / 2, color="lightgrey", linestyle="dashed", alpha=0.4
        )  # draw central line for reco vertex position.
        ax.set_xlabel(None)
        ax.set_ylabel(None)
        ax.yaxis.set_ticks([])

        ax.text(
            0.95,
            0.95,
            f"{threshold} ADC\n{window} wires\n{islands} islands",
            transform=ax.transAxes,
            ha="right",
            va="top",
            size="small",
            c="white",
        )

        # ax.set_title(f"Collection Plane ({threshold} ADC, {window} wires)", loc="left")
        ax.yaxis.set_ticks(
            range(0, window + 1, 25), [str(x) for x in range(0, window + 1, 25)]
        )
        rse = (arr.run[nevent], arr.subrun[nevent], arr.event[nevent])

    fig.supxlabel("Time ticks")
    fig.supylabel("Relative Wire Numbers")
    # fig.suptitle(f"Collection Plane {rse}", x=0.125, y=0.905, ha="left")
    fig.suptitle(f"plane1 {rse}")
    # plt.savefig(
    #     f"plots/ct_scanned_heatmaps_{rse[0]}-{rse[1]}-{rse[2]}.png", bbox_inches="tight"
    # )
    plt.show()


if __name__ == "__main__":
    # scan_both()
    scan_adc()

# a = ak.Array(isles)

# h = (
#     hist.Hist.new.Regular(10, 0, 6, name="threshold")
#     .IntCategory([100], growth=True, name="window")
#     .Integer(0, 10, growth=True, name="islands")
#     .Double()
# )

# h.fill(
#     threshold=a.threshold,
#     window=a.window,
#     islands=a.islands
# )

# h.sort("window").project("window", "islands").plot()

# plt.show()
