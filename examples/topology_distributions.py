#!/usr/bin/env python3
from pathlib import Path

import hist
import matplotlib.pyplot as plt

import sigmazerosearch.utils as utils
from sigmazerosearch.general import Config
from sigmazerosearch.selection import (
    Cut,
    ParameterSet,
    Sample,
    SampleSet,
    SampleType,
    Selection,
    signal_def,
)

plt.style.use("plots/sigmazerosearch.tex.mplstyle")

pset = ParameterSet(
    pid_cut=0.6,
    min_length=10,
    max_separation=1,
    proton_pid_cut=0.2,
    pion_pid_cut=0.2,
    separation_cut=3,
    w_lambda_min=1.1,
    w_lambda_max=1.20,
)

conf = Config(
    plot_save=True,
    plot_dir=Path("plots/"),
    plot_format=["png", "svg"],
    perf_percent=True,
    iterate=True,
    # iterate_step="50 MB",
    iterate_step="40 MB",
    branch_list=[
        "flash_match_nu_slice_ID",
        "reco_primary_vtx_inFV",
        "reco_primary_vtx_x",
        "reco_primary_vtx_y",
        "reco_primary_vtx_z",
        "pfp_trk_shr_score",
        "mc_nu_pdg",
        "mc_hyperon_pdg",
        "mc_nu_pos_x",
        "mc_nu_pos_y",
        "mc_nu_pos_z",
        "mc_decay_pdg",
        "trk_llrpid",
        "trk_length",
        "trk_start_x",
        "trk_start_y",
        "trk_start_z",
    ],
)

sel = Selection(
    config=conf,
    cuts=[
        Cut(
            "has-flash-match",
            lambda arr: arr["flash_match_nu_slice_ID"] != -1,
        ),
        Cut(
            "fv",
            lambda arr: arr["reco_primary_vtx_inFV"],
        ),
        # Cut(
        #     "showers",
        #     lambda arr: ak.sum(arr["pfp_trk_shr_score"] < 0.5, axis=1) >= 1,
        # ),
        # Cut(
        #     "tracks",
        #     lambda arr: ak.sum(arr["pfp_trk_shr_score"] >= 0.5, axis=1) >= 3,
        # ),
        # Cut(
        #     "muon-id",
        #     lambda arr: ak.sum(select_mu_candidate(arr, pset), axis=1) >= 1,
        # ),
        # Cut(
        #     "ct-test",
        #     lambda arr: island.count_event_islands(arr, pset).ct_test_islands_plane0 <= 10,
        # )
        # Cut(
        #     "photon-bdt",
        #     # NOTE: do we need to require it
        #     lambda arr: ak.sum(arr.pfp_photon_bdt_score >= -0.0179, axis=1) >= 1,
        #     # lambda arr: ak.sum(
        #     #     np.logical_and.reduce(
        #     #         [
        #     #             arr.pfp_photon_bdt_score >= -0.0179,
        #     #             arr.pfp_trk_shr_score >= 0.5,
        #     #         ]
        #     #     ),
        #     #     axis=1,
        #     # )
        #     # >= 1,
        # ),
        # Cut("lambda-select-box", lambda arr: select_p_pi_candidates_box(arr, pset)),
        # Cut(
        #     "has-photon",
        #     lambda arr: np.logical_and.reduce([
        #         ak.sum(arr["trk_llrpid"] >= 0.0, axis=1) >= 1
        #     ])
        # )
    ],
    samples=SampleSet(
        Sample(
            "hyperon",
            # "/home/niam/phd/data/hyperons/multiSlice/run3b_RHC/analysisOutputRHC_mSlice_cthorpe_make_hyperon_events_numi_rhc_run3b_hyperon_reco2_reco2.root",
            "/home/niam/phd/data/hyperons/multiSlice/run3b_RHC/cttest/analysisOutputRHC_mSlice_cttest_nuslice_Hyperon_cthorpe_make_hyperon_events_numi_run3b_hyperon_reco2_reco2.root",
            # "/home/niam/phd/data/hyperons/multiSlice/run3b_RHC/cttest/analysisOutputRHC_mSlice_Hyperon_cthorpe_make_hyperon_events_numi_run3b_hyperon_reco2_reco2.root",
            SampleType.Hyperon,
            None,
        ),
        Sample(
            "background_all",
            # "/home/niam/phd/data/hyperons/multiSlice/run3b_RHC/analysisOutputRHC_mSlice_Background_all.root",
            "/home/niam/phd/data/hyperons/multiSlice/run3b_RHC/cttest/analysisOutputRHC_mSlice_cttest_nuslice_Background_all_prodgenie_numi_uboone_overlay_rhc_mcc9_run3b.root",
            # "/home/niam/phd/data/hyperons/multiSlice/run3b_RHC/cttest/analysisOutputRHC_mSlice_Background_prodgenie_numi_uboone_overlay_rhc_mcc9_run3b.root",
            SampleType.Background,
            None,
        ),
        Sample(
            "dirt_all",
            "/home/niam/phd/data/hyperons/multiSlice/run3b_RHC/cttest/analysisOutputRHC_Dirt_prodgenie_numi_uboone_overlay_dirt_rhc_mcc9_run3b_v28_all.root",
            SampleType.Dirt,
            None,
        ),
        base_dir="/home/niam/phd/data/hyperons/run3b_RHC",
        target_POT=4.9e20,
    ),
    params=pset,
    label="nominal",
)

sel.samples.as_table()
sel.open_files()

arr = sel.apply_cut(sel.cuts, accumulate=True)

sel.cut_summary(header=True, format="md")

h = (
    hist.Hist.new.Integer(0, 6, name="ntracks", label="\# Tracks")
    .Integer(0, 6, name="nshowers", label="\# Showers")
    .StrCategory([x.name for x in SampleType], name="sampletype", label="Sample")
    .IntCategory([3212], name="hyperonpdg", label="Hyperon PDG")
    .Bool(name="signal", label="Signal")
    .Weight()
)

h.fill(
    ntracks=utils.npfp(arr, "track"),
    nshowers=utils.npfp(arr, "shower"),
    sampletype=arr["sample"],
    hyperonpdg=arr["mc_hyperon_pdg"],
    signal=signal_def(arr),
    weight=arr["weight"],
)

h_signal = h[{"signal": 1}]
h_background = h[{"signal": 0}]


def two_dimensional_plots():
    fig, ax = plt.subplots()

    h_signal.project("ntracks", "nshowers").plot2d(ax=ax)

    ax.set_title("NuMI Run 3b (Signal)")

    utils._save_plot(sel.config, fig, "topology_distributions_signal")

    plt.show()

    fig, ax = plt.subplots()

    h_background.project("ntracks", "nshowers").plot2d(ax=ax)

    ax.set_title("NuMI Run 3b (Background)")

    utils._save_plot(sel.config, fig, "topology_distributions_background")

    plt.show()


def track_plots():
    fig, ax = plt.subplots(layout="constrained")
    h_signal.project("ntracks", "sampletype").stack("sampletype").plot(
        ax=ax,
        histtype="fill",
        stack=True,  # , lw=0.5, ec="grey"
    )
    h_signal.project("ntracks").plot(ax=ax, histtype="band")

    ax.set_title("NuMI Run 3b (Signal)")
    ax.set_ylabel("Events")

    ax.legend(loc="upper right")
    # fig.legend(loc="outside upper left", ncols=5)
    utils._save_plot(sel.config, fig, "topology_tracks_signal")
    plt.show()

    fig, ax = plt.subplots(layout="constrained")

    h_background.project("ntracks", "sampletype").stack("sampletype").plot(
        ax=ax, histtype="fill", stack=True
    )
    h_background.project("ntracks").plot(ax=ax, histtype="band")

    ax.set_title("NuMI Run 3b (Background)")
    ax.set_ylabel("Events")
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)

    # ax.legend()
    ax.legend(loc="upper right")
    utils._save_plot(sel.config, fig, "topology_tracks_background")
    plt.show()


def shower_plots():
    fig, ax = plt.subplots(layout="constrained")
    h_signal.project("nshowers", "sampletype").stack("sampletype").plot(
        ax=ax, histtype="fill", stack=True
    )
    h_signal.project("nshowers").plot(ax=ax, histtype="band")

    ax.set_title("NuMI Run 3b (Signal)")
    ax.set_ylabel("Events")

    # ax.legend()
    ax.legend(loc="upper right")
    utils._save_plot(sel.config, fig, "topology_showers_signal")
    plt.show()

    fig, ax = plt.subplots(layout="constrained")

    h_background.project("nshowers", "sampletype").stack("sampletype").plot(
        ax=ax,
        histtype="fill",
        stack=True,
    )
    h_background.project("nshowers").plot(ax=ax, histtype="band")

    ax.set_title("NuMI Run 3b (Background)")
    ax.set_ylabel("Events")
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)

    # ax.legend()
    ax.legend(loc="upper right")
    utils._save_plot(sel.config, fig, "topology_showers_background")
    plt.show()


if __name__ == "__main__":
    two_dimensional_plots()
    track_plots()
    shower_plots()
