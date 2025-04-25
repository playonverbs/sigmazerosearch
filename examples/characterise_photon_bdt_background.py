#!/usr/bin/env python3
from pathlib import Path

import awkward as ak
import hist
import matplotlib.pyplot as plt

from sigmazerosearch.alg.muon import select_mu_candidate
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
from sigmazerosearch.truth import GenEventType
from sigmazerosearch.utils import _save_plot, npfp

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
        "mc_ccnc",
        "mc_mode",
        "mc_lepton_pdg",
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
        "pfp_photon_bdt_score",
    ],
)

sel = Selection(
    config=conf,
    cuts=[
        Cut.null(),
        Cut(
            "has-flash-match",
            lambda arr: arr["flash_match_nu_slice_ID"] != -1,
        ),
        Cut(
            "fv",
            lambda arr: arr["reco_primary_vtx_inFV"],
        ),
        Cut(
            "showers",
            lambda arr: ak.sum(arr["pfp_trk_shr_score"] < 0.5, axis=1) >= 1,
        ),
        Cut(
            "tracks",
            lambda arr: ak.sum(arr["pfp_trk_shr_score"] >= 0.5, axis=1) >= 3,
        ),
        Cut(
            "muon-id",
            lambda arr: ak.sum(select_mu_candidate(arr, pset), axis=1) >= 1,
        ),
        Cut(
            "photon-bdt",
            # NOTE: do we need to require it
            lambda arr: ak.sum(arr.pfp_photon_bdt_score >= -0.0179, axis=1) >= 1,
            # lambda arr: ak.sum(
            #     np.logical_and.reduce(
            #         [
            #             arr.pfp_photon_bdt_score >= -0.0179,
            #             arr.pfp_trk_shr_score >= 0.5,
            #         ]
            #     ),
            #     axis=1,
            # )
            # >= 1,
        ),
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
            # "/home/niam/phd/data/hyperons/multiSlice/run3b_RHC/analysisOutputRHC_mSlice_Background_all.root",
            "/home/niam/phd/data/hyperons/multiSlice/run3b_RHC/cttest/analysisOutputRHC_Dirt_prodgenie_numi_uboone_overlay_dirt_rhc_mcc9_run3b_v28_all.root",
            # "/home/niam/phd/data/hyperons/multiSlice/run3b_RHC/cttest/analysisOutputRHC_mSlice_Background_prodgenie_numi_uboone_overlay_rhc_mcc9_run3b.root",
            SampleType.Dirt,
            None,
        ),
        # Sample(
        #     "EXT_all",
        #     "/home/niam/phd/data/hyperons/multiSlice/run3b_RHC/cttest/analysisOutputRHC_EXT_prod_extnumi_mcc9_v08_00_00_45_run3_run3b_reco2_all_reco2.root",
        #     SampleType.EXT,
        #     None
        # ),
        base_dir="/home/niam/phd/data/hyperons/run3b_RHC",
        target_POT=4.9e20,
    ),
    params=pset,
    label="nominal",
)

# logging.getLogger().setLevel(logging.INFO)

# logger = logging.getLogger("sigmazerosearch")
# logger.setLevel(logging.INFO)
# logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
# logging.info("Opening files")
sel.samples.as_table()
sel.open_files()

# logging.info("Applying cuts")

array = sel.apply_cut(sel.cuts, accumulate=True)

sel.cut_summary(header=True, format="md")

h = (
    hist.Hist.new.Integer(3, 10, name="ntracks", label=r"\# Tracks")
    .Integer(1, 10, name="nshowers", label=r"\# Showers")
    .StrCategory([t.name for t in GenEventType], name="mode", label="Mode")
    .StrCategory(["CC", "NC"], name="ccnc", label="CC/NC")
    .IntCategory([-12, -14, -16, 12, 14, 16], name="nupdg", label=r"$\nu$ type")
    .IntCategory(
        [-12, -14, -16, 12, 14, 16, -11, -13, -15, 11, 13, 15],
        name="leptonpdg",
        label="Lepton type",
    )
    .StrCategory([x.name for x in SampleType], name="sampletype", label="Sample")
    .Bool(name="signal")
    .Weight()
)

h.fill(
    ntracks=npfp(array, "track"),
    nshowers=npfp(array, "shower"),
    mode=array.mc_mode,
    ccnc=array.mc_ccnc,
    nupdg=array.mc_nu_pdg,
    leptonpdg=array.mc_lepton_pdg,
    sampletype=array.sample,
    signal=signal_def(array),
    weight=array.weight,
)


def topo_plots():
    fig, ax = plt.subplots()

    hsb: hist.BaseHist = h[{"signal": False}].project("ntracks")

    hsb.plot(ax=ax, label="Background")
    h[{"signal": True}].project("ntracks").plot(ax=ax, label="Signal")

    ax.set_title("NuMI Run 3b")
    ax.legend()
    plt.show()

    fig, ax = plt.subplots()

    hsb = h[{"signal": False}].project("nshowers")

    hsb.plot(ax=ax, label="Background")
    h[{"signal": True}].project("nshowers").plot(ax=ax, label="Signal")

    ax.set_title("NuMI Run 3b")
    ax.legend()
    plt.show()


fig, ax = plt.subplots(layout="constrained", figsize=(11, 5))

h[{"signal": False}].project("ccnc", "mode").stack("ccnc").plot(
    ax=ax, stack=True, histtype="fill"
)
h[{"signal": False}].project("mode").plot(ax=ax, histtype="band")

ax.set_title("NuMI Run 3b")
ax.legend()
_save_plot(sel.config, fig, "background_sel_after_photon_bdt_mode")
plt.show()

fig, ax = plt.subplots(layout="constrained", figsize=(11, 5))

h[{"signal": False}].project("mode", "ntracks").plot(ax=ax)
# hist_exp_label(ax, {"NuMI RHC": 4.9e20}, "Simulation", c="white")

ax.set_title("NuMI Run 3b")
_save_plot(sel.config, fig, "background_sel_after_photon_bdt_mode_tracks")
plt.show()

fig, ax = plt.subplots(layout="constrained", figsize=(11, 5))

h[{"signal": False}].project("mode", "nshowers").plot(ax=ax)
# hist_exp_label(ax, {"NuMI RHC": 4.9e20}, "Simulation", c="white")

ax.set_title("NuMI Run 3b")
_save_plot(sel.config, fig, "background_sel_after_photon_bdt_mode_showers")
plt.show()


if __name__ == "__main__":
    topo_plots()
