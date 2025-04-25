#!/usr/bin/env python3
from pathlib import Path

import awkward as ak
import hist
import matplotlib.pyplot as plt
import numpy as np
import particle as part

import sigmazerosearch.utils as utils
from sigmazerosearch.general import PDG, Config
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
        "mc_ccnc",
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
        "pfp_true_pdg",
        "pfp_true_length",
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
    hist.Hist.new.Regular(50, 0, 700, name="trklength", label="Length [cm]", flow=True)
    .StrCategory([x.name for x in SampleType], name="sampletype", label="Sample")
    .IntCategory([1] + list(PDG), name="truepdg", label="True PDG")
    .Bool(name="signal", label="Signal")
    .Weight()
)

arr = arr[(np.abs(arr.mc_nu_pdg) == 14) & (arr.mc_ccnc == "CC")]

trk_cond = arr.pfp_trk_shr_score >= 0.5

# index only the longest track, apply to fields later
maxl = ak.from_regular(
    ak.drop_none(ak.argmax(arr["trk_length"][trk_cond], axis=-1, keepdims=True))
)

barr = ak.broadcast_arrays(
    arr["trk_length"][trk_cond][maxl],
    arr["sample"],
    np.abs(arr["pfp_true_pdg"])[trk_cond][maxl],
    signal_def(arr),
    arr["weight"],
)

h.fill(
    trklength=ak.flatten(barr[0]),
    sampletype=ak.flatten(barr[1]),
    truepdg=ak.flatten(barr[2]),
    signal=ak.flatten(barr[3]),
    weight=ak.flatten(barr[4]),
)

h_signal = h[{"signal": 1}]
h_background = h[{"signal": 0}]

fig, ax = plt.subplots(layout="constrained")

hs = h.project("trklength", "truepdg").stack("truepdg")
hs = hist.Stack.from_iter(filter(lambda s: not s.empty(flow=True), hs))

for i in hs:
    i.name = (
        f"${part.Particle.from_pdgid(i.name).latex_name}$" if i.name != 1 else "N/A"
    )

hs.plot(ax=ax, stack=True, histtype="fill")
h.project("trklength").plot(ax=ax, histtype="band")

ax.axvline(10, color="grey", linestyle="dashed", alpha=0.4, zorder=-1)
ax.set_title("NuMI Run 3b")
ax.set_title("Leading Tracks", loc="left")
utils.hist_bin_label(ax, h, "Tracks", "cm")

ax.legend()
utils._save_plot(sel.config, fig, "muon/muon_leading_track_length")
plt.show()

h_all_tracks = (
    hist.Hist.new.Regular(100, 0, 700, name="trklength", label="Length [cm]", flow=True)
    .Regular(50, -1, 1, name="trkllrpid", label="Track LLRPID", flow=True)
    .Regular(50, 0, 5, name="nusep", label=r"$\nu$ separation [cm]", flow=True)
    .StrCategory([x.name for x in SampleType], name="sampletype", label="Sample")
    .IntCategory([1] + list(PDG), name="truepdg", label="True PDG")
    .Bool(name="signal", label="Signal")
    .Weight()
)

barr = ak.broadcast_arrays(
    arr["trk_length"][trk_cond],
    arr["trk_llrpid"][trk_cond],
    ak.drop_none(
        utils.displacement(arr, "trk_start_x", "trk_start_y", "trk_start_z")[trk_cond]
    ),
    arr["sample"],
    np.abs(arr["pfp_true_pdg"])[trk_cond],
    signal_def(arr),
    arr["weight"],
)

h_all_tracks.fill(
    trklength=ak.flatten(barr[0]),
    trkllrpid=ak.flatten(barr[1]),
    nusep=ak.flatten(barr[2]),
    sampletype=ak.flatten(barr[3]),
    truepdg=ak.flatten(barr[4]),
    signal=ak.flatten(barr[5]),
    weight=ak.flatten(barr[6]),
)

fig, ax = plt.subplots(figsize=(10, 6), layout="constrained")

hsa = h_all_tracks.project("trkllrpid", "truepdg").stack("truepdg")
hsa = hist.Stack.from_iter(filter(lambda s: not s.empty(flow=True), hsa))

for i in hsa:
    i.name = (
        f"${part.Particle.from_pdgid(i.name).latex_name}$" if i.name != 1 else "N/A"
    )

hsa.plot(ax=ax, yerr=False, density=True)

ax.set_title("NuMI Run 3b")
utils.hist_bin_label(ax, h_all_tracks.project("trkllrpid"), "Arb.")
ax.legend()

utils._save_plot(sel.config, fig, "muon/numu_cc_track_llrpid")
plt.show()


fig, ax = plt.subplots(layout="constrained")

hsa = h_all_tracks.project("nusep", "truepdg").stack("truepdg")
hsa = hist.Stack.from_iter(filter(lambda s: not s.empty(flow=True), hsa))

for i in hsa:
    i.name = (
        f"${part.Particle.from_pdgid(i.name).latex_name}$" if i.name != 1 else "N/A"
    )

hsa.plot(ax=ax, stack=True, histtype="fill")
h_all_tracks.project("nusep").plot(ax=ax, histtype="band")

ax.set_title("NuMI Run 3b")
utils.hist_bin_label(ax, h_all_tracks.project("nusep"), "Tracks", "cm")
ax.legend()

utils._save_plot(sel.config, fig, "muon/numu_cc_track_separation")
plt.show()


# def track_plots():
#     fig, ax = plt.subplots()
#     h_signal.project("ntracks", "sampletype").stack("sampletype").plot(
#         ax=ax, histtype="fill", stack=True
#     )
#     h_signal.project("ntracks").plot(ax=ax, histtype="band")

#     ax.set_title("NuMI Run 3b (Signal)")
#     ax.set_ylabel("Events")

#     ax.legend()
#     utils._save_plot(sel.config, fig, "topology_tracks_signal")
#     plt.show()

#     fig, ax = plt.subplots()

#     h_background.project("ntracks", "sampletype").stack("sampletype").plot(
#         ax=ax, histtype="fill", stack=True
#     )
#     h_background.project("ntracks").plot(ax=ax, histtype="band")

#     ax.set_title("NuMI Run 3b (Background)")
#     ax.set_ylabel("Events")

#     ax.legend()
#     utils._save_plot(sel.config, fig, "topology_tracks_background")
#     plt.show()


# if __name__ == "__main__":
#     track_plots()
