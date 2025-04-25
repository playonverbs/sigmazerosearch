#!/usr/bin/env python3
from pathlib import Path

import awkward as ak
import hist
import matplotlib.pyplot as plt
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
        "pfp_true_pdg",
        "pfp_true_length",
        "pfp_completeness",
        "pfp_purity",
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
            "/home/niam/phd/data/hyperons/multiSlice/run3b_RHC/cttest/analysisOutputRHC_mSlice_cttest_nuslice_Hyperon_cthorpe_make_hyperon_events_numi_run3b_hyperon_reco2_reco2.root",
            SampleType.Hyperon,
            None,
        ),
        Sample(
            "background_all",
            "/home/niam/phd/data/hyperons/multiSlice/run3b_RHC/cttest/analysisOutputRHC_mSlice_cttest_nuslice_Background_all_prodgenie_numi_uboone_overlay_rhc_mcc9_run3b.root",
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
    hist.Hist.new.Regular(
        10, 0, 1, name="completeness", label="PFP Completeness", flow=True
    )
    .Regular(10, 0, 1, name="purity", label="PFP Purity")
    .Regular(50, 0, 1, name="trkshrscore", label="Track/Shower Score")
    .StrCategory([x.name for x in SampleType], name="sampletype", label="Sample")
    .IntCategory([1] + list(PDG), name="truepdg", label="True PDG")
    .Bool(name="signal", label="Signal")
    .Weight()
)

barr = ak.broadcast_arrays(
    arr.pfp_completeness,
    arr.pfp_purity,
    arr.pfp_trk_shr_score,
    arr.sample,
    # np.abs(arr.pfp_true_pdg),
    arr.pfp_true_pdg,
    signal_def(arr),
    arr.weight,
)

# [b.show() for b in barr]

h.fill(
    completeness=ak.flatten(barr[0]),
    purity=ak.flatten(barr[1]),
    trkshrscore=ak.flatten(barr[2]),
    sampletype=ak.flatten(barr[3]),
    truepdg=ak.flatten(barr[4]),
    signal=ak.flatten(barr[5]),
    weight=ak.flatten(barr[6]),
)

# h_track = h[:, :, hist.loc(0.5):, ...]
# h_shower = h[:, :, :hist.loc(0.5), ...]

# fig, ax = plt.subplots()

# h_track.project("trkshrscore").plot(ax=ax, label="tracks")
# h_shower.project("trkshrscore").plot(ax=ax, label="showers")
# ax.legend()
# plt.show()

fig, ax = plt.subplots(layout="constrained")

hss = h[{"signal": True}].project("trkshrscore", "truepdg").stack("truepdg")
hss = hist.Stack.from_iter(filter(lambda s: not s.empty(), hss))
hsb = h[{"signal": False}].project("trkshrscore", "truepdg").stack("truepdg")
hsb = hist.Stack.from_iter(filter(lambda s: not s.empty(), hsb))

for i in hss:
    i.name = (
        f"${part.Particle.from_pdgid(i.name).latex_name}$" if i.name != 1 else "N/A"
    )

for i in hsb:
    i.name = (
        f"${part.Particle.from_pdgid(i.name).latex_name}$" if i.name != 1 else "N/A"
    )

hss.plot(ax=ax, stack=True, histtype="fill")
ax.set_yscale("log")
ax.set_title("NuMI Run 3b (Signal)")
ax.legend(ncols=4, fontsize="small")
utils._save_plot(sel.config, fig, "pfps/signal_events_pfp_trk_shr_score")
plt.show()

fig, ax = plt.subplots(layout="constrained")

hsb.plot(ax=ax, stack=True, histtype="fill")
ax.set_yscale("log")
ax.set_title("NuMI Run 3b (Background)")
ax.legend(ncols=4, fontsize="small")
utils._save_plot(sel.config, fig, "pfps/background_events_pfp_trk_shr_score")
plt.show()


fig, axs = plt.subplots(ncols=2, figsize=(10, 6), layout="constrained")

hs_shower = h[:, :, 0j:0.5j:sum, sum, :, sum].stack("truepdg").project("completeness")
hs_shower = hist.Stack.from_iter(filter(lambda s: not s.empty(flow=True), hs_shower))

for i in hs_shower:
    i.name = (
        f"${part.Particle.from_pdgid(i.name).latex_name}$" if i.name != 1 else "N/A"
    )

hs_shower.plot(ax=axs[0])

axs[0].set_title("Showers", loc="left")

hs_track = h[:, :, 0.5j:1j:sum, sum, :, sum].stack("truepdg").project("completeness")
hs_track = hist.Stack.from_iter(filter(lambda s: not s.empty(flow=True), hs_track))

for i in hs_track:
    i.name = (
        f"${part.Particle.from_pdgid(i.name).latex_name}$" if i.name != 1 else "N/A"
    )

hs_track.plot(ax=axs[1], label="_")

axs[1].set_title("NuMI Run 3b")
axs[1].set_title("Tracks", loc="left")

fig.legend(loc="outside right center", frameon=False)
utils._save_plot(sel.config, fig, "pfps/pfp_completeness_tracks_showers_fv")
plt.show()

fig, axs = plt.subplots(ncols=2, figsize=(10, 6), layout="constrained")

hs_shower = h[:, :, 0j:0.5j:sum, sum, :, sum].stack("truepdg").project("purity")
hs_shower = hist.Stack.from_iter(filter(lambda s: not s.empty(flow=True), hs_shower))

for i in hs_shower:
    i.name = (
        f"${part.Particle.from_pdgid(i.name).latex_name}$" if i.name != 1 else "N/A"
    )

hs_shower.plot(ax=axs[0])

axs[0].set_title("Showers", loc="left")

hs_track = h[:, :, 0.5j:1j:sum, sum, :, sum].stack("truepdg").project("purity")
hs_track = hist.Stack.from_iter(filter(lambda s: not s.empty(flow=True), hs_track))

for h in hs_track:
    h.name = (
        f"${part.Particle.from_pdgid(h.name).latex_name}$" if h.name != 1 else "N/A"
    )

hs_track.plot(ax=axs[1], label="_")

axs[1].set_title("NuMI Run 3b")
axs[1].set_title("Tracks", loc="left")

fig.legend(loc="outside right center", frameon=False)
utils._save_plot(sel.config, fig, "pfps/pfp_purity_tracks_showers_fv")
plt.show()
