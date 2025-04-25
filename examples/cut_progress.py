#!/usr/bin/env python3

import awkward as ak
import matplotlib.pyplot as plt

from sigmazerosearch import defaults
from sigmazerosearch.alg.muon import select_mu_candidate
from sigmazerosearch.selection import (
    Cut,
    ParameterSet,
    Sample,
    SampleSet,
    SampleType,
    Selection,
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

sel = Selection(
    config=defaults.config,
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
        # Cut(
        #     "ct-test",
        #     lambda arr: island.count_event_islands(arr, pset).ct_test_islands_plane0 <= 10,
        # )
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

sel.apply_cut(sel.cuts)

sel.cut_summary(header=True, format="md")

# arr = sel.apply_cut(sel.cuts, accumulate=True)

# sel.cut_summary(header=True, format="md")

# h = (
#     hist.Hist.new.Integer(0, 50, name="ct_islands_plane0", flow=True)
#     # .Integer(0, 100, name="ct_islands_plane1")
#     # .Integer(0, 100, name="ct_islands_plane2")
#     .Bool(name="signal")
#     .Double()
# )

# isls = island.count_event_islands(arr, pset)

# h.fill_flattened(
#     ct_islands_plane0=isls.ct_test_islands_plane0,
#     signal=signal_def(arr)
# )

# fig, ax = plt.subplots()

# h.stack("signal").plot(ax=ax)

# plt.show()

# with open("ct_island_count_noshowers.pkl", "wb") as f:
#     pickle.dump(h, f)

# sel.plot_eff_pur()

# h = sel.hist_at(sel.cuts[-1:], dim="event")

# h.project("Tracks", "Showers").plot2d()
