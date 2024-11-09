#!/usr/bin/env python3
import logging

import awkward as ak
import matplotlib.pyplot as plt

from sigmazerosearch.alg.muon import select_mu_candidate
from sigmazerosearch.selection import (
    Cut,
    ParameterSet,
    Sample,
    SampleSet,
    SampleType,
    Selection,
)

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
    cuts=[
        # Cut(
        #     "has-reco-nu",
        #     lambda arr: arr["reco_primary_vtx_x"] > -999,
        # ),
        Cut(
            "fv",
            lambda arr: arr["reco_primary_vtx_inFV"],
        ),
        Cut(
            "tracks",
            lambda arr: ak.sum(arr["pfp_trk_shr_score"] >= 0.5, axis=1) >= 3,
        ),
        Cut(
            "showers",
            lambda arr: ak.sum(arr["pfp_trk_shr_score"] < 0.5, axis=1) >= 1,
        ),
        Cut(
            "muon-id",
            lambda arr: ak.sum(select_mu_candidate(arr, pset), axis=1) >= 1,
        ),
        # Cut("lambda-analysis", "LambdaBDT > 0.4"),
        # Cut("photon-bdt", "PhotonBDT > 0.2"),
    ],
    samples=SampleSet(
        Sample(
            "hyperon",
            "/home/niam/phd/data/hyperons/multiSlice/run3b_RHC/analysisOutput_mSlice_cthorpe_make_hyperon_events_numi_rhc_run3b_hyperon_reco2_reco2.root",
            SampleType.Hyperon,
            None,
        ),
        base_dir="/home/niam/phd/data/hyperons/multiSlice/run3b_RHC",
    ),
    params=pset,
    label="nominal",
)


def main():
    # TODO: wrap Selection with a context manager `with open()`-like
    # needs: sel.__enter__(...) and sel.__exit(...)
    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
    print("Opening files")
    sel.open_files()

    for sam in sel.samples:
        print(sam.file_name.split("/")[-1], sam.type.name, sam.POT, sep="\t")

    sel.apply_cut(sel.cuts)

    sel.cut_summary(header=True, format="md")

    # sel.plot_slice_info(signal=False)
    # sel.plot_slice_info(signal=True)
    # sel.plot_reco_effs(signal=True)
    # sel.plot_eff_pur()

    arr = sel.samples[0].df.arrays()

    fig, ax = plt.subplots()

    h2 = ax.hist2d(
        arr["flash_match_nu_slice_ID"].to_numpy(),
        arr["true_nu_slice_ID"].to_numpy(),
    )
    fig.colorbar(h2[3], ax=ax, label="Number of slices")
    plt.show()

    print("Closing files")
    sel.close_files()


if __name__ == "__main__":
    main()
