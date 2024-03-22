#!/usr/bin/env python3
from sigmazerosearch.types.selection import (
    Selection,
    Cut,
    Sample,
    SampleSet,
    SampleType,
)
import awkward as ak

sel = Selection(
    cuts=[
        Cut(
            "fv",
            "reco_primary_vtx == True",
            lambda arr: arr["reco_primary_vtx_inFV"] == True,
        ),
        Cut(
            "tracks",
            "NPrimaryTracklikeDaughters >= 3",
            lambda arr: ak.count(arr["pfp_trk_shr_score"] > 0.5, axis=1) >= 3,
        ),
        Cut(
            "showers",
            "NPrimaryShowerlikeDaughters >= 1",
            lambda arr: ak.count(arr["pfp_trk_shr_score"] < 0.5, axis=1) >= 1,
        ),
        # Cut("muon-id", "..."),
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
        base_dir="/home/niam/phd/data/hyperons/run3b_RHC",
    ),
    params=[2.0],
)


def main():
    # TODO: wrap Selection with a context manager `with open()`-like
    # needs: sel.__enter__(...) and sel.__exit(...)
    print("Opening files")
    sel.open_files()

    for sam in sel.samples:
        print(sam.file_name.split("/")[-1])
    for n, cut in enumerate(sel.cuts):
        print(f"Applying {n}. {cut.name} Cut")
        sel.apply_cut(cut.name)
        print("...Applied cut")

    sel.cut_summary(header=True)

    # sel.plot_eff_pur()
    print("Closing files")
    sel.close_files()


if __name__ == "__main__":
    main()
