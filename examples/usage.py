#!/usr/bin/env python3
from sigmazerosearch.types.selection import (
    Selection,
    Cut,
    Sample,
    SampleSet,
    SampleType,
)

sel = Selection(
    cuts=[
        Cut("fv", "InActiveTPC == True"),
        Cut("tracks", "NPrimaryTracklikeDaughters >= 3"),
        Cut("showers", "NPrimaryShowerlikeDaughters >= 1"),
        Cut("muon-id", "..."),
        Cut("lambda-analysis", "LambdaBDT > 0.4"),
        Cut("photon-bdt", "PhotonBDT > 0.2"),
    ],
    samples=SampleSet(
        Sample(
            "hyperon",
            "/home/niam/phd/data/hyperons/run3b_RHC/analysisOutputRHC_Overlay_GENIE_Hyperon_cthorpe_make_hyperon_events_numi_rhc_run3b_hyperon_reco2_reco2.root",
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
        print(sam.file_name)
    for n, cut in enumerate(sel.cuts[0:1]):
        print(f"Applying {n} {cut.name} Cut")
        sel.apply_cut(cut.name)
        print("...Applied cut")

    # sel.plot_eff_pur()
    print("Closing files")
    sel.close_files()


if __name__ == "__main__":
    main()
