from ..sigmazerosearch.types.selection import (
    Selection,
    Cut,
    Sample,
    SampleSet,
    SampleType,
)

sel = Selection(
    cuts=[
        Cut("fv", "inFV"),
        Cut("tracks", "NPrimaryTracklikeDaughters >= 3"),
        Cut("showers", "NPrimaryShowerlikeDaughters >= 1"),
        Cut("muon-id", "..."),
        Cut("lambda-analysis", "LambdaBDT > 0.4"),
        Cut("photon-bdt", "PhotonBDT > 0.2"),
    ],
    samples=SampleSet(
        [
            Sample(
                "hyperon",
                "analysisOutputRHC_Overlay_GENIE_Hyperon_cthorpe_make_hyperon_events_numi_rhc_run3b_hyperon_reco2_reco2.root",
                SampleType.Hyperon,
                None,  # TODO: fix this to allow getting POT from the file MetaTree
            ),
        ]
    ),
)

for cut in sel.cuts:
    print(f"Applying {cut.name} Cut")
    sel.apply_cut(cut.name)

sel.plot_eff_pur()

sel.close_files()
