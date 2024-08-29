from sigmazerosearch.loader import load_ntuple


def test_loader():
    ntuple = load_ntuple(
        "/home/niam/phd/data/hyperons/run3b_RHC/analysisOutputRHC_Overlay_GENIE_Hyperon_cthorpe_make_hyperon_events_numi_rhc_run3b_hyperon_reco2_reco2.root:ana/OutputTree"
    )

    print(ntuple.keys())
