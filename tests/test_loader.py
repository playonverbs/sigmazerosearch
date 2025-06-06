import matplotlib.pyplot as plt
import numpy as np
import uproot as up

from sigmazerosearch import loader


def test_loader():
    ntuple = loader.load_ntuple(
        "/home/niam/phd/data/hyperons/run3b_RHC/analysisOutputRHC_Overlay_GENIE_Hyperon_cthorpe_make_hyperon_events_numi_rhc_run3b_hyperon_reco2_reco2.root:ana/OutputTree"
    )

    assert isinstance(ntuple, up.TTree)


def test_load_dead_wire_regions():
    linear = loader.load_dead_wire_map("data/MCC9_channel_list.txt", transpose=False)
    transposed = loader.load_dead_wire_map("data/MCC9_channel_list.txt", transpose=True)

    assert np.min(linear) >= 0 and np.max(linear) < 8256

    fig, ax = plt.subplots(layout="constrained")

    [ax.axvline(x=i, c="lightgrey", alpha=0.3) for i in transposed[2]]
    [ax.axvline(x=i, c="lightblue", alpha=0.3) for i in transposed[1]]

    plt.show()
