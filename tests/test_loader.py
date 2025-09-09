import matplotlib.pyplot as plt
import numpy as np
import uproot as up

from sigmazerosearch import loader


def test_loader():
    ntuple = loader.load_ntuple(
        "/home/niam/phd/data/hyperons/multiSlice/run3b_RHC/analysisOutputRHC_cttest_Hyperon_ana_hyperon_sigmazero_mSlice_prodgenie_run3b.root:ana/OutputTree"
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
