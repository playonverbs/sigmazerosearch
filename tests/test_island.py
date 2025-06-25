import matplotlib.pyplot as plt
import pytest
import uproot as up

from sigmazerosearch.alg import island


@pytest.fixture(
    scope="session",
    params=[
        (0, 100),
        (100, 125),
    ],
)
def event_sample(request):
    entry_low, entry_high = request.param
    with up.open(
        "/home/niam/phd/today/2025-05-29/analysisOutputRHC_mSlice_newflux_Hyperon_cthorpe_make_hyperon_events_numi_run3b_hyperon_reco2_reco2.root"
    ) as fd:
        yield fd.get("ana/OutputTree").arrays(
            entry_start=entry_low, entry_stop=entry_high
        )


def test_window_to_map_single(event_sample, benchmark):
    windows = event_sample.ct_test_window_plane0
    benchmark(island._window_to_map, windows[0], 250, 100)


def test_window_to_map_numpy(event_sample, benchmark):
    windows = event_sample.ct_test_window_plane0
    map = benchmark(island._window_to_map_numpy, windows[0], 250, 100)

    # assert np.array_equal(map, island._window_to_map(windows[0], 250, 100))
    plt.matshow(map)
    plt.gca().set_title("_window_to_map")
    plt.matshow(island._window_to_map(windows[0], 250, 100))
    plt.gca().set_title("_window_to_map_numpy")
    plt.show()
