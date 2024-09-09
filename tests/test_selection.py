import awkward as ak
import numpy as np
import pytest
import uproot as up

from sigmazerosearch.selection import Cut, EventCategory, Sample, SampleType, Selection


@pytest.fixture(
    params=[
        (0, 100),
        (100, 125),
    ]
)
def event_sample(request):
    entry_low, entry_high = request.param
    with up.open(
        "/home/niam/phd/data/hyperons/multiSlice/run3b_RHC/analysisOutputRHC_mSlice_ctslice_cthorpe_make_hyperon_events_numi_rhc_run3b_hyperon_reco2_reco2.root"
    ) as fd:
        return fd.get("ana/OutputTree").arrays(
            entry_start=entry_low, entry_stop=entry_high
        )


@pytest.fixture
def rand_bool_arr():
    return np.random.random_integers(0, 1, 100)


def test_EventCategory(event_sample):
    # get range of EventCategory enum for check.
    min, max = np.amin(list(EventCategory)), np.amax(list(EventCategory))
    res = EventCategory.from_arr(event_sample)

    assert not ak.any(res < min), "an array value is snaller than the enum minimum"
    assert not ak.any(res > max), "an array value is larger than the enum maximum"


@pytest.fixture
def eg_Cut():
    return Cut("bar", lambda: "hello")


def test_Cut(eg_Cut):
    c1 = Cut("this", lambda: 2)
    assert c1.name == "this"
    assert c1() == 2

    c2 = Cut("that", lambda: 3)
    assert c1.name != c2.name
    assert c2() == 3

    assert eg_Cut() == "hello"


def test_Sample():
    s1 = Sample(
        "hyperon",
        "analysisOutputRHC_Hyperon_numi_rhc_run3b.root",
        SampleType.Hyperon,
        1e20,
    )
    s2 = Sample(
        "dirt", "analysisOutputRHC_Dirt_numi_rhc_run3b.root", SampleType.Dirt, 2e30
    )

    assert s1.name != s2.name


def test_Selection():
    s = Selection(
        params={}, samples={}, cuts=[Cut("cut1", lambda: 2), Cut("cut2", lambda: 3)]
    )

    s.cuts[0].n_signal = (60, 10)
    s.cuts[0].n_passing = (40, 5)
    s.cuts[0].n_background = (1e6, 3e4)
    s.cuts[1].n_signal = (30, 10)
    s.cuts[1].n_passing = (15, 5)
    s.cuts[1].n_background = (1e4, 1e2)
