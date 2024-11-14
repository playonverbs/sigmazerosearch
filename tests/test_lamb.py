import pytest
import uproot as up

from sigmazerosearch.alg.lamb import select_p_pi_candidates_box
from sigmazerosearch.selection import ParameterSet


@pytest.fixture(
    params=[
        (0, 100),
        (100, 125),
    ]
)
def event_sample(request):
    entry_low, entry_high = request.param
    with up.open(
        "/home/niam/phd/data/hyperons/multiSlice/run3b_RHC/analysisOutputRHC_mSlice_cthorpe_make_hyperon_events_numi_rhc_run3b_hyperon_reco2_reco2.root"
    ) as fd:
        return fd.get("ana/OutputTree").arrays(
            entry_start=entry_low, entry_stop=entry_high
        )


@pytest.fixture()
def pset():
    return ParameterSet(
        pid_cut=0.6,
        min_length=10,
        max_separation=1,
        proton_pid_cut=0.2,
        pion_pid_cut=0.2,
        separation_cut=3,
        w_lambda_min=1.1,
        w_lambda_max=1.20,
    )


def test_select_p_pi_candidates_box(event_sample, pset):
    what = select_p_pi_candidates_box(event_sample, pset)

    print(what.show())


# def test_reconstruct_p_pi(event_sample, pset):
#     indices = np.random.randint(2, size=(len(event_sample), 2))
#     reconstruct_p_pi(event_sample, pset, ak.Array(indices))
