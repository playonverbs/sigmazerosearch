import awkward as ak
import pytest
import uproot as up

from sigmazerosearch import defaults, mva, selection, truth, utils


@pytest.fixture
def Selection_array():
    pset = selection.ParameterSet(
        pid_cut=0.6,
        min_length=10,
        max_separation=1,
        proton_pid_cut=0.2,
        pion_pid_cut=0.2,
        separation_cut=3,
        w_lambda_min=1.1,
        w_lambda_max=1.20,
    )

    sel = selection.Selection(
        config=selection.Config(
            iterate=True,
            iterate_step="50MB",
            branch_list=(
                defaults.config_no_bdt.branch_list + ["pfp_x", "pfp_y", "pfp_z"]  # type: ignore
            ),
        ),
        cuts=[
            selection.Cut(
                "has-flash-match",
                lambda arr: arr["flash_match_nu_slice_ID"] != -1,
            ),
            selection.Cut("fv", lambda arr: arr["reco_primary_vtx_inFV"]),
        ],
        samples=selection.SampleSet(
            selection.Sample(
                "hyperon",
                "/home/niam/phd/data/hyperons/multiSlice/run3b_RHC/analysisOutputRHC_cttest_Hyperon_ana_hyperon_sigmazero_mSlice_prodgenie_run3b.root",
                selection.SampleType.Hyperon,
                None,
            ),
            target_POT=4.9e20,
        ),
        params=pset,
        label="bdt_prepare",
    )
    sel.open_files()

    return sel.apply_cut(sel.cuts, accumulate=True)


@pytest.fixture
def make_MVAHandler(Selection_array):
    return mva.MVAHandler(
        data=Selection_array,
        target_branches={
            "shr_open_angle",
            "shr_length",
            "pfp_trk_shr_score",
            "trk_llrpid",
            "pfp_true_pdg",
            "pfp_true_origin",
            "nusep",  # calculated in this script
        },
        extra_fields={
            "nusep": lambda arr: utils.displacement(arr, "pfp_x", "pfp_y", "pfp_z")
        },
        precondition=mva.precond_has_pfps,
        signal_condition=(
            lambda output: output.pfp_true_origin == truth.OriginType.SigmaZero.value
        ),
        background_condition=(
            lambda output: output.pfp_true_origin != truth.OriginType.SigmaZero.value
        ),
    )


def test_MVAHandler(make_MVAHandler, tmp_path):
    print(make_MVAHandler)

    for field in make_MVAHandler.extra_fields.keys():
        assert field in ak.fields(make_MVAHandler.data), f"{field} not found in array"

    output_path = tmp_path / "test.root"

    make_MVAHandler.save_bdt_trees(output_path)

    fd = up.open(output_path)

    signal_key, background_key = (
        f"{make_MVAHandler.mva_method}/SignalTree",
        f"{make_MVAHandler.mva_method}/SignalTree",
    )

    _ = fd.get(signal_key)
    _ = fd.get(background_key)
