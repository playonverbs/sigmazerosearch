"""
Provides default definitions for parameters and configuration options.

Currently the only additions here are the equivalent frozen 'tunes' created in the
Lambda analysis. These are prepended with `cthorpe_` to label as them being
from that analysis.
"""

from pathlib import Path

from sigmazerosearch.general import Config
from sigmazerosearch.selection import ParameterSet

cthorpe_FHC_Tune_325 = ParameterSet(
    min_length=10,
    max_separation=1,
    separation_cut=3,
    pid_cut=0.6,
    proton_pid_cut=0.1,
    pion_pid_cut=-0.1,
    w_lambda_min=1,
    w_lambda_max=2,
)
"""
Corresponds to [Frozen FHC Tune 325](https://github.com/cthorpe123/HyperonAnalysis/blob/master/ParameterSets/Parameters.h#L29)
"""

cthorpe_RHC_Tune_397 = ParameterSet(
    min_length=10,
    max_separation=1,
    separation_cut=3,
    pid_cut=0.6,
    proton_pid_cut=0.1,
    pion_pid_cut=-0.1,
    w_lambda_min=1,
    w_lambda_max=2,
)
"""
Corresponds to [Frozen RHC Tune 397](https://github.com/cthorpe123/HyperonAnalysis/blob/master/ParameterSets/Parameters.h#L53)
"""

config = Config(
    plot_save=True,
    plot_dir=Path("plots/"),
    data_dir=Path("data/"),
    plot_format=["png", "svg"],
    perf_percent=True,
    iterate=True,
    iterate_step="40 MB",
    branch_list=[
        "run",
        "subrun",
        "event",
        "flash_match_nu_slice_ID",
        "reco_primary_vtx_inFV",
        "reco_primary_vtx_x",
        "reco_primary_vtx_y",
        "reco_primary_vtx_z",
        "pfp_trk_shr_score",
        "mc_nu_pdg",
        "mc_hyperon_pdg",
        "mc_nu_pos_x",
        "mc_nu_pos_y",
        "mc_nu_pos_z",
        "mc_decay_pdg",
        "trk_llrpid",
        "trk_length",
        "trk_start_x",
        "trk_start_y",
        "trk_start_z",
        "pfp_photon_bdt_score",
        "ct_test_window_plane0",
        "ct_test_window_plane1",
        "ct_test_window_plane2",
        "ct_test_primary_vtx_wires",
    ],
)

config_no_window = config.with_(
    {
        "branch_list": [
            "run",
            "subrun",
            "event",
            "flash_match_nu_slice_ID",
            "reco_primary_vtx_inFV",
            "reco_primary_vtx_x",
            "reco_primary_vtx_y",
            "reco_primary_vtx_z",
            "pfp_trk_shr_score",
            "mc_nu_pdg",
            "mc_hyperon_pdg",
            "mc_nu_pos_x",
            "mc_nu_pos_y",
            "mc_nu_pos_z",
            "mc_decay_pdg",
            "trk_llrpid",
            "trk_length",
            "trk_start_x",
            "trk_start_y",
            "trk_start_z",
            "pfp_photon_bdt_score",
        ]
    }
)

config_no_bdt = config.with_(
    {
        "branch_list": [
            "run",
            "subrun",
            "event",
            "flash_match_nu_slice_ID",
            "reco_primary_vtx_inFV",
            "reco_primary_vtx_x",
            "reco_primary_vtx_y",
            "reco_primary_vtx_z",
            "pfp_trk_shr_score",
            "mc_nu_pdg",
            "mc_hyperon_pdg",
            "mc_nu_pos_x",
            "mc_nu_pos_y",
            "mc_nu_pos_z",
            "mc_decay_pdg",
            "trk_llrpid",
            "trk_length",
            "trk_start_x",
            "trk_start_y",
            "trk_start_z",
            "ct_test_window_plane0",
            "ct_test_window_plane1",
            "ct_test_window_plane2",
            "ct_test_primary_vtx_wires",
        ]
    }
)
