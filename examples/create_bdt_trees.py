import logging

import awkward as ak
import numpy as np
import uproot as up

import sigmazerosearch.utils as utils
from sigmazerosearch.general import PDG, Config
from sigmazerosearch.selection import (
    Cut,
    ParameterSet,
    Sample,
    SampleSet,
    SampleType,
    Selection,
    signal_def,
)
from sigmazerosearch.truth import OriginType

logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)

pset = ParameterSet(
    pid_cut=0.6,
    min_length=10,
    max_separation=1,
    proton_pid_cut=0.2,
    pion_pid_cut=0.2,
    separation_cut=3,
    w_lambda_min=1.1,
    w_lambda_max=1.20,
)

sel = Selection(
    config=Config(iterate=True, iterate_step="50 MB"),
    cuts=[
        Cut(
            "has-flash-match",
            lambda arr: arr["flash_match_nu_slice_ID"] != -1,
        ),
        Cut("fv", lambda arr: arr["reco_primary_vtx_inFV"]),
    ],
    samples=SampleSet(
        Sample(
            "hyperon",
            "/home/niam/phd/data/hyperons/multiSlice/run3b_RHC/analysisOutputRHC_mSlice_cthorpe_make_hyperon_events_numi_rhc_run3b_hyperon_reco2_reco2.root",
            SampleType.Hyperon,
            None,
        ),
        Sample(
            "background_all",
            "/home/niam/phd/data/hyperons/multiSlice/run3b_RHC/analysisOutputRHC_mSlice_Background_all.root",
            SampleType.Background,
            None,
        ),
        target_POT=4.9e20,
    ),
    params=pset,
    label="bdt_prepare",
)


def photon_signal(arr: ak.Array) -> ak.Array:
    return np.logical_and.reduce(
        [
            signal_def(arr),
            ak.sum(arr.pfp_true_pdg == PDG.Photon.value, axis=1) == 1,
            # arr.pfp_true_pdg == PDG.Photon
        ]
    )  # type: ignore


# define only the branches we need for signal/background differentiation and
# the bdt training
BDT_BRANCHES = {
    "shr_open_angle",
    "shr_length",
    "pfp_trk_shr_score",
    "trk_llrpid",
    # "trk_three_plane_dedx",
    "pfp_true_pdg",
    "pfp_true_origin",
    "nusep",  # calculated in this script
}

OUTPUT_PATH = "~/phd/data/hyperons/multiSlice/BDT/photonBDT/photonbdt_train_run3b_hyperon_background.root"


def main():
    sel.samples.as_table()
    sel.open_files()

    logging.info(
        f"Applying the following cuts consecutively {', '.join([c.name for c in sel.cuts])}"
    )
    array = sel.apply_cut(sel.cuts, accumulate=True)

    if array is None:
        print("Could not get selected array")
        return

    sel.cut_summary(header=True, format="md")

    bdt_dir = up.recreate(OUTPUT_PATH)

    # Write signal and background trees via applying the bdt signal and it's negation.
    # Trees intended for TMVA usage must be flat trees

    array["nusep"] = utils.displacement(array, "pfp_x", "pfp_y", "pfp_z")

    # cond = np.logical_and(
    #     array['trk_three_plane_dedx'] >= 0,
    #     array['trk_three_plane_dedx'] <= 999,
    # )

    d = dict(
        zip(
            ak.fields(array[BDT_BRANCHES]),
            ak.unzip(
                array[BDT_BRANCHES][
                    ak.num(
                        # array[BDT_BRANCHES][array.pfp_trk_shr_score <= 0.5].pfp_true_pdg
                        array[BDT_BRANCHES].pfp_true_pdg
                    )
                    > 0
                ]
            ),
        )
    )

    for k, v in d.items():
        d[k] = ak.flatten(v)

    output = ak.zip(d)

    logging.info(f"writing trees to {OUTPUT_PATH}")

    bdt_dir["bdt/SignalTree"] = output[
        output.pfp_true_origin == OriginType.SigmaZero.value
    ]  # [BDT_BRANCHES - {"pfp_true_pdg", "pfp_true_origin"}]

    bdt_dir["bdt/BackgroundTree"] = output[
        output.pfp_true_origin != OriginType.SigmaZero.value
    ]  # [BDT_BRANCHES - {"pfp_true_pdg", "pfp_true_origin"}]

    logging.info(
        f"Saved {bdt_dir['bdt/SignalTree'].num_entries} Signal entries and {bdt_dir['bdt/BackgroundTree'].num_entries} Background entries"
    )

    sel.close_files()


if __name__ == "__main__":
    main()
