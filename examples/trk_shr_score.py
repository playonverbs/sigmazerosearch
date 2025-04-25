import logging
from pathlib import Path

import hist
import matplotlib.pyplot as plt

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
    config=Config(
        iterate=True,
        iterate_step="50 MB",
        plot_save=True,
        plot_dir=Path("plots"),
        plot_format="svg",
    ),
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

    h = (
        hist.Hist.new.IntCategory(
            list(PDG), name="truepdg", label="True Backtracked PDG"
        )
        .Regular(10, 0, 1, name="cmp", label="PFP Completeness")
        .Regular(10, 0, 1, name="pur", label="PFP Purity")
        .Regular(50, -1, 1, name="llrpid", label="Track LLRPID")
        .Regular(50, 0, 1, name="trkshrscore", label="Track/Shower Score")
        .StrCategory([cut.name for cut in sel.cuts], name="cut", label="Cut Number")
        .IntCategory(list(OriginType), name="trueorigin", label="True Origin")
        .Boolean(name="signalevt", label="Signal Event")
        .Double()
    )

    # h.fill_flattened(
    #     truepdg=array[signal_def(array)].pfp_true_pdg,
    #     trueorigin=array[signal_def(array)].pfp_true_origin,
    #     cut=array[signal_def(array)].cut,
    #     cmp=array[signal_def(array)].pfp_completeness,
    #     pur=array[signal_def(array)].pfp_purity,
    #     llrpid=array[signal_def(array)].trk_llrpid,
    # )

    h.fill_flattened(
        truepdg=array.pfp_true_pdg,
        trueorigin=array.pfp_true_origin,
        cut=array.cut,
        cmp=array.pfp_completeness,
        pur=array.pfp_purity,
        llrpid=array.trk_llrpid,
        trkshrscore=array.pfp_trk_shr_score,
        signalevt=signal_def(array),
    )

    plt.style.use("plots/sigmazerosearch.tex.mplstyle")

    scale = sel.samples.target_POT / sel.samples[0].POT

    _ = len(array[signal_def(array)])

    h *= scale

    fig, ax = plt.subplots()
    st = (
        (h[{"cut": "fv", "signalevt": False}] / scale)
        .project("trkshrscore", "truepdg")
        .stack("truepdg")
    )

    st = hist.Stack.from_iter(filter(lambda s: s.sum() != 0.0, st))

    for hi in st:
        hi.name = PDG(hi.name).name

    ax.set_yscale("log")
    st.plot(ax=ax, stack=True, histtype="fill")

    ax.set_title("Track/Shower Scores in\nBackground Events", loc="left")
    ax.set_title("NuMI Run 3b")
    ax.set_ylabel(rf"\# PFPs / {st.axes.widths[0][0]}")
    ax.legend(ncols=3, fontsize="small")
    utils._save_plot(sel.config, fig, "background_events_pfp_trk_shr_score")

    # ax.set_xticklabels(
    #     [PDG(int(tick.get_text())).name for tick in ax.get_xticklabels()]
    # )
    # ax.axhline(
    #     n_signal, color="lightgrey", linestyle="dashed"
    # )
    # ax.set_title("NuMI Run 3b")
    # ax.set_title("Signal PFPs (Back-tracked)", loc="left")
    # ax.set_ylabel(fr"\# PFPs (raw)")
    # utils._save_plot(sel.config, fig, "signal_pfp_reco_ids")
    plt.show()


if __name__ == "__main__":
    main()
