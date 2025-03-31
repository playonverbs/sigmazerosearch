import logging
from pathlib import Path

import awkward as ak
import hist
import matplotlib.pyplot as plt

import sigmazerosearch.utils as utils
from sigmazerosearch.alg.muon import select_mu_candidate
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
        Cut(
            "showers",
            lambda arr: utils.npfp(arr, "shower") >= 1,
        ),
        Cut(
            "tracks",
            lambda arr: utils.npfp(arr, "track") >= 3,
        ),
        Cut(
            "muon-id",
            lambda arr: ak.sum(select_mu_candidate(arr, pset), axis=1) >= 1,
        ),
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

    plt.style.use("plots/sigmazerosearch.tex.mplstyle")
    plt.style

    signal_cond = (array.pfp_trk_shr_score > 0.5) & (
        array.pfp_true_origin == OriginType.SigmaZero.value
    )

    background_cond = (array.pfp_trk_shr_score > 0.5) & (
        array.pfp_true_origin != OriginType.SigmaZero.value
    )

    n_signal = ak.sum(signal_cond)
    n_background = ak.sum(background_cond)

    print(signal_cond)

    print(f"final cut: {sel.cuts[-1].name}")
    print(f"num signal pfps: {n_signal}")
    print(f"num background pfps: {n_background}")

    h = (
        hist.Hist.new.Regular(
            50, -0.4, 0.2, name="photonbdt", label="Photon BDT Response"
        )
        .Regular(10, 0, 1, name="trkshrscore", label="Track/Shower Score")
        .IntCategory(list(PDG), name="truepdg", label="Backtracked PDG")
        .IntCategory(list(OriginType), name="trueorigin", label="Backtracked Origin")
        .Bool(name="signalevt", label="Signal Event")
        .Bool(name="signalpfp", label="Signal PFP")
        .Double()
    )

    h.fill_flattened(
        photonbdt=array.pfp_photon_bdt_score,
        trkshrscore=array.pfp_trk_shr_score,
        truepdg=array.pfp_true_pdg,
        trueorigin=array.pfp_true_origin,
        signalevt=signal_def(array),
        signalpfp=signal_cond,
    )

    fig, ax = plt.subplots()

    st = (
        h[{"trkshrscore": slice(0j, 0.5j)}]
        .project("photonbdt", "truepdg")
        .stack("truepdg")
    )
    st = hist.Stack.from_iter(filter(lambda s: s.sum() != 0.0, st))

    for hi in st:
        hi.name = PDG(hi.name).name

    st.plot(ax=ax, histtype="fill", stack=True)

    ax.set_title("Photon BDT Response", loc="left")
    ax.set_title("NuMI Run 3b")
    ax.set_ylabel(rf"\# PFPs / {st.axes.widths[0][0]:.2f}")
    ax.axvline(-0.0179, color="gray")

    ax.legend(title=f"Showers after {sel.cuts[-1].name}", ncols=2, fontsize="small")
    utils._save_plot(sel.config, fig, "all_events_pfp_photon_bdt_score")
    plt.show()

    fig, ax = plt.subplots()

    signal_scale = 100
    st2 = (
        h[{"trkshrscore": slice(0j, 0.5j)}]
        .project("photonbdt", "signalpfp")
        .stack("signalpfp")
    )

    st2[True] *= signal_scale

    st2[True].name = rf"Signal (BDT) $\times {signal_scale}$"
    st2[False].name = r"Background (BDT)"

    st2.plot(ax=ax, histtype="fill", stack=True)
    ax.set_title("Photon BDT Response", loc="left")
    ax.set_title("NuMI Run 3b")
    ax.set_ylabel(rf"\# PFPs / {st.axes.widths[0][0]:.2f}")
    ax.axvline(-0.0179, color="gray")
    ax.legend(title=f"Showers after {sel.cuts[-1].name}", fontsize="small")
    utils._save_plot(sel.config, fig, "all_events_by_bdt_cat_pfp_photon_bdt_score")
    plt.show()


if __name__ == "__main__":
    main()
