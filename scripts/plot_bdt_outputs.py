"""
Reads a ROOT file created via running either a standard TMVA workflow or the
`bdt_prep` executable.

Various options are available to plot histograms stored in that file and save
those plots.
"""

import argparse
from dataclasses import dataclass
from typing import Literal

import hist
import matplotlib.pyplot as plt
import numpy as np
import uproot as up


@dataclass
class Variable:
    name: str
    unit: str | None = None

    def __str__(self):
        if self.unit:
            return f"{self.name} [{self.unit}]"
        else:
            return self.name

    def fmt_unit(self) -> str:
        if self.unit:
            return self.unit
        else:
            return "arb"


variable_map = {
    "trk_llrpid": Variable("Track LLRPID"),
    "pfp_trk_shr_score": Variable("Track/Shower Score"),
    "shr_length": Variable("Shower Length", "cm"),
    "shr_open_angle": Variable("Shower Opening Angle", "rad"),
    "Separation from nu vertex": Variable("nu separation", "cm"),
    "trk_three_plane_dedx": Variable("Track Three-Plane dE/dX", "MeV/cm"),
}

# def test_train_regex(tt, sb):
#     re.compile()
#     pass


def get_pairs(fd, filter_name1: str, filter_name2: str, filter_classname):
    return list(
        zip(
            list(
                fd.itervalues(
                    filter_name=filter_name1, filter_classname=filter_classname
                )
            ),
            list(
                fd.itervalues(
                    filter_name=filter_name2, filter_classname=filter_classname
                )
            ),
        )
    )


def get_binwidth(hist):
    return hist.axes[0].widths[0]


def conv_hist(pair):
    return (pair[0].to_hist(), pair[1].to_hist())


def input_hists(fd):
    pairs = get_pairs(
        fd,
        filter_name1="*Input*Signal*",
        filter_name2="*Input*Background*",
        filter_classname="TH1*",
    )

    pairs = list(map(conv_hist, pairs))
    return pairs


def bdt_response_hists(fd):
    # pairs = list(
    #     zip(
    #         # fd.itervalues(
    #         #     filter_name=lambda x: "dd"
    #         # )
    #         fd.itervalues(
    #             filter_name="dataset/Method_*/*/MVA_*_Train_S", filter_classname="TH1*"
    #         ),
    #         fd.itervalues(
    #             filter_name="dataset/Method_*/*/MVA_*_Train_B", filter_classname="TH1*"
    #         ),
    #     )
    # )

    pairs = list(
        (
            fd.get("dataset/Method_BDT/BDT/MVA_BDT_S"),
            fd.get("dataset/Method_BDT/BDT/MVA_BDT_Train_S"),
            fd.get("dataset/Method_BDT/BDT/MVA_BDT_B"),
            fd.get("dataset/Method_BDT/BDT/MVA_BDT_Train_B"),
        )
    )

    pairs = list(map(lambda x: x.to_hist(), pairs))
    return pairs


def mva_performance_hists(fd):
    pairs = list(
        zip(
            fd.itervalues(
                filter_name="dataset/Method_*/*/MVA_*_effS", filter_classname="TH1*"
            ),
            fd.itervalues(
                filter_name="dataset/Method_*/*/MVA_*_effB", filter_classname="TH1*"
            ),
        )
    )
    pairs = list(map(conv_hist, pairs))
    return pairs


def corr_matrices(fd):
    pairs = get_pairs(
        fd,
        filter_name1="dataset/CorrelationMatrixS",
        filter_name2="dataset/CorrelationMatrixB",
        filter_classname="TH2*",
    )

    pairs = list(map(conv_hist, pairs))
    return pairs


def main(args):
    fd = up.open(args.filename)

    try:
        plt.style.use("../plots/sigmazerosearch.tex.mplstyle")
    except Exception:
        pass
    finally:
        plt.style.use("plots/sigmazerosearch.tex.mplstyle")

    if args.inputs:
        inputs = input_hists(fd)

        for sig, bkg in inputs:
            alpha = 0.2
            fig, ax = plt.subplots(layout="constrained")

            bkg.plot(histtype="band", ax=ax, fc="tab:orange", alpha=alpha, hatch=None)
            bkg.plot(ax=ax, label="Background", color="tab:orange", yerr=False)

            sig.plot(histtype="band", ax=ax, fc="tab:blue", alpha=alpha, hatch=None)
            sig.plot(ax=ax, label="Signal", color="tab:blue", yerr=False)

            var = variable_map.get(sig.label, Variable(sig.label))
            ax.set_title("NuMI Run 3b")
            ax.set_title("Photon BDT Inputs", loc="left")
            ax.set_xlabel(str(var))
            ax.set_ylabel(f"PFPs / {get_binwidth(sig):.2f} {var.fmt_unit()}")
            ax.legend()
            if args.save:
                fig.savefig(
                    f"MVA_input_{var.name.replace(' ', '_').replace('/', '_')}.{args.ext}",
                    bbox_inches="tight",
                )

        plt.show()

    if args.response:
        responses = bdt_response_hists(fd)

        sig, sig_train, bkg, bkg_train = responses
        alpha = 0.2
        fig, ax = plt.subplots(layout="constrained")

        # bkg.plot(histtype="band", ax=ax, fc="tab:orange", alpha=alpha, hatch=None)
        bkg.plot(ax=ax, label="Background (Test)", color="tab:orange", yerr=False)
        bkg_train.plot(
            ax=ax,
            label="Background (Train)",
            histtype="errorbar",
            xerr=True,
            c="tab:orange",
            markersize=5,
        )

        # sig.plot(histtype="band", ax=ax, fc="tab:blue", alpha=alpha, hatch=None)
        sig.plot(ax=ax, label="Signal (Test)", color="tab:blue", yerr=False)
        sig_train.plot(
            ax=ax,
            label="Signal (Train)",
            histtype="errorbar",
            xerr=True,
            c="tab:blue",
            markersize=5,
        )

        method = sig.name.split("_")[1]

        ax.set_title("NuMI Run 3b")
        ax.set_title("Photon BDT Response", loc="left")
        ax.set_xlabel("BDT Response")
        ax.legend()
        if args.save:
            fig.savefig(f"MVA_{method}_response.{args.ext}", bbox_inches="tight")

        # for sig, sig_train, bkg, bkg_train in responses:
        #     alpha = 0.2
        #     fig, ax = plt.subplots()

        #     bkg.plot(histtype="band", ax=ax, fc="tab:orange", alpha=alpha, hatch=None)
        #     bkg.plot(ax=ax, label="Background", color="tab:orange", yerr=False)

        #     sig.plot(histtype="band", ax=ax, fc="tab:blue", alpha=alpha, hatch=None)
        #     sig.plot(ax=ax, label="Signal", color="tab:blue", yerr=False)

        #     method = sig.name.split("_")[1]

        #     ax.set_title("NuMI Run 3b")
        #     ax.set_title("Photon BDT Response", loc="left")
        #     ax.set_xlabel("BDT Response")
        #     ax.legend(title=method)
        #     if args.save:
        #         fig.savefig(f"MVA_{method}_response.{args.ext}", bbox_inches="tight")

        plt.show()

    if args.corr:
        sig_corr, bkg_corr = corr_matrices(fd)[0]

        fig, ax = plt.subplots()
        sig_corr.plot2d(ax=ax, labels=True)
        ax.set_title("Photon BDT", loc="left")
        ax.set_title("Correlation (Signal)")
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_xticklabels(
            list(map(lambda x: variable_map[x.get_text()].name, ax.get_xticklabels()))
        )
        ax.set_yticklabels(
            list(map(lambda x: variable_map[x.get_text()].name, ax.get_yticklabels()))
        )
        # ax.tick_params(axis="x", labelrotation=30)
        plt.xticks(rotation=30, ha="right")
        if args.save:
            fig.savefig(f"MVA_input_corr_matrix_signal.{args.ext}", bbox_inches="tight")

        fig, ax = plt.subplots()
        bkg_corr.plot2d(ax=ax, labels=True)
        ax.set_title("Photon BDT", loc="left")
        ax.set_title("Correlation (Background)")
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_xticklabels(
            list(map(lambda x: variable_map[x.get_text()].name, ax.get_xticklabels()))
        )
        ax.set_yticklabels(
            list(map(lambda x: variable_map[x.get_text()].name, ax.get_yticklabels()))
        )
        # ax.tick_params(axis="x", labelrotation=30)
        plt.xticks(rotation=30, ha="right")
        if args.save:
            fig.savefig(
                f"MVA_input_corr_matrix_background.{args.ext}", bbox_inches="tight"
            )

        plt.show()

    if args.roc:
        plot_roc(fd, args=args)

    if args.perf:
        n_signal, n_bkg = (265, 171389)  # XXX: Hardcoded values
        for sig, bkg in mva_performance_hists(fd):
            fig, ax = plt.subplots(layout="constrained")
            sig.plot(ax=ax, label="Signal", yerr=False, w2method="sqrt")
            bkg.plot(
                ax=ax, label="Background", yerr=False, w2method="sqrt"
            )  # yerr=False)
            method = sig.name.split("_")[1]

            ax.set_title("NuMI Run 3b")
            ax.set_title("Photon BDT Efficiencies", loc="left")
            ax.set_xlabel("Response Score Cut")
            ax.axhline(1.0, color="lightgrey", linestyle="dashed")
            h1, l1 = ax.get_legend_handles_labels()
            ax2 = ax.twinx()
            sign, optimum = calc_significance(sig, bkg, n_signal=n_signal, n_bkg=n_bkg)
            sign.plot(
                ax=ax2,
                label=f"Significance (S={n_signal}, B={n_bkg})",
                yerr=False,
                ec="tab:red",
            )
            ax.axvline(optimum, color="lightgrey", linestyle="dashed")
            h2, l2 = ax2.get_legend_handles_labels()
            ax.legend(h1 + h2, l1 + l2, title=method)

            if args.save:
                fig.savefig(f"MVA_{method}_performance.{args.ext}", bbox_inches="tight")

        plt.show()


def plot_roc(
    ntuple,
    args,
    method: Literal[
        "BDT",
        "BDTG",
    ] = "BDT",
    diag=True,
):
    fig, ax = plt.subplots(layout="constrained")
    h = ntuple.get(f"dataset/Method_{method}/{method}/MVA_{method}_rejBvsS").to_hist()

    if h is None:
        raise IndexError("could not find appropriate")

    integ = sum(h.values() * h.axes[0].widths)

    if diag:
        ax.plot([0.0, 1.0], [1.0, 0.0], color="lightgrey", linestyle="dashed")
    h.plot(ax=ax, yerr=False, label=f"BDT (AUC: {integ:.2f})")
    ax.set_title("Photon BDT", loc="left")
    ax.set_title("ROC Curve")
    ax.set_xlabel("Signal Efficiency")
    ax.set_ylabel("Background Rejection")
    ax.legend()
    if args.save:
        fig.savefig(f"MVA_ROC_curve.{args.ext}")
    plt.show()


def calc_significance(
    signal_eff: hist.BaseHist, bkg_eff: hist.BaseHist, n_signal=1000, n_bkg=1000
) -> hist.BaseHist:
    H_sig, _ = signal_eff.to_numpy()
    H_bkg, _ = bkg_eff.to_numpy()

    p = (H_sig * n_signal) / np.sqrt((H_bkg * n_bkg) + (H_sig * n_signal))

    H_significance = hist.Hist.new.Regular(
        signal_eff.axes[0].size,
        signal_eff.axes[0].edges[0],
        signal_eff.axes[0].edges[-1],
    ).Double()
    H_significance[...] = p

    optimum = H_significance.to_numpy()[1][np.argmax(H_significance.to_numpy()[0])]

    return H_significance, optimum


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="plot_bdt_outputs")
    parser.add_argument("filename", type=str)
    parser.add_argument("--roc", action="store_true", help="plot ROC curve")
    parser.add_argument(
        "--corr", action="store_true", help="plot input variable correlation matrices"
    )
    parser.add_argument(
        "--inputs", action="store_true", help="plot input variable distributions"
    )
    parser.add_argument(
        "--response", action="store_true", help="plot MVA response distributions"
    )
    parser.add_argument(
        "--perf",
        action="store_true",
        help="plot selection performance as a function of MVA response score",
    )
    parser.add_argument(
        "--save", action="store_true", help="save plot to the extension defined in ext"
    )
    parser.add_argument(
        "--ext", type=str, default="svg", help="file extension to save plot as"
    )

    args = parser.parse_args()

    main(args)
