import matplotlib.pyplot as plt
import numpy as np
import uproot as up

from sigmazerosearch import defaults, utils

plt.style.use("./plots/sigmazerosearch.tex.mplstyle")


def main():
    fd = up.open("./scripts/calc_xs.root")

    gr_posterior = fd.get("posterior").values()
    gr_ci = fd.get("conf_interval").values()

    n_obs = fd.get("param_N_obs").value
    n_bkg = fd.get("param_N_bkg").value
    efficiency = fd.get("param_efficiency").value
    cl = fd.get("param_CL").value * 100

    lul_bayes = fd.get("param_lUL_Bayes").value / 1e-40
    lul_fc = fd.get("param_lUL_FC").value / 1e-40
    mc_ref = fd.get("param_MC_ref").value / 1e-40

    fig, ax = plt.subplots(layout="constrained")

    ax.plot(*gr_posterior, c="tab:blue")
    ax.fill_between(
        gr_ci[0],
        np.interp(gr_ci[0], *gr_posterior),
        fc="tab:blue",
        alpha=0.4,
        label=rf"{cl:.0f}\% CL region",
    )

    ax.axvline(
        lul_bayes,
        alpha=1,
        label=rf"$\sigma_\mathrm{{UL}}^\mathrm{{Bayes}} = {lul_bayes:0.2f} \times 10^{{-40}} \mathrm{{cm}}^2/\mathrm{{Ar}}$ ({cl:.0f}\% CL)",
        c="tab:blue",
    )

    ax.axvline(
        lul_fc,
        alpha=0.7,
        label=rf"$\sigma_\mathrm{{UL}}^\mathrm{{FC}} = {lul_fc:0.2f} \times 10^{{-40}} \mathrm{{cm}}^2/\mathrm{{Ar}}$ (stat.)",
        c="red",
        ls="dotted",
    )

    ax.axvline(
        mc_ref,
        alpha=0.7,
        label=rf"$\sigma_\mathrm{{MC\;Ref}} = {mc_ref:0.2f} \times 10^{{-40}} \mathrm{{cm}}^2/\mathrm{{Ar}}$",
        c="orange",
        ls="dashed",
    )

    ax.set_xlabel(
        r"$\sigma ( \overline{\nu}_\mu + \mathrm{Ar} \rightarrow \mu^+ + \Sigma^0 + X )$ "
        r"[$10^{-40} \mathrm{cm}^2/\mathrm{Ar}$]"
    )
    ax.set_ylabel(r"Posterior PDF [$10^{40}\;\mathrm{cm}^2 \mathrm{Ar}^{-1}$]")

    ax.set_title("NuMI Run 3b (RHC)")

    ax.set_ylim(0)

    ax.legend(
        title=rf"$N_\mathrm{{obs}} = {n_obs}$, $N_\mathrm{{bkg}} = {n_bkg}$, $\epsilon = {efficiency:0.1f}\%$",
        title_fontsize=13,
        fontsize=13,
    )
    utils._save_plot(defaults.config, fig, "xs/xs_estimation")
    plt.show()


if __name__ == "__main__":
    main()
