from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np

from sigmazerosearch import defaults, utils

plt.style.use("./plots/sigmazerosearch.tex.mplstyle")


@dataclass
class Params:
    # p_thresh_p: float = 300.0 # MeV
    p_thresh_p: float = 0.3  # GeV
    # p_thresh_pi: float = 100.0 # MeV
    p_thresh_pi: float = 0.1  # GeV


def lambda_fraction(mom_lambda: float, params: Params = Params()) -> float:
    """
    Calculate the modified phase space of lambda production as viewed with
    detection thresholds on the proton and pion.
    """

    m_lambda = 1.115683  # GeV
    m_p = 0.938272  # GeV
    m_pi = 0.139570  # GeV
    p = 0.101  # GeV
    p_thresh_p = params.p_thresh_p
    p_thresh_pi = params.p_thresh_pi
    E_p = np.sqrt((m_p**2) + (p**2))
    E_pi = np.sqrt((m_pi**2) + (p**2))

    beta = mom_lambda / np.sqrt((mom_lambda**2) + (m_lambda**2))
    gamma = 1 / np.sqrt(1 - (beta**2))

    # print(f"p_lambda = {mom_lambda:.2f}\tbeta = {beta:.2f}\tgamma = {gamma:.2f}")

    A = max(
        (
            (np.sqrt((m_p**2) + (abs(p_thresh_p) ** 2)) - (gamma * E_p))
            / (beta * gamma * p)
        ),
        -1.0,
    )

    B = min(
        (
            (-np.sqrt((m_pi**2) + (abs(p_thresh_pi) ** 2)) + (gamma * E_pi))
            / (beta * gamma * p)
        ),
        1.0,
    )

    # print(f"A = {A:.2f}\tB = {B:.2f}")

    if A > B:
        return 0
    else:
        return 0.5 * (B - A)


def main():
    inputs = np.arange(0, 2.5, 0.001)

    outputs = []

    for i in inputs:
        outputs.append(lambda_fraction(i))

    fig, ax = plt.subplots(layout="constrained")

    ax.axhline(1.0, c="lightgray", linestyle="dashed")
    ax.plot(inputs, outputs)

    ax.set_title(r"Modified $\Lambda$ phase-space")

    ax.text(
        0.975,
        0.05,
        f"$p^\mathrm{{thresh}}_p = {Params().p_thresh_p}$ GeV\n$p^\mathrm{{thresh}}_\pi = {Params().p_thresh_pi}$ GeV",
        transform=ax.transAxes,
        ha="right",
    )

    ax.set_xlabel(r"$p_{\Lambda}$ [GeV]")
    ax.set_ylabel(r"$f(p_{\Lambda})$")

    utils._save_plot(defaults.config, fig, "partial_phase_space_lambda")

    plt.show()


if __name__ == "__main__":
    main()
