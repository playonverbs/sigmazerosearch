import matplotlib.pyplot as plt
import numpy as np
import pytest

from sigmazerosearch.alg import syst


@pytest.fixture
def p_vecs():
    return np.random.uniform(low=-5.0, high=20.0, size=(1000000, 3))


def test_plot_get_numi_angle(p_vecs):
    angles = syst.get_numi_angle(p_vecs[:, 0], p_vecs[:, 1], p_vecs[:, 2], "beam")
    target_angles = syst.get_numi_angle(
        p_vecs[:, 0], p_vecs[:, 1], p_vecs[:, 2], "target"
    )

    binscheme = [0, 10, 20, 110, 160]
    _, ax = plt.subplots(layout="constrained")

    ax.hist(angles, binscheme, histtype="step", label="Beam")
    ax.hist(target_angles, binscheme, histtype="step", label="Target")

    ax.set_xlabel(r"$\theta_{\nu - X}$")
    ax.set_title(r"Detector $\rightarrow$ NuMI angle")
    ax.legend()
    plt.show()

    print(syst.get_numi_angle([0.0, 1.0], [1.0, 0.0], [2.0, 4.0], "beam"))
    print(syst.get_numi_angle([0.0, 1.0], [1.0, 0.0], [2.0, 4.0], "target"))
