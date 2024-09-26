import pathlib

import pytest

from sigmazerosearch.general import PDG, Config


@pytest.fixture
def default_Config():
    return Config(plot_save=False, plot_dir=None, branch_list=None)


def test_Config(default_Config):
    assert Config.default().branch_list == default_Config.branch_list
    assert Config.default().plot_dir == default_Config.plot_dir
    assert Config.default().plot_save == default_Config.plot_save

    assert Config.default() == default_Config

    with pytest.raises(Exception):
        Config(plot_save=True, plot_dir=None)

    with pytest.raises(Exception):
        Config(plot_save=True, plot_dir=pathlib.Path("./test_utils.py"))

    with pytest.raises(Exception):
        Config(iterate=False, iterate_step=100)

    with pytest.raises(Exception):
        Config(iterate=True, iterate_step=None)

    with pytest.raises(Exception):
        Config(iterate=False, iterate_step="1 MB")


def test_PDG():
    assert PDG.Lambda.anti == 3122

    # neutral particles should return their own value
    assert PDG.Kaon0.anti == PDG.Kaon0.value
    assert PDG.Sigma0.anti == PDG.Sigma0.value
    assert PDG.Lambda.anti == PDG.Lambda.value
    assert ~PDG.Kaon == -321
    assert -PDG.Kaon == -321
    assert PDG.Sigma0 == 3212
