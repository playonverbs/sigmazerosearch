import pathlib

import pytest

from sigmazerosearch.general import PDG, Config
from sigmazerosearch.selection import Cut, Sample, SampleType, Selection

# TODO: split this file into separate files for each large type for more
# detailed testing, e.g. multiple tests for Selection


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


def test_PDG():
    assert PDG.Lambda.anti == 3122

    # neutral particles should return their own value
    assert PDG.Kaon0.anti == PDG.Kaon0.value
    assert PDG.Sigma0.anti == PDG.Sigma0.value
    assert PDG.Lambda.anti == PDG.Lambda.value
    assert ~PDG.Kaon == -321
    assert -PDG.Kaon == -321
    assert PDG.Sigma0 == 3212


@pytest.fixture
def eg_Cut():
    return Cut("bar", lambda: "hello")


def test_Cut(eg_Cut):
    c1 = Cut("this", lambda: 2)
    assert c1.name == "this"
    assert c1() == 2

    c2 = Cut("that", lambda: 3)
    assert c1.name != c2.name
    assert c2() == 3

    assert eg_Cut() == "hello"


def test_Sample():
    s1 = Sample(
        "hyperon",
        "analysisOutputRHC_Hyperon_numi_rhc_run3b.root",
        SampleType.Hyperon,
        1e20,
    )
    s2 = Sample(
        "dirt", "analysisOutputRHC_Dirt_numi_rhc_run3b.root", SampleType.Dirt, 2e30
    )

    assert s1.name != s2.name


def test_Selection():
    s = Selection(
        params={}, samples={}, cuts=[Cut("cut1", lambda: 2), Cut("cut2", lambda: 3)]
    )

    s.cuts[0].n_signal = (60, 10)
    s.cuts[0].n_passing = (40, 5)
    s.cuts[0].n_background = (1e6, 3e4)
    s.cuts[1].n_signal = (30, 10)
    s.cuts[1].n_passing = (15, 5)
    s.cuts[1].n_background = (1e4, 1e2)
