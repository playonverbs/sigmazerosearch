from sigmazerosearch.types.selection import Selection, Cut, Sample, SampleType
from sigmazerosearch.types.general import PDG
from sigmazerosearch.types.truth import GenType, GenEventType

# TODO: split this file into separate files for each large type for more
# detailed testing, e.g. multiple tests for Selection


def test_PDG():
    assert PDG.Lambda.anti == 3122

    # neutral particles should return their own value
    assert PDG.Kaon0.anti == PDG.Kaon0.value
    assert PDG.Sigma0.anti == PDG.Sigma0.value
    assert PDG.Lambda.anti == PDG.Lambda.value
    assert ~PDG.Kaon == -321
    assert -PDG.Kaon == -321


def test_Cut():
    c1 = Cut("this", "null")
    assert c1.name == "this"
    assert c1.cutexpr == "null"

    c2 = Cut("that", "!null")
    assert c1.name != c2.name
    assert c1.cutexpr != c2.cutexpr


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
    s = Selection()
    s.cuts = [Cut("cut1", "inFV"), Cut("cut2", "inFV && NuMu")]

    s.cuts[0].n_signal = (60, 10)
    s.cuts[0].n_passing = (40, 5)
    s.cuts[0].n_background = (1e6, 3e4)
    s.cuts[1].n_signal = (30, 10)
    s.cuts[1].n_passing = (15, 5)
    s.cuts[1].n_background = (1e4, 1e2)

    t = Selection.with_cuts([Cut("cut1", "inFV")])
    assert t.cuts[0].name == s.cuts[0].name
