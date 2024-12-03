```{eval-rst}
.. sigmazerosearch documentation master file, created by
   sphinx-quickstart on Sat Aug 31 01:13:25 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
```

# sigmazerosearch

A framework for performing hyperon-argon cross-section analyses using
data simulated and recorded from the
[MicroBooNE](https://microboone.fnal.gov "Micro Booster Neutrino Experiment")
experiment.

This project aims to:

- Support the selection of hyperon production events as seen inside Liquid
  Argon Time Projection Chambers, with a preference towards $\Sigma^0$ direct
  production events.
- Characterise both track and shower-like objects identified with the Pandora
  multi-algorithm pattern recognition framework.
- Provide reusuability beyond MicroBooNE-specific tooling -- looking
  towards near and future experiments using the
  {abbr}`LArSoft (Liquid Argon Software)` framework.

```{card}
:link: https://github.com/cthorpe123/HyperonAnalysis

The author is grateful for and takes inspiration from
<github:cthorpe123/HyperonAnalysis>.
```

---

This project currently supports Python 3.10, more information on the supporting
packages is available in <project:#dependencies>.

```{toctree}
:maxdepth: 2
:caption: Contents:
:hidden:

usage
dependencies
api/index
larsoft
ntuples
```

```{toctree}
:caption: Development
:hidden:

contributing/index
```
