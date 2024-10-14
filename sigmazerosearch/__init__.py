"""
A framework for operating hyperon-Ar cross-section analyses on
[MicroBooNE](https://microboone.fnal.gov/) data.

This module encapsulates a framework for operating a particle physics selection
focussed on measuring a cross-section for Hyperon-Argon interactions using
MicroBooNE Monte-Carlo and Data files.

Many references are made to objects produced from the
[LArSoft](https://code-doc.larsoft.org/docs/) software framework for
Liquid-Argon TPC experiments. Readers are also assumed to have knowledge of
well-known particle physics terminology.

> The main entry point for the framework is the
> `sigmazerosearch.selection.Selection` object.

The input `Sample` files to this framework are produced by the
HyperonProduction LArSoft analyser located at
<https://github.com/playonverbs/HyperonProduction>.

An occasionally-updated spec for the input files is available
[here](https://hep.lancs.ac.uk/hedgedoc/s/FkkYg14yW).

:::{warning}
This project is under very active development. Expect breaking API changes.
:::
"""
