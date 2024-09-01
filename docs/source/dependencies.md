# Dependencies

Major project dependencies are as follows:

[Uproot](<inv:uproot#index>)
: Uproot provides Python-native I/O for ROOT objects. In particular we rely on
  reading the input ntuple files that constitute processed MicroBooNE data.

[Awkward Array](<inv:awkward-array#index>)
: Handles multi-dimensional arrays with ragged dimensions. Essential for an
  array where the outermost dimension corresponds to an event and inner
  dimensions can correspond to the number of truth-level particles,
  reconstructed objects, etc.

[Matplotlib](<inv:matplotlib#index>)
: The de-facto plotting library for Python.

[Hist](<inv:hist#index>)
: Provides helpful objects to handle histogramming analysis values. Simplifies
  a lot of boilerplate code that might be needed for histogramming.

[Tabulate](<pypi:tabulate>)
: Pretty looking tables printed on the command-line. Allows for outputting in
  multiple formats including Markdown and LaTeX.
