"""
Methods that operate on physics information to aid in selections.

Many of the methods here accept
{py:obj}`sigmazerosearch.selection.ParameterSet` as an argument which contains
values that affect the physics decisions of this framework.

Submodules are organised by what problems they aim to solve:

- `fv`: provides a fiducial volume cut and some other definitions regarding the
  detector geometry.
- `muon`: provides methods of selecting a muon-like object.
- `lamb`: methods for selecting the $\Lambda \\rightarrow p + \pi^-$ subsystem.
"""
