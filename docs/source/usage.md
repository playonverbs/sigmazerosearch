# Usage

The framework revolves around declaring a <project:#Selection>. This object
encapsulates the sample files that comprise a selection and the state at each
stage (or Cut) as the selection progresses.

The object has several inputs for the selection:

:::{list-table}
:align: left
:widths: 25 30 15 50
:header-rows: 1

*   - Input
    - Collection
    - Single
    - Note
*   - Sample
    - <project:#SampleSet>
    - <project:#Sample>
    - Represents an NTuple produced by <github:playonverbs/HyperonProduction>.
      Also manages sample type and POT.
*   - Cut
    - <inv:python:py:*:#list>[<project:#Cut>]
    - <project:#Cut>
    - Stores data required for selection performance calculations.
*   - ParameterSet
    - <project:#ParameterSet>
    - %
    - Holds selection-dependent parameters, such as the values that certain
      cuts depend on.
:::

Alongside these the <project:#Config> class keeps options that concern how the
selection operates.

## Creating a Selection

:::{attention}
A <project:#ParameterSet> must be defined before the main body of the Selection
object.

This is because while defining parameters within a class instance, the other
parameters within the class cannot be referenced until the entire object has
been resolved.
:::

```{code-block} python
:linenos:

pset = ParameterSet.from_dict({...})

sel = Selection(
    params=pset,
    cuts=[
        Cut(
            "fv",
            lambda arr: fv.in_active_tpc(arr)
        )
    ],
    samples=SampleSet(
        Sample(
            "hyperon",
            "~/mydata/hyperon.root",
            SampleType.Hyperon,
            None
        ),
        Sample(
            "background",
            "~/mydata/background.root",
            SampleType.Background,
            None
        ),
        target_POT=2E20
    )
)
```

### Defining your own cuts

The <project:#Cut> objects used in the framework support the already defined
selection cut methods using `lambda`-currying. This is needed as many selection
functions take parameters via <project:#ParameterSet>; however the object
provided to the Cut object must have the function signature:
{external:class}`ak.Array` &rarr; {external:class}`ak.Array`. Wrapping
functions allows this.

When defining your own cuts the primary requirement is that the cut function
returns a 1-dimensional <inv:#ak.Array> consisting of <inv:#bool> with a length
matching the number of events. In the
[Datashape](https://datashape.readthedocs.io/en/latest/overview.html) language,
for 100,000 events this is represented as `100000 * bool`.

```python
Cut(
    "foo",
    lambda arr: pass
)
```

### Extending the Available Parameters

Your own cuts may require additional parameters that are not supplied by
<project:#ParameterSet>. To add your own parameters you should create a new
frozen dataclass that inherits from the existing one like so:

```python
from dataclasses import dataclass

from sigmazerosearch.selection import ParameterSet

@dataclass(frozen=True)
class MyParameterSet(ParameterSet):
    a_new_param: float
    """Used in <x> cut, units are GeV"""
    flag_param: bool
    """Flag for <cut>"""
```

Here we can note a few things:
1. Fields should not be defined with default values.
2. Fields should be documented with a minimum mention of the field's physical
  units (if any).

Now the extended ParameterSet can be used in your own and pre-existing cut
functions like so:

```python
import awkward as ak

def select_foos(arr: ak.Array, my_pset: MyParameterSet):
    if my_pset.flag_param:
        return arr < my_pset.a_new_param

    return arr > my_pset.a_new_param
```
