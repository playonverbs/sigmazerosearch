# Data Structures

This document details data structures that exist semantically in the framework,
but are not structured as classes or dataclasses for reasons of program speed.

These are primarily within <inv:#numpy.array> or <inv:#ak.Array> objects.

## Island Window

An island window is a 2-dimensional array of integer values that aims to
represent pixellated ionisation activity on a single wire plane. Coordinates
correspond to wire, time coordinates. $(x,y) \equiv
(\mathrm{time},\mathrm{wire})$.
The background of no ionisation activity is represented by a `0` value, while
any other positive value represents ionisation.

```
╭───────────────────────╮
│ 0 | 1 | 0 | 0 | 0 | 0 │
│───────────────────────│
│ 0 | 1 | 1 | 0 | 0 | 0 │
│───────────────────────│
│ 0 | 0 | 1 | 0 | 0 | 0 │
│───────────────────────│
│ 0 | 0 | 1 | 0 | 0 | 0 │
│───────────────────────│
│ 0 | 0 | 1 | 0 | 0 | 0 │
╰───────────────────────╯
```

The diagram above as an island window would encompass a window spanning `x = 5`
wires and the entire time axis binned into `y = 6` bins.

### Construction

Island windows start as a more compact `std::vector<std::vector<float>>`. The
dimensions of the outermost `vector` represent the number of wires that
comprise the relative wire window placed around the reconstructed neutrino
vertex position. The inner vector dimensions represent the number of hits on
that particular wire. Floats within the innermost vector are the resolved peak
hit times of `recob::Hits` in the neutrino slice.

This initial structure was chosen due to limitations in outputting more complex
objects to ROOT TTrees in the LArSoft framework used by the upstream
<github:playonverbs/HyperonProduction> module.

The functions <project:#alg.island._window_to_map> and
<project:#alg.island._window_to_map_numpy> both construct a window from the
compact format.

### Variants

#### Binary Island Map

The binary island window is an island window with elements that are only either
`0` or `1`. This represents the simplest level of information stating if a
location contains ionisation activity.

- `0` represents a pixel with no activity -- or the 'background'.
- `1` represents any ionisation activity.

Binary island maps are produced by the <project:#alg.island._window_to_map> and
<project:#alg.island._window_to_map_numpy> functions. The functions that take
them as inputs are <project:#alg.island._find_map_islands>.

#### Segmented Island Map
The segmented island map is produced after the primary island finding stage and
represents ionisation activity separated by their proximity to other connected
islands of activity. It has the same shape as a standard island window but
instead elements representing activity can be any positive, non-zero integer.

- `0` represents a pixel with no activity -- the 'background'.
- `1..n` represents a pixel with activity. Elements with the same number are
  considered adjacent to one another as per the procedure used by the
  underlying <inv:#scipy.ndimage.label>.
