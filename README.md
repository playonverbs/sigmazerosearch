# sigmazerosearch

> A loose python framework to help perform a particle physics selection for
> hyperons using MicroBooNE data files.

This project has dependencies managed by Poetry, use it to create a new
virtualenv to house dependencies for this project. 

Example usages (and what aims to be a source of top-down development) of the
framework are listed in the `examples/` directory.

Currently a `Selection` is only intended to be instantiated and run once in a
file, unexpected side-effects may occur when trying to define and use multiple
in the same script.

Tests are located in the `test/` directory and are run via `poetry run pytest
-v`. The test suite currently composes of unit tests only -- purpose-built
samples are required for full integration tests and as such integration tests
are not planned.

<details>
<summary>Algorithm Implementation Details</summary>

Most functions are written in as much of a vectorised way as possible, and will
return full arrays of values matching the maximum outer dimension of the
operand array.

</details>
