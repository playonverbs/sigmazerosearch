# Contributing

## Tools

- [uv](https://docs.astral.sh/uv/) manages dependencies and the Python
  virtualenv.
- [pytest](https://pytest.org/) is the test handler for those in `tests/`
- [ruff](https://docs.astral.sh/ruff) is the primary linter and formatter.
- [pre-commit](https://github.com/pre-commit/pre-commit) runs the linting and
  formatting steps before a commit.

## Writing Documentation

Documentation is generated using [Sphinx](https://sphinx-doc.org). Markup files
are read from the `docs/source/` directory and API documentation is built from
the `sigmazerosearch/` directory.

You can generate documentation from markup and Python docstrings and serve it
locally using `make docs`.
