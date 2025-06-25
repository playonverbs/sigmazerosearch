.PHONY: *

test:
	uv run pytest -v

docs:
	uv run sphinx-autobuild docs/source docs/build/html --watch sigmazerosearch/

lint:
	uv run ruff check .

format:
	uv run ruff format .
