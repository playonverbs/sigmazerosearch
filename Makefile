.PHONY: *

test:
	poetry run pytest -v

docs:
	poetry run sphinx-autobuild docs/source docs/build/html

lint:
	poetry run ruff check .

format:
	poetry run ruff format .
