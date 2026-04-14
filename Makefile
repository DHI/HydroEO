LIB = HydroEO

check: lint typecheck test doctest

build: typecheck test
	python -m build

lint:
	uvx ruff check $(LIB)

format:
	uvx ruff format $(LIB)

test:
	uv run pytest --disable-warnings

typecheck:
	uvx mypy $(LIB)/ --config-file pyproject.toml

doctest:
	uvx pytest --doctest-modules $(LIB)

coverage: 
	uvx pytest --cov-report html --cov=$(LIB) tests/

docs: FORCE
	mkdocs build

clean:
	python -c "import shutil; shutil.rmtree('dist', ignore_errors=True)"
	python -c "import shutil; shutil.rmtree('htmlcov', ignore_errors=True)"
	python -c "import os; os.remove('.coverage') if os.path.exists('.coverage') else None"
	python -c "import shutil; shutil.rmtree('site', ignore_errors=True)"

FORCE: