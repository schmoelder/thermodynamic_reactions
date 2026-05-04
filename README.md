# Thermodynamic Reactions

A thermodynamically consistent reaction library for CADET, developed alongside an educational JupyterBook.

## Prerequisites

- Python ≥ 3.9

## Install

```bash
pip install -e .
```

With optional dependency groups:

```bash
pip install -e . --group test   # adds pytest
pip install -e . --group docs   # adds jupyter-book
```

Run the test suite:

```bash
pytest tests/
```

## Build the book

```bash
cd jupyterbook
myst build --html
```

The output is written to `jupyterbook/_build/html/`.

## Repository layout

```
reactions/          # installable Python library
jupyterbook/        # MyST JupyterBook source
tests/              # test suite
pyproject.toml      # package metadata
```
