#* Variables
SHELL := /usr/bin/env bash
PYTHON_VERSION := 3.11
PYTHONPATH := `pwd`
DOCS_URL := https://reptar.oasci.org
REPO_PATH := $(shell git rev-parse --show-toplevel)
CONDA := conda run -p $(REPO_PATH)/.venv

#* Setup
.PHONY: conda-setup
conda-setup:
	conda create -y --prefix $(REPO_PATH)/.venv python=$(PYTHON_VERSION)
	conda install -y conda-lock conda-libmamba-solver -p $(REPO_PATH)/.venv
	$(CONDA) conda config --set solver libmamba
	conda install -y -c conda-forge poetry pre-commit -p $(REPO_PATH)/.venv

.PHONY: conda-dependencies
conda-dependencies:
	$(CONDA) conda config --add channels conda-forge
	$(CONDA) conda config --add channels conda-forge/label/libint_dev
	$(CONDA) conda install psi4 python=$(PYTHON_VERSION) -c conda-forge/label/libint_dev -c conda-forge -y
	$(CONDA) conda install xtb -c conda-forge

.PHONY: pre-commit-install
pre-commit-install:
	$(CONDA) pre-commit install

.PHONY: from-conda-lock
from-conda-lock:
	$(CONDA) conda-lock install -p $(REPO_PATH)/.venv $(REPO_PATH)/conda-lock.yml

.PHONY: write-conda-lock
write-conda-lock:
	$(CONDA) conda env export --from-history | grep -v "^prefix" > environment.yml
	$(CONDA) conda-lock -f environment.yml -p linux-64 -p osx-64



#* Installation
.PHONY: install
install:
	$(CONDA) poetry lock --no-interaction
	$(CONDA) poetry export --without-hashes > requirements.txt
	$(CONDA) poetry install --no-interaction
	-$(CONDA) mypy --install-types --non-interactive ./reptar

.PHONY: update-poetry
update-poetry:
	$(CONDA) poetry update --no-interaction

.PHONY: update
update: update-poetry install


.PHONY: validate
validate:
	$(CONDA) pre-commit run --all-files



#* Formatters
.PHONY: codestyle
codestyle:
	$(CONDA) pyupgrade --exit-zero-even-if-changed --py311-plus **/*.py
	$(CONDA) isort --settings-path pyproject.toml ./
	$(CONDA) black --config pyproject.toml ./

.PHONY: formatting
formatting: codestyle




#* Linting
.PHONY: test
test:
	PYTHONPATH=$(PYTHONPATH) $(CONDA) pytest -c pyproject.toml --cov=reptar --cov-report=xml tests/

.PHONY: check-codestyle
check-codestyle:
	$(CONDA) isort --diff --check-only --settings-path pyproject.toml ./
	$(CONDA) black --diff --check --config pyproject.toml ./
	-$(CONDA) pylint reptar

.PHONY: mypy
mypy:
	$(CONDA) mypy --config-file pyproject.toml ./

.PHONY: check-safety
check-safety:
	$(CONDA) poetry check
	$(CONDA) safety check --full-report
	$(CONDA) bandit -ll --recursive reptar tests

.PHONY: lint
lint: test check-codestyle mypy check-safety



#* Cleaning
.PHONY: pycache-remove
pycache-remove:
	find . | grep -E "(__pycache__|\.pyc|\.pyo$$)" | xargs rm -rf

.PHONY: dsstore-remove
dsstore-remove:
	find . | grep -E ".DS_Store" | xargs rm -rf

.PHONY: mypycache-remove
mypycache-remove:
	find . | grep -E ".mypy_cache" | xargs rm -rf

.PHONY: ipynbcheckpoints-remove
ipynbcheckpoints-remove:
	find . | grep -E ".ipynb_checkpoints" | xargs rm -rf

.PHONY: pytestcache-remove
pytestcache-remove:
	find . | grep -E ".pytest_cache" | xargs rm -rf

.PHONY: build-remove
build-remove:
	rm -rf build/

.PHONY: cleanup
cleanup: pycache-remove dsstore-remove mypycache-remove ipynbcheckpoints-remove pytestcache-remove



#* Documentation
.PHONY: docs
docs:
	rm -rf ./docs/html/
	$(CONDA) sphinx-build -nT ./docs/source/ ./docs/html/
	touch ./docs/html/.nojekyll

.PHONY: open-docs
open-docs:
	xdg-open ./docs/html/index.html 2>/dev/null

.PHONY: docs-multiversion
docs-multiversion:
	rm -rf ./docs/html/
	$(CONDA) sphinx-multiversion -nT ./docs/source/ ./docs/html/
	touch ./docs/html/.nojekyll

	# Create html redirect to main
	echo "<head>" > ./docs/html/index.html
	echo "  <meta http-equiv='refresh' content='0; URL=$(DOCS_URL)/main/index.html'>" >> ./docs/html/index.html
	echo "</head>" >> ./docs/html/index.html
