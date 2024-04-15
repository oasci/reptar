SHELL := /usr/bin/env bash
PYTHON_VERSION := 3.11
PYTHON_VERSION_CONDENSED := 311
PACKAGE_NAME := reptar
PACKAGE_PATH := $(PACKAGE_NAME)/
TESTS_PATH := tests/
CONDA_NAME := $(PACKAGE_NAME)-dev
CONDA := conda run -n $(CONDA_NAME)
CONDA_LOCK_OPTIONS := -p linux-64 -p osx-64 --channel conda-forge --channel conda-forge/label/libint_dev

###   ENVIRONMENT   ###

# See https://github.com/pypa/pip/issues/7883#issuecomment-643319919
export PYTHON_KEYRING_BACKEND := keyring.backends.null.Keyring

.PHONY: conda-create
conda-create:
	- conda deactivate
	conda remove -y -n $(CONDA_NAME) --all
	conda create -y -n $(CONDA_NAME)
	$(CONDA) conda install -y python=$(PYTHON_VERSION)
	$(CONDA) conda install -y conda-lock

# Default packages that we always need.
.PHONY: conda-setup
conda-setup:
	$(CONDA) conda install -y -c conda-forge poetry
	$(CONDA) conda install -y -c conda-forge pre-commit
	$(CONDA) conda install -y -c conda-forge tomli tomli-w
	$(CONDA) conda install -y -c conda-forge conda-poetry-liaison

# Conda-only packages specific to this project.
.PHONY: conda-dependencies
conda-dependencies:
	$(CONDA) conda install -y -c conda-forge xtb
	$(CONDA) conda install -y -c conda-forge/label/libint_dev -c conda-forge psi4

.PHONY: nodejs-dependencies
nodejs-dependencies:
	$(CONDA) conda install -y -c conda-forge nodejs
	$(CONDA) npm install markdownlint-cli2 --global

# Add custom command to add unique channel
.PHONY: conda-lock
conda-lock:
	- rm conda-lock.yml
	$(CONDA) conda env export --from-history | grep -v "^prefix" > environment.yml
	$(CONDA) conda-lock -f environment.yml $(CONDA_LOCK_OPTIONS)
	$(CONDA) cpl-deps pyproject.toml --env_name $(CONDA_NAME)
	$(CONDA) cpl-clean --env_name $(CONDA_NAME)

.PHONY: from-conda-lock
from-conda-lock:
	$(CONDA) conda-lock install -n $(CONDA_NAME) conda-lock.yml
	$(CONDA) cpl-clean --env_name $(CONDA_NAME)

.PHONY: pre-commit-install
pre-commit-install:
	$(CONDA) pre-commit install

# Reads `pyproject.toml`, solves environment, then writes lock file.
.PHONY: poetry-lock
poetry-lock:
	$(CONDA) poetry lock --no-interaction

.PHONY: install
install:
	$(CONDA) poetry install --no-interaction
	- mkdir .mypy_cache
	- $(CONDA) mypy --install-types --non-interactive --explicit-package-bases $(PACKAGE_NAME)

.PHONY: environment
environment: conda-create nodejs-dependencies from-conda-lock pre-commit-install install

.PHONY: locks
locks: conda-create conda-setup conda-dependencies nodejs-dependencies conda-lock pre-commit-install poetry-lock install


###   TESTING   ###

.PHONY: download-test-files
download-test-files:
	$(CONDA) git clone https://github.com/oasci/reptar-data examples/reptar-data/

.PHONY: test
test:
	$(CONDA) pytest -c pyproject.toml --cov=$(PACKAGE_NAME) --cov-report=xml --junit-xml=report.xml --color=yes $(TESTS_PATH)

.PHONY: coverage
coverage:
	$(CONDA) coverage report


###   FORMATTING   ###

.PHONY: validate
validate:
	- $(CONDA) markdownlint-cli2 "**.md" --config ./.markdownlint.yaml --fix
	- $(CONDA) pre-commit run --all-files

.PHONY: formatting
formatting:
	- $(CONDA) isort --settings-path pyproject.toml ./
	- $(CONDA) black --config pyproject.toml ./



###   LINTING   ###

.PHONY: check-codestyle
check-codestyle:
	$(CONDA) isort --diff --check-only $(PACKAGE_PATH)
	$(CONDA) black --diff --check --config pyproject.toml $(PACKAGE_PATH)
	- $(CONDA) pylint --rcfile pyproject.toml $(PACKAGE_PATH)

.PHONY: mypy
mypy:
	- $(CONDA) mypy --config-file pyproject.toml $(PACKAGE_PATH)

.PHONY: lint
lint: check-codestyle mypy




###   CLEANING   ###

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

.PHONY: psi-remove
psi-remove:
	find . | grep -E ".clean" | xargs rm -rf
	find . | grep -E "timer.dat" | xargs rm -rf

.PHONY: coverage-remove
coverage-remove:
	find . | grep -E ".coverage" | xargs rm -rf

.PHONY: build-remove
build-remove:
	rm -rf build/

.PHONY: cleanup
cleanup: pycache-remove dsstore-remove mypycache-remove ipynbcheckpoints-remove pytestcache-remove psi-remove coverage-remove


###   BUILD   ###

.PHONY: build
build:
	$(CONDA) poetry build



###   DOCS   ###

mkdocs_port := $(shell \
	start_port=3000; \
	max_attempts=100; \
	for i in $$(seq 0 $$(($$max_attempts - 1))); do \
		current_port=$$(($$start_port + i)); \
		if ! lsof -i :$$current_port > /dev/null; then \
			echo $$current_port; \
			break; \
		fi; \
		if [ $$i -eq $$(($$max_attempts - 1)) ]; then \
			echo "Error: Unable to find an available port after $$max_attempts attempts."; \
			exit 1; \
		fi; \
	done \
)

.PHONY: serve
serve:
	echo "Served at http://127.0.0.1:$(mkdocs_port)/"
	$(CONDA) mkdocs serve -a localhost:$(mkdocs_port)

.PHONY: docs
docs:
	$(CONDA) mkdocs build -d public/
	- rm -f public/gen_ref_pages.py

.PHONY: open-docs
open-docs:
	xdg-open public/index.html 2>/dev/null
