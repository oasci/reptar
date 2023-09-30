(setup-environment)=
# Setup environment

A crucial aspect of consistent development is the standardization of the computational environment.
We use a combination of [`conda`](https://conda.io/) and [`poetry`](https://python-poetry.org/) to accomplish this.

However, mixing environment managers like `conda` and `poetry` must be done carefully.
This usually involves installing all desired `conda` packages and using only poetry afterward.
If you want to use a new `conda` package down the road, you usually need to recreate the environment from scratch.

## Steps

### Installing conda

If you do not have `conda` installed, follow the instructions [here](https://docs.conda.io/projects/miniconda/en/latest/#quick-command-line-install).

:::{note}
We recommend using [`libmamba`](https://conda.github.io/conda-libmamba-solver/getting-started/) instead of [`mamba`](https://mamba.readthedocs.io/en/latest/) or [classic `conda`](https://conda.github.io/conda-libmamba-solver/libmamba-vs-classic/).
:::

### Conda environment

First, we setup a `conda` environment inside the repository (`.venv`).

```bash
make conda-setup
```

::::{tab-set}

:::{tab-item} From `conda-lock.yml`

```bash
make from-conda-lock
```

```{warning}
Only `linux-64` and `osx-64` are supported.
If you must use windows, it is recommended that you build from scratch.
```

:::

:::{tab-item} From scratch

Install the required conda packages.

```bash
make conda-dependencies
```

If you desire more conda packages, activate the conda environment first.

```bash
conda activate ./.venv
```

Add all other conda channels so they are exported to `environment.yml`.
`conda-forge` is already included in `make conda-dependencies`.

```bash
conda config --add channels conda-forge
```

Install all desired conda packages; for example,

```bash
conda install openmm -c conda-forge
```

If needed, write a new `conda-lock` file.

```bash
make write-conda-lock
```

:::

::::

<!-- conda list -e -p .venv/ | sed '1,3d ; s/=[a-z][A-Za-z0-9].*$//g ; s/^_[A-Za-z0-9].*$//g ; s/=/==/g' -->

### Poetry-tracked packages

After installing all `conda` packages, we switch over to exclusively using `poetry`.
The following command uses `poetry` to install all packages specified in `pyproject.toml`.

```bash
make install
```

### `pre-commit`

TODO:

```bash
make pre-commit-install
```
