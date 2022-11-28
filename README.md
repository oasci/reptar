<h1 align="center">reptar</h1>

<h4 align="center">A tool for storing and analyzing manuscript-scale computational chemistry data</h4>

<h4 align="center" style="padding-bottom: 0.5em;"><a href="https://www.aalexmmaldonado.com/reptar/">Documentation</a></h4>

<p align="center">
    <a href="https://github.com/aalexmmaldonado/reptar/actions/workflows/python-package.yml">
        <img src="https://github.com/aalexmmaldonado/reptar/actions/workflows/python-package.yml/badge.svg" alt="Build Status ">
    </a>
    <a href="https://codecov.io/gh/aalexmmaldonado/reptar">
        <img src="https://codecov.io/gh/aalexmmaldonado/reptar/branch/main/graph/badge.svg?token=74wLrsOMTD" alt="codecov">
    </a>
    <a href="https://github.com/aalexmmaldonado/reptar/releases">
        <img src="https://img.shields.io/github/v/release/aalexmmaldonado/reptar" alt="GitHub release (latest by date)">
    </a>
    <a href="https://doi.org/10.5281/zenodo.6508586">
        <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.6508586.svg" alt="DOI">
    </a>
    <a href="https://github.com/aalexmmaldonado/reptar/blob/main/LICENSE" target="_blank">
        <img src="https://img.shields.io/github/license/aalexmmaldonado/reptar" alt="License">
    </a>
    <a href="https://github.com/aalexmmaldonado/reptar/" target="_blank">
        <img src="https://img.shields.io/github/repo-size/aalexmmaldonado/reptar" alt="GitHub repo size">
    </a>
    <a href="https://github.com/psf/black" target="_blank">
        <img src="https://img.shields.io/badge/code%20style-black-000000.svg" alt="Black style">
    </a>
</p>

<p align="center">
    <a href="#motivation">Motivation</a> •
    <a href="#installation">Installation</a> •
    <a href="#file-types">File Types</a> •
    <a href="#key-value-pairs">Key-value pairs</a> •
    <a href="#workflow">Workflow</a> •
    <a href="#license">License</a>
</p>

<p align="center"><b>Disclaimer:</b> reptar is under active development and is not suitable for production.</p>

# Motivation

The computational chemistry community often fails to openly provide raw and/or processed data used to draw their scientific conclusions.

For large projects, frameworks such as [QCArchive](https://qcarchive.molssi.org/), [Materials Project](https://materialsproject.org/), [Pitt Quantum Repository](https://pqr.pitt.edu/), [ioChem-BD](https://www.iochem-bd.org/) and many others provide great storage solutions.
This approach would not be practical for fluid data pipelines and small-scale projects such as a single manuscript.

Alternatively, you could use individual files in formats such as JSON, XML, YAML, npz, etc.
These are great options for customizable data storage with their own advantages and disadvantages.
However, you often must choose between (1) a standardized parser that might not support your workflow or (2) writing your own.

Reptar is designed for easy data storage and analysis for individual projects.
Customizable parsers provide a simple way to extract new data without submitting issues and pull requests (although this is highly encouraged).
While files are the heart of reptar, it strives to be file-type agnostic by providing the same interface for all supported file types.
The result is a user-specified file streamlined for analysis in Python and archival on places such as [GitHub](https://github.com/) and [Zenodo](https://zenodo.org/).

# Installation

You can install reptar from [PyPI](https://pypi.org/project/reptar/) by using `pip install reptar`.
Or, the latest development version can be installed directly from the [GitHub repository](https://github.com/aalexmmaldonado/reptar) or from [TestPyPI](https://test.pypi.org/project/reptar/).

```bash
git clone https://github.com/aalexmmaldonado/reptar
cd reptar
pip install .
```

# File types

Reptar supports three file types with a single interface: exdir, JSON, and npz.
JSON is a text file for storing key-value pairs with few dimensions (i.e., no large arrays).
NumPy's npz format is useful for arrays; however, no nesting is possible and loading data often requires postprocessing for 0D arrays (e.g., ``np.array('data')``).

Exdir is a simple, yet powerful open file format that mimics the [HDF5](https://www.hdfgroup.org/solutions/hdf5/) format with metadata and data stored in directories with YAML and npy files instead of a single binary file.
This provides several advantages such as mixing human-readable YAML and binary NumPy files, being easier for version control, and only loading requested portions of datasets into memory.
For more detailed information, please read this [*Front. Neuroinform.* article about exdir](https://doi.org/10.3389/fninf.2018.00016).

# Key-value pairs

All data is stored under a ``key``-``value`` pair within the reptar framework.
The ``key`` tells reptar where the data is stored and is conceptually related to standard file paths (without file extensions).
Nested data is specified by separating the nested keys with a ``/``.
For example, ``energy_pot``, ``md_run/geometry``, and ``entity_ids`` are all valid keys.
Note that ``gradients`` and ``/gradients`` would translate to the same value (``/`` species the "root" of the file).

# Workflow

## Storing data

We refer to a "reptar file" as any file that can be used with the ``reptar.File`` class.
Creating a reptar file starts by having a set of data files generated from some calculation.
Paths to these data files are passed into ``reptar.Creator.from_calc`` that extracts information using a ``reptar.parser`` class.
Information parsed from these files, ``parsed_info``, is then used to populate a ``reptar.File`` object.

Data can also be manually added by using ``File.put(key, data)`` where ``key`` is a string specifying where to store the data.

## Accessing data

Data can be added or retrieved using the same interface regardless of the underlying file format (e.g., exdir, JSON, and npz).
The only thing required is the respective ``key`` specifying where it is stored.
Then, ``File.get(key)`` can retrieve the data.

When working with JSON and npz files, ``File.save()`` must be explicitly called after any modification.

## Writing to other formats

Other packages often require data to be formatted in their own specific way.
Reptar provides ways to extract data from reptar files using ``File.get(key)`` and passing it into the desired ``reptar.writer`` function.
Reptar currently automates the creation of:

- [Atomic simulation environment (ASE) databases](https://wiki.fysik.dtu.dk/ase/tutorials/tut06_database/database.html),
- [Gaussian approximate potentials (GAP) extended XYZ files](https://libatoms.github.io/GAP/gap_fit.html#data),
- [Protein data bank (PDB) files](https://www.wwpdb.org/documentation/file-format),
- [Schnetpack databases](https://schnetpack.readthedocs.io/en/stable/tutorials/tutorial_01_preparing_data.html),
- [XYZ files](https://en.wikipedia.org/wiki/XYZ_file_format).

# License

Distributed under the MIT License. See `LICENSE` for more information.
