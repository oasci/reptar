<h1 align="center">Reptar</h1>

<h4 align="center">A tool for storing and analyzing computational chemistry data in a file-agnostic way</h4>

<p align="center">
    <a href="https://app.travis-ci.com/github/aalexmmaldonado/reptar" target="_blank">
        <img src="https://app.travis-ci.com/aalexmmaldonado/reptar.svg?branch=main" alt="Travis Build">
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
    <a href="https://lgtm.com/projects/g/aalexmmaldonado/reptar/context:python">
        <img src="https://img.shields.io/lgtm/grade/python/g/aalexmmaldonado/reptar.svg?logo=lgtm&logoWidth=18" alt="Language grade: Python">
    </a>
    <a href="https://github.com/aalexmmaldonado/reptar/blob/main/LICENSE" target="_blank">
        <img src="https://img.shields.io/github/license/aalexmmaldonado/reptar" alt="License">
    </a>
    <a href="https://github.com/aalexmmaldonado/reptar/" target="_blank">
        <img src="https://img.shields.io/github/repo-size/aalexmmaldonado/reptar" alt="GitHub repo size">
    </a>
</p>

<p align="center">
    <a href="#motivation">Motivation</a> •
    <a href="#about">About</a> •
    <a href="#license">License</a>
</p>

# Motivation

Computational chemistry is falling behind in providing open raw and processed data used to draw scientific conclusions.
Often it is the lack of time, expertise, and options that push researchers to overlook the importance of reproducible findings.
Projects such as [QCArchive](https://qcarchive.molssi.org/), [Materials Project](https://materialsproject.org/), [Pitt Quantum Repository](https://pqr.pitt.edu/), [ioChem-BD](https://www.iochem-bd.org/) and many others provide a rigid data framework for a specific purpose (e.g., quantum chemistry and material properties).
In other words, data that does not directly fit into their paradigm are incompatible (for good reason).

Alternatively, you could use file formats such as JSON, XML, YAML, npz, etc. for a specific project.
These are great options for customizable data storage with their own advantages and disadvantages.
However, you often must choose between (1) a standardized parser that might not support your workflow or (2) writing your own.

Reptar provides customizable parsers and data storage frameworks for whatever an individual project demands.
Data is stored in one of the supported file types and generalized routines are used to access and store data.
All data is stored in a key-value pair format where users can use predetermined definitions or include their own.
Regardless if you are running nudged elastic band calculations in VASP, free energy perturbation simulations in GROMACS, or gradient calculations in Psi4, you can store data easily with reptar by selecting a parser a specifying the desired file type.
The result is a user-specified data file streamlined for analysis in Python and optimized for archival on places such as [GitHub](https://github.com/) and [Zenodo](https://zenodo.org/).

# File types

Reptar supports three file types with a single interface: exdir, JSON, and npz.
JSON is a text file for storing key-value pairs with few dimensions (i.e., no large arrays).
NumPy's npz format is useful for arrays; however, no nesting is possible and loading data often requires postprocessing for 0D arrays (e.g., ``np.array('data')``).

Exdir is a simple, yet powerful open file format that mimics the [HDF5](https://www.hdfgroup.org/solutions/hdf5/) format with metadata and data stored in directories with YAML and npy files instead of a single binary file.
This provides several advantages such as mixing human-readable YAML and binary NumPy files, being easier for version control, and only loading requested portions of datasets into memory.
For more detailed information, please read this [*Front. Neuroinform.* article about exdir](https://doi.org/10.3389/fninf.2018.00016).

# Installation

You can install reptar using ``pip install reptar`` or install the latest version directly from the [GitHub repository](https://github.com/aalexmmaldonado/reptar).

    git clone https://github.com/aalexmmaldonado/reptar
    cd reptar
    pip install .

# License

Distributed under the MIT License. See `LICENSE` for more information.
