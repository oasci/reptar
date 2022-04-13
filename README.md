<h1 align="center">reptar</h1>

<h4 align="center"><b>Rep</b>roducible, <b>t</b>ailorable <b>ar</b>chive</h4>

<p align="center">
    <a href="https://github.com/aalexmmaldonado/reptar/blob/main/LICENSE" target="_blank">
        <img src="https://img.shields.io/github/license/aalexmmaldonado/reptar" alt="License">
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
In other words, data that does not directly fit into their paradigm are incompatible&mdash;and for good reason.

Reptar is different in that it is a customizable data storage framework for whatever an individual project demands.
Data is converted into reptar files where users can use predetermined key value pair definitions or define their own.
Regardless if you are running nudged elastic band calculations in VASP, free energy perturbation simulations in GROMACS, or gradient calculations in Psi4, you can create reptar files by using a parser that retrieves data and stores them in the appropriate fashion.
The result is a user-specified data file streamlined for analysis in Python and optimized for archival on places such as GitHub and [Zenodo](https://zenodo.org/).

# About

Reptar is essentially a collection of tools and routines for using the [exdir](https://exdir.readthedocs.io/en/latest/) framework with typical computational chemistry workflows.
Exdir is a simple, yet powerful open file format that mimics the [HDF5](https://www.hdfgroup.org/solutions/hdf5/) format with metadata and data stored in directories with YAML and npy files instead of a single binary file.
This provides several advantages such as mixing human-readable YAML and binary NumPy files, being easier for version control, and only loading requested portions of datasets into memory.
For more detailed information, please read this [*Front. Neuroinform* article about exdir](https://doi.org/10.3389/fninf.2018.00016).

# TODO

- [ ] Write scripts that create the typical files (i.e., xyz and pbd) from groups.
- [ ] Figure out way to document user-specified data in reptars so we can customize each file.

# License

Distributed under the MIT License. See `LICENSE` for more information.
