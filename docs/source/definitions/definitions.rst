===========
Definitions
===========

One challenging aspect of repositories is communicating the meaning of each label we give our data.
Reptar tries to suggest common data types encountered in the computational chemistry community.
Labels and their respective descriptions are in YAML files that could be stored directly in the data.
We consider three types of data:

- **identification:** unique labels that describe the within.
- **system_info:** specifies the atomistic system.
- **runtime_info:** results and information from calculations.

We further separate definitions into the following categories.

.. toctree::
    :maxdepth: 1

    base
    md
    molprop
    pes
    qc
    sampling
    solv
    xtb
