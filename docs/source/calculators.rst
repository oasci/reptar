===========
Calculators
===========

Reptar provides several driver and workers for calculations on data stored in supported file formats.
The :class:`~reptar.calculators.Driver` class manages all calculations and workers for reptar.

- `Psi4 <https://psicode.org/psi4manual/master/index.html>`__: quantum chemical methods such as DFT and wave function methods.
    - :func:`~reptar.calculators.psi4_workers.psi4_worker`
- `xtb <https://xtb-docs.readthedocs.io/en/latest/contents.html>`__ and `xtb-python <https://xtb-python.readthedocs.io/en/latest/>`__: a semiempirical quantum mechanics method.
    - :func:`~reptar.calculators.xtb_workers.xtb_worker`, :func:`~reptar.calculators.xtb_workers.xtb_python_worker`

We use `ray <https://docs.ray.io/en/latest/ray-overview/installation.html>`__ to parallelize our calculations across one or multiple nodes.





Drivers and workers
===================

Reptar uses a driver (supervisor) and worker workflow where the results are stored in a reptar :class:`~reptar.File`.
When running calculations, an array is initialized for a desired property (e.g., energy) where all values are ``NaN``.
``NaN`` values represents a calculation that still needs to run.
The driver then spawns workers with batches of calculations to run.
Results are stored as soon as a worker finishes.

.. mermaid::

    flowchart LR
        start((Start)) --> props[/Property<br/>arrays/]
        props -- store --> file[(File)]
        file --> idxs[/NaN<br/>indices/]
        idxs --> driver[Driver]
        subgraph workers [Workers]
        worker1[Worker 1]
        worker2[Worker 2]
        worker3[. . .]
        end
        driver --> batchInfo[/Batch<br/>info/]
        batchInfo --> workers
        workers --> batchResults[/Batch<br/>results/]
        batchResults -- update --> file


Data
====

In order to track and manage calculations, we store relevant data using :class:`~reptar.calculators.Data`.
It contains all supported properties such as atomic positions, energies, gradients, and grid data.

Most calculations involve taking :attr:`~reptar.calculators.Data.Z` and :attr:`~reptar.calculators.Data.R` from a ``source`` group and then storing calculations in a destination (``dest``) group.
:attr:`~reptar.calculators.Data.prepare_tasks` can assist you with copying over data to the destination and initializing arrays for computed data.

Tasks
=====

Each worker has a number of ``tasks`` they can compute based on the external code and reptar support.


``E``
-----

Computes the total energy of a structure with respect to the provided worker and options.

**Required data:** :attr:`~reptar.calculators.Data.E`

``G``
-----

TODO

**Required data:** :attr:`~reptar.calculators.Data.E` and :attr:`~reptar.calculators.Data.G`


``opt``
-------

TODO

**Required data:** :attr:`~reptar.calculators.Data.E`, :attr:`~reptar.calculators.Data.conv_opt`, and :attr:`~reptar.calculators.Data.R_opt`

``cube``
--------

TODO

**Required data:** :attr:`~reptar.calculators.Data.cube_R` and :attr:`~reptar.calculators.Data.cube_V`


Examples
========


H2O energy+gradient (Psi4)
--------------------------

This provides a script of computing energy and gradient (`engrad`) calculations with Psi4 of water molecules.
We use :download:`this zarr file<./files/data/1h2o.zarr.zip>` (make sure to extract the file first).

.. tab-set::

    .. tab-item:: Script

        .. literalinclude:: ./files/scripts/1h2o_psi4_engrads.py
            :language: python
            :linenos:


    .. tab-item:: Output

        .. literalinclude:: ./files/scripts/1h2o-psi4-engrads.txt
            :language: text
            :linenos:

