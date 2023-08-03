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

Reptar uses a driver/supervisor and worker workflow where the results are directly stored in a reptar :class:`~reptar.File`.
When running calculations, you first create a property (e.g., energy) array where all values are ``NaN``.
Each ``NaN`` value represents a calculation that still needs to run.
The driver then spawns workers with batches of calculations to run.
Once a worker is finished, we use the :meth:`~reptar.calculators.Data.save` method to store the results in case the job terminates early.

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






Examples
========


H2O energy+gradient (Psi4)
--------------------------

This provides a script of computing energy and gradient (`engrad`) calculations with Psi4 of water molecules.
We use :download:`this zarr file<./files/data/1h2o.zarr.zip>` (make sure to extract the file first).

.. tab-set::

    .. tab-item:: Script

        .. literalinclude:: ./files/scripts/1h2o-psi4-engrads.py
            :language: python
            :linenos:


    .. tab-item:: Output

        .. literalinclude:: ./files/scripts/1h2o-psi4-engrads.txt
            :language: text
            :linenos:

