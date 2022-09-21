===========
Calculators
===========

Reptar provides a framework for computing common properties and electronic energies and gradients.
`Ray <https://docs.ray.io/en/latest/index.html>`__ is used to parallelize computations when possible.

The additional dependencies are:

* `ray <https://docs.ray.io/en/latest/ray-overview/installation.html>`__ (required)
* `psi4 <https://psicode.org/installs/v15/>`__
* `xtb-python <https://github.com/grimme-lab/xtb-python#installation>`__

Energy and gradients
====================

Reptar provides a framework with a ray driver and premade ray task/worker functions.

.. autoclass:: reptar.calculators.drivers.driverENGRAD
    :noindex:

``worker_func`` calls the desired program and computes energies and gradients.

.. note::

    ``worker_func`` should not already have the ``ray.remote`` decorator.
    This is automatically included when :class:`~reptar.calculators.drivers.driverENGRAD` is initialized.

Reptar currently has a worker for `Psi4 <https://psicode.org/psi4manual/master/index.html>`__, :func:`~reptar.calculators.psi4_workers.psi4_engrad`, and `xTB <https://xtb-docs.readthedocs.io/en/latest/contents.html>`__, :func:`~reptar.calculators.xtb_workers.xtb_engrad`. 

.. autofunction:: reptar.calculators.psi4_workers.psi4_engrad
    :noindex:

.. _psi4-engrad-example:

Psi4 example
------------

The following example demonstrates a parallelized workflow for computing MP2/cc-pVTZ energies and gradients for twenty (H\ :sub:`2`\ O)\ :sub:`3` structures.

.. code-block:: python

    import numpy as np
    import ray
    from reptar import File
    from reptar.calculators.drivers import driverENGRAD
    from reptar.calculators.psi4_workers import psi4_engrad

    # Initialize ray and worker CPUs.
    n_cpus = 8  # Total number of CPUs to use.
    n_cpus_worker = 2  # Will use four parallel workers.
    ray.init(num_cpus=n_cpus)

    # Get exdir data.
    exdir_path = '3h2o.exdir'
    f = File(exdir_path, mode='a')
    Z = f.get('atomic_numbers')
    R = f.get('geometry')

    # Create energy and gradient arrays.
    E = np.empty(R.shape[0])
    E[:] = np.nan
    G = np.empty(R.shape)
    G[:] = np.nan

    # Specify Psi4 properties
    charge = 0
    mult = 1
    psi4_method = 'mp2'
    psi4_options = {
        'basis': 'cc-pVTZ',
        'reference': 'rhf',
        'e_convergence': '1e-10',
        'd_convergence': '1e-10',
        'freeze_core': 'false'
    }
    psi4_threads = n_cpus_worker
    psi4_mem = '2 GB'
    psi4_args = (charge, mult, psi4_method, psi4_options, psi4_threads, psi4_mem)

    # Creating engrad driver.
    max_calcs = 20  # Manually trims R to the first 20 structures (can remove).
    engrads = driverENGRAD(
        Z, R, E, G, psi4_engrad, psi4_args, n_cpus,
        n_cpus_worker=n_cpus_worker, end_slice=max_calcs
    )

    # Run parallel Psi4 calculations
    E, G = engrads.run()

    # Add data to exdir file.
    f.put('energy_ele_mp2.ccpvtz_psi4', E)
    f.put('grads_mp2.ccpvtz_psi4', G)

