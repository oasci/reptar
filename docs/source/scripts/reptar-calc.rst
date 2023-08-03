===============
``reptar-calc``
===============

This script provides a user-friendly interface to run calculations using reptar drivers and workers.
YAML files are used to specify how to run the calculations.

.. code-block:: text

    $ reptar-calc --help
    usage: reptar-calc [-h] [--ray_address [RAY_ADDRESS]] [--log_level [LOG_LEVEL]] [config_path]

    Run calculations using reptar

    positional arguments:
    config_path           Path to YAML configuration file

    options:
    -h, --help            show this help message and exit
    --ray_address [RAY_ADDRESS]
                            Desired ray address (will override config file)
    --log_level [LOG_LEVEL]
                            Desired logging level

Example YAML files
==================

xTB
---

.. dropdown:: Optimization

  .. code-block:: yaml

      # Information about the physical system.
      system:

        # Total charge.
        charge: 0

        # Multiplicity.
        multiplicity: 1

      # Information about the reptar file containing structures to optimize and where
      # to store the results.
      rfile:

        # File path to reptar file.
        path: data.zarr

        # Key to the source of structures to optimize.
        source_key: ./

        # Label for the atomic numbers data in ``source``.
        Z_label: atomic_numbers

        # Label for the Cartesian coordinate data in ``source``.
        R_label: geometry

        # Key to store the optimization results. We recommend this be different than
        # ``source_key``.
        dest_key: opt_gfn2

        # Label to store if an optimization converged in ``destination``.
        conv_opt_label: conv_opt

        # Label to store the Cartesian coordinates of the last optimization step in
        # ``destination``.
        R_opt_label: geometry

        # Label to store the electronic energy of the last optimization step in
        # ``destination``.
        E_opt_label: energy_ele_gfn2

      # Driver for the desired calculation.
      driver:

        # Import string for the desired calculation driver.
        class: reptar.calculators.drivers.DriverOpt

        # Keyword arguments for the specified driver in ``class``.
        # Please refer to the specific driver's API in the documentation.
        # https://www.aalexmmaldonado.com/reptar/main/index.html
        kwargs:

          use_ray: true

          ray_address: auto

          n_workers: 2

          n_cpus_per_worker: 4

          chunk_size: 10

          start_slice: null

          end_slice: null

      # Worker to be used with the driver.
      worker:

        # Import string for the desired worker.
        function: reptar.calculators.xtb_workers.xtb_opt

        # Blocks specified in xcontrol that will be used to generate ``input_lines``.
        # Each key and value should be specified in the same format as xTB desires.
        # https://github.com/grimme-lab/xtb/blob/main/man/xcontrol.7.adoc
        blocks:

          opt:
            
            optlevel: verytight

        # Special handling of geometry constraints.
        # https://xtb-docs.readthedocs.io/en/latest/xcontrol.html#constraining-potentials
        constrain: null
        # Alternatively, below is an example where we add a constraint to xTB. These formats
        # create the nested lists format.
        # constrain:
        #   - 
        #     - distance
        #     - 
        #       - 0
        #       - 3
        #       - 1.4
        #   - 
        #     - angle
        #     - 
        #       - 4
        #       - 5
        #       - 6
        #       - auto

        # Keyword arguments for the ``xtb_opt`` worker. ``input_lines`` is handled separately.
        # These are kept separate as they must be specified as a command line argument.
        # Please refer to the specific workers's API in the documentation.
        # https://www.aalexmmaldonado.com/reptar/main/index.html
        kwargs:

          acc: 0.01
          
          n_cores: 4

          xtb_path: xtb

          log_dir: /home/alex/repos/resiparm-examples/gfp-cro/5-reptar-xtb-opt/logs

Psi4
----

.. dropdown:: Energy and gradient

  .. code-block:: yaml

    # Information about the physical system.
    system:

      # Total charge.
      charge: 0

      # Multiplicity.
      multiplicity: 1

    # Computations to complete in the order they are provided.
    tasks:
      - E
      - G

    # Information about the reptar file containing structures to optimize and where
    # to store the results.
    rfile:

      # File path to reptar file.
      path: ./tmp/calculators/1h2o-psi4-script.zarr

      # Information for data sources
      source:

        # Key to group to retrieve data.
        key: /1h2o

        # Only used to get data.
        labels:
          Z: atomic_numbers
          R: geometry
      
      # Information for where to store calculated data.
      # This location can be the same as ``source``
      destination:

        # Key to group to store data.
        key: 1h2o

        # Used to put data.
        labels:
          Z: atomic_numbers
          R: geometry
          E: energy_ele_mp2.def2tzvp
          G: grads_mp2.def2tzvp

    # Driver for the desired calculation.
    driver:

      # Keyword arguments for the specified driver in ``class``.
      # Please refer to the specific driver's API in the documentation.
      # https://www.aalexmmaldonado.com/reptar/main/index.html
      kwargs:

        use_ray: false
        n_workers: 1
        n_cpus_per_worker: 2
        chunk_size: 1
        ray_address: auto
      
      start_slice: null

      end_slice: null

    # Worker to be used with the driver.
    worker:

      # Import string for the desired worker.
      path: reptar.calculators.psi4_workers.psi4_worker

      # Keyword arguments for the worker.
      # These are kept separate as they must be specified as a command line argument.
      # Please refer to the specific workers's API in the documentation.
      # https://www.aalexmmaldonado.com/reptar/main/index.html
      kwargs:

        charge: 0

        mult: 1

        method: mp2

        options:

          basis: def2-TZVPPD
          df_basis_scf: def2-universal-jkfit
          df_basis_mp2: def2-tzvppd-ri
          reference: rhf
          scf_type: df
          e_convergence: 10
          d_convergence: 10
          mp2_type: df
          qc_module: dfmp2
          freeze_core: true

        threads: 2

        mem: 1 GB



