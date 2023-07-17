=======
Scripts
=======

Reptar has a few helper scripts to make life easier.

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

.. dropdown:: Example YAML configuration files
    :animate: fade-in-slide-down

    .. tab-set::

        .. tab-item:: xTB optimization

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

                    n_cpus_per_worker: 3

                    chunk_size: 5

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
                    
                    n_cores: 3

                    xtb_path: xtb

                    work_dir: ~/repos/resiparm-examples/gfp-cro/5-reptar-xtb-opt/logs

          
        .. tab-item:: Psi4 optimization

          .. code-block:: yaml

            # Information about system.
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
              dest_key: opt_bp86

              # Label to store if an optimization converged in ``destination``.
              conv_opt_label: conv_opt

              # Label to store the Cartesian coordinates of the last optimization step in
              # ``destination``.
              R_opt_label: geometry

              # Label to store the electronic energy of the last optimization step in
              # ``destination``.
              E_opt_label: energy_ele_bp86.def2svp

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

                n_cpus_per_worker: 3

                chunk_size: 5

                start_slice: null

                end_slice: null

            # Worker to be used with the driver.
            worker:

              # Import string for the desired worker.
              function: reptar.calculators.psi4_workers.psi4_opt

              # Keyword arguments for the ``psi4_opt`` worker. ``input_lines`` is handled separately.
              # These are kept separate as they must be specified as a command line argument.
              # Please refer to the specific workers's API in the documentation.
              # https://www.aalexmmaldonado.com/reptar/main/index.html
              kwargs:

                charge: 0

                mult: 1

                method: bp86

                options:

                  basis: def2-svp

                  e_convergence: 6

                  d_convergence: 6

                  geom_maxiter: 200

                  g_convergence: gau



``reptar-write-xyz``
====================

We love to work in XYZ files; however, text files can quickly grow in size.
Binary files, like ``npy``, store the same data in less space.
Most visualization packages require XYZ files, so this script quickly looks up the relevant data in a group and writes an XYZ file.

.. code-block:: text

    $ reptar-write-xyz --help
    usage: reptar-write-xyz [-h] [--group_key [GROUP_KEY]] [--comment_key [COMMENT_KEY]] [--save_dir [SAVE_DIR]]
                        [group_path]

    Write XYZ file from a reptar supported file type.

    positional arguments:
    group_path            Path to group path. This can be to a file or exdir/zarr group.

    options:
    -h, --help            show this help message and exit
    --group_key [GROUP_KEY]
        Manually specify group key for file.
    --comment_key [COMMENT_KEY]
        Label of data to include as xyz comment.
    --save_dir [SAVE_DIR]
        Directory to save XYZ file.

.. admonition:: Example

    As long as we have ``atomic_numbers`` and ``geometry`` data we can run quickly get an XYZ file.

    .. code-block:: bash

        $ ls
        file.exdir
        $ reptar-write-xyz ./file.exdir/parent_key/group_key
        reptar v0.0.3.post41
        Parsing args
        Collecting data
        Writing XYZ file
        Done!
        $ ls
        file.exdir  data.xyz
    
    If we had electronic energies stored in ``energy_ele`` under ``group_key`` we could include them as comments using ``--comment_key energy_ele``.
    Note that we assume the energies are stored in Hartrees and will be automatically converted to kcal/mol.