========================
``reptar-geometry-scan``
========================

We often want to see how the potential energy surface changes with respect to some internal coordinates.
This script lets you specific ``angle`` and ``dihedral`` scans by specifying atom indices and angle ranges.
You can then follow up with a constrained geometry optimization.

.. code-block:: text

    $ reptar-geometry-scan --help
    usage: reptar-geometry-scan [-h] [--ray_address [RAY_ADDRESS]] [--log_level [LOG_LEVEL]]
                            [config_path]

    Scan internal coordinates with an optional geometry optimization

    positional arguments:
    config_path           Path to YAML configuration file

    options:
    -h, --help            show this help message and exit
    --ray_address [RAY_ADDRESS]
                            Desired ray address (will override config file)
    --log_level [LOG_LEVEL]
                            Desired logging level

.. attention::

    Only xTB optimizations have been validated.

Example YAML files
==================

GFP chromophore
---------------

TODO: 

The below YAML config results in the following :download:`zarr file<../files/data/cro-dihedral-scan-script.zarr.zip>`.

.. dropdown:: YAML

  .. literalinclude:: ../files/scripts/config-gfp-cro-scan.yaml
    :language: yaml
    :linenos:
