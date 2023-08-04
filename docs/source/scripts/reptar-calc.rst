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

  .. literalinclude:: ../files/scripts/reptar-calc-xtb-opt.yml
    :language: yaml
    :linenos:

Psi4
----

.. dropdown:: Energy and gradient

  .. literalinclude:: ../files/scripts/reptar-calc-psi4-engrads.yml
    :language: yaml
    :linenos:



