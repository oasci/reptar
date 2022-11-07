=======
Scripts
=======

Reptar has a few helper scripts to make life easier.

``reptar-write-xyz``
====================

We love to work in XYZ files; however, text files can quickly grow in size.
Binary files, like ``npy``, store the same data in less space.
Most visualization packages require XYZ files, so this script quickly looks 

.. code-block:: bash

    $ reptar-write-xyz --help
    usage: reptar-write-xyz [-h] [--comment_key [COMMENT_KEY]] [--save_dir [SAVE_DIR]] [group_key]

    Write XYZ file from exdir file.

    positional arguments:
    group_key             Path to group key. Will take absolute path and split based on `.exdir`.

    options:
    -h, --help            show this help message and exit
    --comment_key [COMMENT_KEY]
                            Label of data to include as xyz comment. If "energy" or "grad" is in the key, then we convert to
                            kcal/mol.
    --save_dir [SAVE_DIR]
                            Directory to save XYZ file.

.. admonition:: Example

    As long as we have ``atomic_numbers`` and ``geometry`` data we can run quickly get an XYZ file.

    .. code-block:: bash

        $ ls
        file.exdir
        $ reptar-write-xyz ./file.exdir/parent_key/group_key
        reptar v0.0.2+73.g88e53e0
        Parsing args
        Collecting data
        Writing XYZ file
        Done!
        $ ls
        file.exdir  data.xyz
    
    If we had electronic energies stored in ``energy_ele`` under ``group_key`` we could include them as comments using ``--comment_key energy_ele``.
    Note that we assume the energies are stored in Hartrees and will be automatically converted to kcal/mol.