====================
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
