======================
``reptar-xyz-to-file``
======================

TODO:

.. code-block:: text

    $ reptar-xyz-to-file --help
    usage: reptar-xyz-to-file [-h] [--save_path [SAVE_PATH]] [--group_key [GROUP_KEY]] [--overwrite] [--log_level [LOG_LEVEL]]
                              [xyz_path]

    Make a reptar file from XYZ

    positional arguments:
    xyz_path              Path to XYZ file (default: None)

    options:
    -h, --help            show this help message and exit
    --save_path [SAVE_PATH]
                            File name (default: data.zarr)
    --group_key [GROUP_KEY]
                            Group key to store in file (default: ./)
    --overwrite           If file exists, overwrite all data (default: False)
    --log_level [LOG_LEVEL]
                            Desired logging level (default: info)
