====
File
====

The heart of reptar is the :class:`~reptar.File` class which allows you to interact with your data in a storage-agnostic way.

Storage types
=============

Every file type comes with its own list of advantages and disadvantages leading to best use cases.
We let the user decide what file format to use and have reptar provide a consistent way of interacting with the data.
Reptar currently supports the following file types:

- `exdir <https://exdir.readthedocs.io/en/latest/>`__,
- `zarr <https://zarr.dev/>`__,
- JSON,
- npz.

Any data stored in this format can be interacted with through a :class:`~reptar.File` object.

Concepts
========

All data is stored in key-value pairs.
The key provides a human-readable label to communicate the meaning of the value.
For example, ``energy_pot`` and ``geometry`` refers to potential energy and Cartesian coordinates, respectively.
More information about key definitions can be found :ref:`here<definitions>`.

The key also tells reptar where the data is stored.
Many file types (e.g., exdir, zarr, and JSON) support nested objects which allow us to organize our data into groups.
A group is essentially a collection of related data (e.g., a directory) that could also contain other groups.
We specify keys for nested groups the same way we do directory paths where the file itself is considered the root: ``group1``, ``/group1``, and ``group1/nested_group`` are all valid group keys.

For example, suppose we have a file that contains two molecular dynamics (MD) runs of a box of water.
We can visualize a group tree where we have two groups (``run_1`` and ``run_2``) that contain the ``geometry`` and ``energy_pot`` from both MD simulations.
If we wanted to get the potential energies of both runs we would use ``/run_1/energy_pot`` and ``/run_2/energy_pot``.
Note that the ``/`` prefix is not required (``run/energy_pot`` works as well).

.. code-block:: text

    .
    └── water-box-md
        ├── run_1
        │   ├── energy_pot
        │   ├── geometry
        ├── run_2
        │   ├── energy_pot
        │   └── geometry

Retrieving data
===============

Reptar provides a single method, :meth:`~reptar.File.get`, to retrieve data from any of the supported file types.
This can return any nested data stored in a loaded file including other groups.

You can also get a list of all non-group keys using :meth:`~reptar.File.get_keys`.


Adding data
===========

Reptar provides a single method, :meth:`~reptar.File.put`, to store data in any of the supported file types.

Note that exdir groups must be initialized first with :meth:`~reptar.File.create_group` before adding data.

Saving
======

JSON and npz files need to be explicitly saved using :meth:`~reptar.File.save`.

File conversion
===============

A group from any file type can be directly converted to a dictionary using :meth:`~reptar.File.as_dict`.

Dictionaries and can also be used to populate a :class:`~reptar.File` object with any specified file type.

