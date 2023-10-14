====
base
====

identification
--------------

.. glossary::

    md5
        Unique string determined by all group data.

        **Type:** :obj:`str`

    md5_arrays
        Unique string determined by all group arrays (i.e., data not stored in arrays are not included). This MD5 hash will not change as much when non-array data is added or removed.

        **Type:** :obj:`str`

    md5_structures
        Unique string determined by all group datasets (no attributes are included). This MD5 is only dependent on ``atomic_numbers`` and ``geometry``. Use this for breadcrumb purposes.

        **Type:** :obj:`str`

    readme
        Provides any information about the group.

        **Type:** :obj:`str`

system_info
-----------

.. glossary::

    atomic_numbers
        Number of protons in the nucleus of each atom.

        **Type:** :obj:`int`

    charge
        Number of protons in the nucleus of each atom.

        **Type:** :obj:`int`

    geometry
        Numerical coordinates specifying the three-dimensional positions of each atom.

        **Type:** :obj:`float`

    comp_ids
        Relates ``entity_id`` to a fragment label for chemical components or species. Labels could be WAT or h2o for water, MeOH for methanol, bz for benzene, etc. There are no standardized labels for species. The index of the label is the respective ``entity_id``. For example, a water and methanol molecule could be ``['h2o', 'meoh']``.

        **Type:** :obj:`str`

    comp_ids_num
        Specifies the number of component IDs as an attribute. For example, ``{'H2O': 1, 'MEOH': 120}``

        **Type:** :obj:`dict`

    entity_ids
        A uniquely identifying integer specifying what atoms belong to which entities. Entities can be a related set of atoms, molecules, or functional group. For example, a water and methanol molecule could be ``[0, 0, 0, 1, 1, 1, 1, 1, 1]``.

        **Type:** :obj:`int`

    mult
        System multiplicity.

        **Type:** :obj:`int`

    n_ele
        Total number of electrons in the system.

        **Type:** :obj:`int`

    periodic
        If periodic boundary conditions are enforced in the x, y, or z direction, respectively.

        **Type:** :obj:`bool`

    periodic_cell
        The three cell vectors defining the periodic cell.

        **Type:** :obj:`float`

runtime_info
------------

.. glossary::

    calc_driver
        The purpose of the calculation such as calculate energies, frequencies, gradients, or perform optimizations. Possible values include: ``energy``, ``frequency``, ``gradient``, ``optimization``, ``molecular_dynamics``.

        **Type:** :obj:`str`

    prov
        TODO

        **Type:** :obj:`str`

    prov_version
        TODO

        **Type:** :obj:`str`

    success
        TODO

        **Type:** :obj:`bool`
