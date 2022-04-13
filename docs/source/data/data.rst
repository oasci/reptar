.. role:: raw-html(raw)
   :format: html

====
Data
====





identification
==============

``r_centered``
--------------

*Type:* bool
    If structures (likely sampled) were translated such that the center of mass is at the origin.

``md5``
-------

*Type:* string
    Unique string determined by all group attributes and datasets.

``md5_data``
------------

*Type:* string
    Unique string determined by all group datasets (no attributes are included).
    This is a data-focused MD5 that will not be affected by adding or changing less-important attributes.
    Use this for breadcrumb purposes.

``r_prov_ids``
--------------

*Type:* dict
    Species an ID (``int``) to uniquely identifying labels for each structure if it originated from another reptar file.
    Labels are usually ``md5_data``.
    For example, ``{0: '6038e101da7fc0085978741832ebc7ad', 1: 'eeaf93dec698de3ecb55e9292bd9dfcb'}``.

``r_prov_specs``
----------------

*Type:* array
    Specifies ``structure_prov_id``, the ``structure index`` and original ``entity_id`` of the reptar prov.
    Take the following array for example.

    .. code-block:: python

        [[0, 483,   0, 12],
         [0,  89, 194,  2],
         [1,   0,   0,  1]]

    This indicates that we should have three structures in our ``Group``.
    Structure one and two come from the ``Group`` that has a ``structure_prov_id`` of ``0`` whereas the last one is from a different ``Group`` (``structure_prov_id`` of ``1``).
    The first structure in this ``Group`` comes from the 484th (index of ``483``) structure in the source ``Group``.
    But this structure does not contain every atom from the original ``Group`` structure; it only contains the atoms from ``0`` and ``12`` ``entity_id``.





system_info
===========

``atomic_numbers``
------------------

*Type:* array
    Number of protons in the nucleus of each atom.

``geometry``
------------

*Type:* array
    Numerical coordinates specifying the three-dimensional positions of each atom.

``charge``
----------

*Type:* int
    Number of protons in the nucleus of each atom.

``comp_ids``
-----------------

*Type:* list
    Relates ``entity_id`` to a fragment label for chemical components or species.
    Labels could be WAT or h2o for water, MeOH for methanol, bz for benzene, etc.
    There are no standardized labels for species.
    The index of the label is the respective ``entity_id``.
    For example, a water and methanol molecule could be ``['h2o', 'meoh']``.

``comp_ids_num``
-----------------

*Type:* dict
    Specifies the number of component IDs as an attribute.

``entity_ids``
--------------

*Type:* array
    A uniquely identifying label specifying what atoms belong to which entities.
    Entities can be a related set of atoms, molecules, or functional group.

    For example, a water and methanol molecule could be ``[0, 0, 0, 1, 1, 1, 1, 1, 1]``.

``mult``
----------------

*Type:* int
    System multiplicity.

``n_ele``
-----------------

*Type:* int
    Total number of electrons in the system.





runtime_info
============

``basis_aux``
-------------------

*Type:* string
    TODO

``calc_driver``
---------------

*Type:* string
    The purpose of the calculation such as calculate energies, frequencies, gradients, or perform optimizations.
    Possible values include: ``energy``, ``frequency``, ``gradient``, ``optimization``.

``approx_cosx``
----------------------

*Type:* boolean
    TODO

``ele_frozen``
--------------------

*Type:* int
    TODO

``hamiltonian``
---------------

*Type:* string
    TODO

``hf_type``
-----------

*Type:* string
    TODO

``solv_implicit_name``
-------------------------

*Type:* string
    TODO

``t_step``
----------------------

*Type:* float
    TODO MD time step in fs

``md_duration``
---------------

*Type:* float
    TODO in ps

``md_steps_dump``
-----------------

*Type:* int
    TODO number of steps to dump information

``prov``
--------------

*Type:* string
    TODO

``prov_version``
----------------------

*Type:* string
    TODO

``basis_n_func``
-----------------------

*Type:* int
    Number of basis functions.

``rij_approximation``
---------------------

*Type:* boolean
    TODO

``approx_rijk``
----------------------

*Type:* boolean
    TODO

``scf_grid_level``
------------------

*Type:* string
    TODO

``scf_grid_level_final``
------------------------

*Type:* string
    TODO

``success``
-----------

*Type:* boolean
    TODO

``temp_thermostat``
-------------------

*Type:* float
    TODO





outputs
=======

``grads``
---------

*Type:* array
    Atomic gradients.

``forces``
----------

*Type:* array
    TODO

``energy_ccsd``
---------------------

*Type:* array or float
    TODO

``energy_ccsd(t)``
------------------------

*Type:* array or float
    TODO

``dipole_moment``
-----------------

*Type:* array
    Dipole moment vector in Debyes.

``energy_ele``
---------------------

*Type:* array or float
    Total electronic energy including exchange and correlation.


``correc_enthalpy``
-------------------------

*Type:* array or float
    Enthalpic energy correction to the electronic energy.

``correc_entropy``
------------------------

*Type:* array or float
    TODO

``conv_targ_geo_energy``
----------------------------

*Type:* array
    TODO

``conv_val_geo_energy``
---------------------------

*Type:* array
    TODO

``conv_targ_geo_grad``
---------------------------

*Type:* array
    TODO

``conv_val_geo_grad``
--------------------------

*Type:* array
    TODO

``conv_targ_geo_grad_max``
---------------------------

*Type:* array
    TODO

``conv_val_geo_grad_max``
--------------------------

*Type:* array
    TODO

``charges_lowdin``
-------------------------

*Type:* array
    Estimated partial atomic charges using the L\ :raw-html:`&#xf6;`\ wdin atomic charge scheme.

``charges_mulliken``
---------------------------

*Type:* array
    Estimated partial atomic charges using a Mulliken population analysis.

``energy_mos``
------------------------------

*Type:* array
    Energies of molecular orbitals.
    The array is of length 1 for restricted calculations, but length 2 for unrestricted calculations.

``energy_correl_mp2``
--------------------------

*Type:* array or float
    Correlation energy from MP2.

``normal_modes``
----------------

*Type:* array
    Normalized, mass-weighted normal modes, q.

``energy_nuc_repul``
----------------------------

*Type:* array or float
    TODO

``freq_vib``
---------

*Type:* array
    TODO

``energy_pot``
--------------------

*Type:* array or float
    TODO

``energy_scf_one_ele``
---------------------------

*Type:* array or float
    TODO

``energy_scf``
--------------------

*Type:* array or float
    TODO

``energy_scf_two_ele``
---------------------------

*Type:* array or float
    TODO

``energy_scf_xc``
-----------------

*Type:* array or float
    TODO

``diag_t1``
-----------------

*Type:* array or float
    TODO

``temp_thermochem``
-------------------

*Type:* float
    TODO

``correc_thermal``
------------------------------

*Type:* array or float
    TODO

``spin_sqrd_uhf_ideal``
----------------------------------------

*Type:* array or float
    uhf_ideal_average_total_spin_squared

``spin_sqrd_uhf_calc``
---------------------------------------------

*Type:* array or float
    uhf_calculated_average_total_spin_squared

``correc_zpe``
-------------------------------------

*Type:* array or float
    TODO

.. toctree::
    :hidden:
    :maxdepth: 2
   
    xTB data <xtb-data>