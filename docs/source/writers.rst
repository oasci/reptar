=======
Writers
=======

Each package often defines their own framework to specify and use data.
This makes transitioning to different packages challenging as one has to write scripts to convert their data into an acceptable format.
An updating game is then played every time new data is needed: update the source then rerun your conversion scripts.
Reptar provides these scripts for you along with comprehensive unit tests to ensure data is converted correctly.

ASE
===

The atomic simulation environment (ASE) provides `its own database format <https://wiki.fysik.dtu.dk/ase/tutorials/tut06_database/database.html>`_ for atomic data.
You can use :func:`~reptar.writers.write_ase_db` to write these databases as shown in the example below.

.. code-block:: python

    from reptar import File
    from reptar.writers import write_ase_db

    rfile_path = '1h2o_120meoh_md.exdir'
    db_path = '1h2o_120meoh_md_ase.db'

    rfile = File(rfile_path, mode='r')
    Z = rfile.get('prod_1/atomic_numbers')
    R = rfile.get('prod_1/geometry')  # Angstroms
    E = rfile.get('prod_1/energy_pot')  # Hartree
    E *= 27.21138602  # eV

    db = write_ase_db(db_path, Z, R, energy=E)

GAP
===

The Gaussian approximation potentials (GAP) code requires data in a specific extended XYZ format.
You can use :func:`~reptar.writers.write_xyz_gap` to write these files as shown in the example below.

.. code-block:: python

    from reptar import File
    from reptar.writers import write_xyz_gap

    rfile_path = '1h2o_120meoh_md.exdir'
    xyz_path = '1h2o_120meoh_md_gap.xyz'

    rfile = File(rfile_path, mode='r')
    Z = rfile.get('prod_1/atomic_numbers')
    R = rfile.get('prod_1/geometry')
    E = rfile.get('prod_1/energy_pot')  # Hartree
    E *= 27.21138602  # eV
    lattice = np.array(
        [[200.0, 0.0, 0.0],
         [0.0, 200.0, 0.0],
         [0.0, 0.0, 200.0]]
    )
    
    write_xyz_gap(xyz_path, lattice, Z, R, E)

GDML
====

The gradient-domain machine learning (GDML) packages (`sGDML <http://quantum-machine.org/gdml/>`_ and `mbGDML <https://keithgroup.github.io/mbGDML/>`_) require datasets to be in a npz format.
Since reptar innately supports npz files, you just need to use the ``from_dict`` option in :class:`~reptar.File`.
An example script is shown below.

.. note::

    When retrieving data from reptar files it is important to explicitly open in them ``r`` mode (to avoid accidentally overwriting data).
    Exdir uses `memory maps <https://numpy.org/doc/stable/reference/generated/numpy.memmap.html>`_ to avoid automatically loading full arrays into memory.
    You cannot make any changes to these data in `r` mode.
    Thus, we automatically convert all memmaps to arrays by default to avoid having to create array copies every time.
    You can request memmaps by using ``rfile.get(key, as_memmap=True)``.

.. code-block:: python

    import numpy as np
    from reptar import File

    hartree2kcalmol = 627.5094737775374055927342256  # Psi4 v1.5

    # Loading Exdir file.
    exdir_path = '140h2o-xtb.md-samples.exdir'
    rfile = File(exdir_path, mode='r')
    group_key = '3h2o'

    # Getting GDML-relevant data.
    Z = rfile.get(f'{group_key}/atomic_numbers')
    R = rfile.get(f'{group_key}/geometry')
    E = np.array(rfile.get(f'{group_key}/energy_ele_mp2.def2tzvp_orca'))   # Hartree
    G = np.array(rfile.get(f'{group_key}/grads_mp2.def2tzvp_orca'))  # Hartree/A
    entity_ids = rfile.get(f'{group_key}/entity_ids')
    comp_ids = rfile.get(f'{group_key}/comp_ids')
    r_prov_ids = rfile.get(f'{group_key}/r_prov_ids')
    r_prov_specs = rfile.get(f'{group_key}/r_prov_specs')

    # Unit conversions.
    E *= hartree2kcalmol  # kcal/mol
    G *= hartree2kcalmol  # kcal/(mol A)
    F = -G

    # Manual data.
    r_unit = 'Angstrom'
    e_unit = 'kcal/mol'
    theory = 'mp2/def2tzvp'

    # Specifying npz path and name.
    npz_path = '140h2o-3h2o.sampled-mp2.def2tzvp.npz'
    npz_name = npz_path[:-4]

    # Creating the soon-to-be npz file as dictionary.
    npz_dict = {
        'type': 'd', 'name': npz_name, 'z': Z, 'R': R, 'r_unit': r_unit,
        'E': E, 'e_unit': e_unit, 'F': F, 'entity_ids': entity_ids,
        'comp_ids': comp_ids, 'theory': theory, 'r_prov_ids': r_prov_ids,
        'r_prov_specs': r_prov_specs
    }
    for label,arr in zip(('E', 'F'), (E, F)):
        npz_dict[f'{label}_min'] = np.min(arr.flatten())
        npz_dict[f'{label}_mean'] = np.mean(arr.flatten())
        npz_dict[f'{label}_var'] = np.var(arr.flatten())
        npz_dict[f'{label}_max'] = np.max(arr.flatten())

    # We explicitly turn everything into arrays (as specified by the npz format).
    # This is import to avoid having reptar make r_prov_ids into a nested group.
    for k,v in npz_dict.items():
        npz_dict[k] = np.array(v)

    # Save the npz file.
    npz_file = File(npz_path, mode='w', from_dict=npz_dict)
    npz_file.save()

PDB
===

Protein data bank (PDB) is a standardized file format for atomic coordinates with proteins in mind.
You can use :func:`~reptar.writers.write_pdb` to write these files as shown in the example below.

.. code-block:: python

    from reptar import File
    from reptar.writers import write_pdb

    rfile_path = '1h2o_120meoh_md.exdir'
    pdb_path = '1h2o_120meoh_md.pdb'

    rfile = File(rfile_path, mode='r')
    Z = rfile.get('prod_1/atomic_numbers')
    R = rfile.get('prod_1/geometry')
    entity_ids = rfile.get('prod_1/entity_ids')
    comp_ids = rfile.get('prod_1/comp_ids')

    write_pdb(pdb_path, Z, R, entity_ids, comp_ids)

Schnetpack
==========

Schnetpack provides a framework for training and applying atomistic neural networks.
They have their own "wrapper" around an ASE database with slightly different requirements.
You can use :func:`~reptar.writers.write_schnetpack_db` to write these databases as shown in the example below.

.. code-block:: python

    from reptar import File
    from reptar.writers import write_schnetpack_db

    rfile_path = '1h2o_120meoh_md.exdir'
    db_path = '1h2o_120meoh_md_schnetpack.db'

    rfile = File(rfile_path, mode='r')
    Z = rfile.get('prod_1/atomic_numbers')
    R = rfile.get('prod_1/geometry')
    E = rfile.get('prod_1/energy_pot')  # Hartree
    E *= 27.21138602  # eV
    
    db = write_schnetpack_db(db_path, Z, R, energy=E, centering_function=None)

XYZ
===

Standard XYZ files can be written with :func:`~reptar.writers.write_xyz` as shown in the example below.

.. code-block:: python

    from reptar import File
    from reptar.writers import write_xyz

    rfile_path = '1h2o_120meoh_md.exdir'
    xyz_path = '1h2o_120meoh_md.xyz'

    rfile = File(rfile_path, mode='r')
    Z = rfile.get('prod_1/atomic_numbers')
    R = rfile.get('prod_1/geometry')
    
    write_xyz(xyz_path, Z, R)
