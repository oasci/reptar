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
    
    write_xyz_gap(xyz_path, lattice, Z, R, E)

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
