========
Sampling
========

Reptar provides automated procedures for sampling structures from reptar files.
For example, sampling water trimers from a molecular dynamics simulation.
To directly add structures to a group you can use the driver function :func:`reptar.sampler.add_structures_to_group`.

Examples
========

Water fragments from MD
-----------------------

We can sample :math:`n`-body structures, or fragments, from molecular dynamics (MD) simulations.
For example, below is the input file for a NVT MD simulation of 30 water molecules at 300 K driven by the GFN2-xTB Hamiltonian.

.. raw:: html

    <script src="http://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>

    <div style="height: 500px; width: 500px;"
    class='viewer_3Dmoljs' data-datatype='xyz'
    data-backgroundcolor='0xffffff'
    data-href='./30h2o-gfn2-opt.xyz'
    data-style='stick'>
    </div>

.. code-block:: text

    $md
        temp     = 300.0  # in K
        time     =  50.0  # in ps
        dump     =  10.0  # in fs
        step     =   1.0  # in fs
        velo     = false
        nvt      = true
        hmass    =   0
        shake    =   0
        sccacc   =   2.0
    $end
    $wall
        potential = logfermi
        temp = 300
        beta = 4
        sphere: auto, all
    $end

After running the MD simulation in `xtb <https://xtb-docs.readthedocs.io/en/latest/contents.html>`__, we stored the data in :download:`this exdir file<./files/30h2o-md/30h2o-gfn2-md.exdir.zip>`.
We can then sample water trimers using the following Python script.

.. code-block:: python

    import os
    from reptar import File
    from reptar.sampler import add_structures_to_group

    rfile_path = './30h2o-gfn2-md.exdir'
    group_key = '/30h2o'
    sample_key = f'{group_key}/samples_3h2o'
    sample_comp_ids = ['h2o', 'h2o', 'h2o']  # Component IDs to sample (3 waters in this case).
    quantity = 5000  # Number of structures to sample.
    center_structures = True  # Translate center of mass to origin.
    structure_idxs = None  # None means we can sample from all structures.

    # Ensures we execute from script directory (for relative paths).
    os.chdir(os.path.dirname(os.path.realpath(__file__)))

    # Load exdir file.
    rfile = File(rfile_path, mode='a', allow_remove=False)

    # Initialize group.
    rfile.create_group(sample_key)

    # Sample structures and automatically create the group.
    add_structures_to_group(
        rfile, group_key, rfile, sample_key, quantity,
        sample_comp_ids, structure_idxs=structure_idxs, criteria=None,
        z_slice=[], cutoff=[], center_structures=center_structures,
        sampling_updates=True, copy_EG=False, write=True
    )

