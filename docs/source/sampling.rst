========
Sampling
========

Reptar provides automated procedures for sampling structures from reptar files.
For example, sampling water trimers from a molecular dynamics simulation.
To directly add structures to a group you can use the driver function :func:`~reptar.sampler.add_structures_to_group`.

Examples
========

.. _30h2o sampling tutorial:

Water fragments from MD
-----------------------

One of the most straightforward way to generate and sample structures is with molecular dynamics (MD).
For example, suppose we want to study water trimers.
Running a MD simulation of just three water molecules would work.
However, we can increase the sampled configurational space by running a larger simulation (say with 30 water molecules) and sampling trimers from that trajectory.
There could be reasons for the former, but for the purpose here we will stick with the larger simulation.

Running the MD simulation
^^^^^^^^^^^^^^^^^^^^^^^^^

First, we generate a structure with 30 water molecules using `packmol <http://leandro.iqm.unicamp.br/m3g/packmol/examples.shtml>`__.
An example input and water xyz structure is shown below.

.. tab-set::

    .. tab-item:: packmol.in

        .. code-block:: text

            tolerance 2.0
            output ./30h2o.pm.xyz
            filetype xyz

            structure ./scratch/h2o.xyz
                number 30
                inside sphere 0.0 0.0 0.0 5.0
            end structure

    .. tab-item:: h2o.xyz

        .. code-block:: text

            3

            O   0.000000000   0.000000000   0.216072000
            H   0.000000000   0.761480000  -0.372151000
            H   0.000000000  -0.761480000  -0.372151000
    
    .. tab-item:: 30h2o.pm.xyz

        .. code-block:: text

            90
            Built with Packmol                                             
            O            3.323459        1.249671       -4.065936
            H            3.029345        1.660272       -4.884937
            H            3.354534        0.308457       -4.263441
            O           -0.948160       -1.841487        4.729458
            H           -1.554566       -2.578312        4.606097
            H           -0.951151       -1.680481        5.678103
            O            4.885123       -1.038823       -2.896730
            H            4.173353       -1.288536       -2.299333
            H            5.182515       -0.181897       -2.575648
            O            0.098903       -1.741661       -2.219922
            H           -0.325547       -1.678885       -3.081177
            H            1.002412       -2.010268       -2.413258
            O            3.356237       -0.552282        2.808103
            H            2.953066        0.250765        3.152234
            H            4.264181       -0.523614        3.125394
            O           -2.236816       -0.856927        0.307986
            H           -2.433326       -1.489435       -0.389994
            H           -2.954838       -0.217759        0.265958
            O           -0.895394        1.030811        2.950351
            H           -0.347410        0.525075        3.558465
            H           -1.601780        1.386667        3.498296
            O            2.922327        4.301236       -2.890033
            H            2.397535        4.526351       -2.115582
            H            2.679924        3.392224       -3.092014
            O           -2.617095       -0.192384       -2.318831
            H           -2.200177        0.423711       -1.708532
            H           -2.877875        0.348977       -3.070350
            O           -5.194295       -0.681161        0.600516
            H           -5.040113       -0.680578        1.550298
            H           -5.161570       -1.609498        0.349564
            O           -3.554735       -3.596293        1.901365
            H           -2.770678       -4.148352        1.981015
            H           -3.684480       -3.229973        2.781612
            O           -1.579676        3.644531        1.519741
            H           -1.657722        4.074528        2.376987
            H           -2.462826        3.310782        1.333959
            O            3.462139        3.860784        2.319864
            H            4.031608        3.353893        2.906912
            H            2.893939        4.366371        2.909261
            O           -1.249907        4.404374       -1.213364
            H           -0.431416        4.655915       -1.652291
            H           -1.022925        3.617448       -0.708305
            O            0.778562        0.689632       -3.124996
            H            0.806418        0.933984       -4.055250
            H            1.567071        0.154825       -2.990485
            O            1.244938        2.734734        3.134846
            H            1.875846        2.125228        2.739477
            H            0.979293        2.310396        3.956569
            O           -5.307111       -0.162279       -2.837538
            H           -5.478044        0.682540       -2.409846
            H           -4.775662       -0.652839       -2.202895
            O            1.273082        2.390312       -1.727962
            H            0.380325        2.560121       -2.044205
            H            1.366355        2.950358       -0.951105
            O           -2.218361       -2.193975       -4.087941
            H           -2.363204       -1.354637       -4.535592
            H           -2.127605       -2.836059       -4.798820
            O           -1.087082       -2.002095        2.536885
            H           -1.832486       -1.402039        2.637696
            H           -1.019609       -2.148617        1.588288
            O            0.331207        0.527784       -0.915313
            H            0.983252        0.707597       -0.230943
            H           -0.329177       -0.022684       -0.483185
            O            3.410300        4.849715       -0.204952
            H            2.982225        5.146897        0.603931
            H            3.418472        3.889930       -0.137100
            O           -1.694084        2.353142       -2.751900
            H           -1.473875        2.499691       -3.677043
            H           -2.586207        2.700787       -2.656413
            O            1.479763       -1.737359        4.444174
            H            1.950275       -2.461670        4.020086
            H            2.170549       -1.125912        4.717680
            O           -4.664876        2.814349       -1.837517
            H           -5.188697        2.546851       -1.075995
            H           -4.402468        3.720054       -1.645948
            O           -0.362427       -3.234498       -0.251277
            H           -1.037745       -3.282916       -0.934988
            H            0.470442       -3.316666       -0.726079
            O            2.478084       -2.047483        1.381657
            H            1.709939       -1.777160        1.894235
            H            2.222336       -2.885622        0.984198
            O            0.683742       -1.164626       -5.112082
            H           -0.254481       -0.981727       -5.222282
            H            0.747140       -2.124685       -5.123241
            O            4.375247        1.552407        1.710137
            H            3.863619        1.611366        0.897353
            H            5.002935        0.840667        1.551103
            O            3.551634        1.732402       -1.283292
            H            3.046326        0.936853       -1.477258
            H            4.026980        1.923029       -2.097887

.. raw:: html

    <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>

    <div style="height: 300px; width: 400px; margin: auto;"
    class='viewer_3Dmoljs' data-datatype='xyz'
    data-backgroundcolor='0xffffff'
    data-href='./30h2o.pm.xyz'
    data-style='stick'
    data-spin='axis:y;speed:0.1'>
    </div>
    <!--
        Change data-href from ./30h2o.pm.xyz to
        https://raw.githubusercontent.com/aalexmmaldonado/reptar/main/docs/source/files/30h2o-md/30h2o.pm.xyz
        for local development
    -->

From here, we can run a MD simulation in `xtb <https://xtb-docs.readthedocs.io/en/latest/contents.html>`__ using the GFN2-xTB hamiltonian.
First, we need to optimize the structure.

.. tab-set::

    .. tab-item:: xtb

        .. code-block:: bash

            xtb ./30h2o.pm.xyz --opt normal --gfn 2 --charge 0 --cycles 1000

    .. tab-item:: 30h2o-gfn2-opt.xyz

        .. code-block:: text

            90
             energy: -152.624675495341 gnorm: 0.000945854334 xtb: 6.5.1 (579679a)
            O            2.25120289029939        1.10108143576949       -4.02490867553866
            H            2.66294234216607        1.80984447363591       -3.50821402297169
            H            2.75775172096255        0.28312631261036       -3.81184664931766
            O           -0.95256449623638       -0.57392963214208        3.45183518767458
            H           -0.84818456826293       -1.27362346980566        2.78756646733230
            H           -0.06148540282055       -0.45371698350298        3.84742126645326
            O            3.30300434795651       -1.10659167453055       -2.94516323501285
            H            2.50827286536011       -1.65204688619153       -2.83746926056623
            H            3.41288267550615       -0.65833160418960       -2.08589020357239
            O            0.73524536660530       -2.03389451584690       -2.34705000465398
            H            0.22950049335087       -2.13006373331815       -3.16672324651899
            H            0.60612360671740       -1.09688884285503       -2.09200010262439
            O            3.73547285800594       -1.23663585026649        3.12025014527739
            H            4.10213805514668       -0.56938145960607        2.50152074503866
            H            4.38861131532885       -1.36474584433919        3.80825090078147
            O           -1.81012150395619       -0.00648839549023       -0.00676117393343
            H           -1.83614510967115       -0.15163181021211       -0.97347976974157
            H           -2.76467357464035       -0.01721187915043        0.27755949121795
            O           -1.42257911345813        2.01266847979891        3.20025413175295
            H           -1.29388715955063        1.02857953603097        3.22346862014552
            H           -1.70766857128757        2.26003440915427        4.07956713342920
            O            2.89345484914112        2.75544214354230       -1.88778150860957
            H            3.49657820435273        3.34480693042225       -1.39103087940482
            H            1.99415867360036        3.03010993894717       -1.63790031345880
            O           -2.13936233892799       -0.12779368044202       -2.67146156366043
            H           -1.28782451455160        0.36991604862680       -2.67275109412832
            H           -2.85147367181497        0.55311262426280       -2.62860331932758
            O           -4.35043711206675       -0.18803348186715        0.53666635986242
            H           -4.27276962068296       -0.60328109463225        1.42516151909861
            H           -4.51967369685446       -0.91788464584039       -0.08696820169943
            O           -3.62143015311021       -1.47047047578363        2.76864023592065
            H           -3.05937169383121       -2.10902585947174        2.30288025407426
            H           -3.00285349442243       -0.98591698576825        3.32736636146122
            O           -1.68193979600845        2.63367702238055        0.57481229392192
            H           -1.67455372296018        2.61518917959202        1.54440705361685
            H           -1.71384200336031        1.68913517246345        0.31698659797033
            O            3.82284953073164        3.10403101666555        2.11442089534857
            H            4.20377423070208        2.21886676813607        1.99514649071269
            H            2.87811417737352        2.95230264281379        2.27999405249979
            O           -3.34237010721480        3.82657814108370       -1.20388616826500
            H           -3.92002313187992        4.42292863266960       -0.72756042913908
            H           -2.70465103287727        3.47484593836176       -0.54969323785201
            O            0.37176672598282        0.60411750058628       -2.17437564625062
            H            0.98674660444200        0.78542370102518       -2.91909300329630
            H            0.37934684496257        1.42632334131966       -1.64845216681300
            O            1.24491587793962        2.21358395108296        2.47070962358131
            H            1.14669514832218        1.55928928010686        1.74547653191014
            H            0.33909145285392        2.41758959057187        2.74729607553693
            O           -4.33423096773400       -1.70622627644674       -1.75392089128073
            H           -5.01250701265737       -1.25882346372627       -2.26329010892846
            H           -3.49196054583264       -1.36147163604393       -2.08154871983191
            O            0.26148299558287        3.14141335569352       -1.23464523917818
            H           -0.34338241731918        3.17843142899437       -1.99898527221953
            H           -0.31355250388024        3.12226135805265       -0.45046716615753
            O           -0.96708408513553       -1.62849503994127       -4.55498482503959
            H           -1.54174130743639       -1.14864389802236       -3.92217541480238
            H           -1.53957989665475       -2.12980446671062       -5.13512399545743
            O           -1.49804942404128       -2.41247673586213        1.29292423921557
            H           -1.53848751128852       -1.60415746094457        0.75274685926920
            H           -0.82779260314129       -2.97319022425615        0.88337573828020
            O            0.85719964426450        0.23387363952974        0.71252241446249
            H            1.00602703850156       -0.57499626181501        1.26087915547257
            H           -0.05854986199189        0.18069626264000        0.40605607490960
            O            4.35873952106253        4.25524377690058       -0.22308169160130
            H            4.17625363137259        3.85509230762832        0.66250508446915
            H            4.28203051423339        5.20228556081059       -0.10738214323595
            O           -1.40796522328477        2.58013934165332       -3.29418656305244
            H           -1.02181585877378        2.12508512725308       -4.07381985667279
            H           -2.10423384086937        3.15863084838944       -3.60269440686876
            O            1.51612072082741        0.15779691466020        4.27756411101432
            H            1.54322246274884        0.97299680436681        3.74178189246683
            H            2.22988421482998       -0.39961825892485        3.93863363348384
            O           -4.21724787153347        1.37493265367295       -1.97931408634760
            H           -4.43515295548269        1.03272616616675       -1.10307686429490
            H           -3.93043384737545        2.29743600397428       -1.83575161770192
            O            0.64525601948158       -3.57831029609696       -0.15476829394553
            H            0.81163208347916       -4.47447153418486       -0.44570614881724
            H            0.67871594537531       -3.01317639185311       -0.95803880889053
            O            1.29974399540190       -1.99397445454454        2.03393864565156
            H            2.19537069792012       -1.96133983252224        2.39323993100106
            H            1.26863605070080       -2.67500848155058        1.34772970336554
            O           -0.24319075963198        1.09095727661368       -5.23111787815875
            H            0.69992899646811        1.20311845702257       -5.03476253523092
            H           -0.40827133017173        0.14132675524140       -5.14300596708372
            O            4.73472962548377        0.57590137941139        1.38896193828557
            H            4.13645155041177        0.50702744229083        0.59124869009581
            H            5.62944529076219        0.57647775859613        1.04802769600075
            O            3.14555883678276        0.41175812183446       -0.67199572170493
            H            2.26301861606164        0.31049729366650       -0.25780564887382
            H            3.11430409318576        1.28740326237157       -1.13023450450887

.. raw:: html

    <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>

    <div style="height: 300px; width: 400px; margin: auto;"
    class='viewer_3Dmoljs' data-datatype='xyz'
    data-backgroundcolor='0xffffff'
    data-href='./30h2o-gfn2-opt.xyz'
    data-style='stick'
    data-spin='axis:y;speed:0.1'>
    </div>
    <!--
        Change data-href from ./30h2o-gfn2-opt.xyz to
        https://raw.githubusercontent.com/aalexmmaldonado/reptar/main/docs/source/files/30h2o-md/30h2o-gfn2-opt.xyz
        for local development
    -->

Now we can run MD with this optimized structure.
The simulation will run at 300 K for 50 ps with a time step of 1.0 fs.
On 24 cores this takes about an hour.

.. tab-set::

    .. tab-item:: xtb

        .. code-block:: bash

            xtb ./30h2o-gfn2-opt.xyz --md --input ./md.inp --gfn 2 --charge 0 --verbose > 30h2o-gfn2-md.out

    .. tab-item:: md.inp

        .. code-block:: text

            $md
                temp   = 300.0  # Temperature set point in K.
                time   =  50.0  # Duration in ps.
                dump   =  10.0  # Store coordinates every __ fs.
                step   =   1.0  # Time step in fs.
                velo   = false  # Write out velocities?
                nvt    = true   # Run NVT?
                hmass  =   0    # Scale hydrogen mass. 0 means normal H mass.
                shake  =   0    # Use SHAKE algorithm to constrain bonds.
                sccacc =   2.0  # SCC accuracy. Defaults to 1.0. Lower is better.
            $end
            # Confine the system to avoid dissociation.
            $wall
                potential = logfermi
                temp = 300  # logfermi temperature. Scales the strength.
                beta = 4  # Specifies the steepness of the potential.
                sphere: auto, all  # Automatically determine radius and confine all atoms.
            $end

Store data in exdir file
^^^^^^^^^^^^^^^^^^^^^^^^

After running the MD simulation, we stored the data in :download:`this exdir file<./files/30h2o-md/30h2o-gfn2-md.exdir.zip>` using the following script.

.. code-block:: python

    import os
    from reptar import File
    from reptar.creator import creator
    from reptar.utils import gen_entity_ids, gen_comp_ids

    rfile_path = './30h2o-gfn2-md.exdir'
    group_key = '/30h2o'

    out_path = '30h2o-gfn2-md.out'
    geom_path = '30h2o-gfn2-opt.xyz'
    traj_path = 'xtb.trj'

    charge = 0
    multiplicity = 1

    # Prepare entity and comp IDs
    atoms_per_molecule = 3
    num_molecules = 30
    comp_id = 'h2o'
    entity_ids = gen_entity_ids(
        atoms_per_molecule, num_molecules
    )
    comp_ids = gen_comp_ids(
        label=comp_id, num_mol=num_molecules, entity_ids=entity_ids
    )

    rfile = File(rfile_path, mode='a', allow_remove=False)
    create = creator(rfile=rfile)

    create.from_calc(
        group_key, out_path=out_path, geom_path=geom_path, traj_path=traj_path
    )
    rfile.put(f'{group_key}/entity_ids', entity_ids)
    rfile.put(f'{group_key}/comp_ids', comp_ids)
    rfile.put(f'{group_key}/charge', charge)
    rfile.put(f'{group_key}/mult', multiplicity)

With the data now in a useable format we can be sampling.

Water trimers
^^^^^^^^^^^^^

The following script will sample and store 5000 random trimers from the entire MD simulation.

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
        sample_comp_ids, structure_idxs=structure_idxs,
        center_structures=center_structures, sampling_updates=True,
        copy_EG=False, write=True
    )

Water dimers with criteria
^^^^^^^^^^^^^^^^^^^^^^^^^^

Sometimes we are only interested in structures that meet some criteria.
For example, suppose we only wanted dimers where the center of mass distance is below 5 Angstroms.
We can use the :func:`~reptar.descriptors.com_distance_sum` descriptor and the :class:`~reptar.descriptors.Criteria` class to accept or reject structures.

.. code-block:: python

    import os
    from reptar import File
    from reptar.sampler import add_structures_to_group
    from reptar.descriptors import Criteria, com_distance_sum

    rfile_path = './30h2o-gfn2-md.exdir'
    group_key = '/30h2o'
    sample_key = f'{group_key}/samples_2h2o'
    sample_comp_ids = ['h2o', 'h2o']  # Component IDs to sample (3 waters in this case).
    quantity = 5000  # Number of structures to sample.
    center_structures = True  # Translate center of mass to origin.
    structure_idxs = None  # None means we can sample from all structures.

    # Ensures we execute from script directory (for relative paths).
    os.chdir(os.path.dirname(os.path.realpath(__file__)))

    # Load exdir file.
    rfile = File(rfile_path, mode='a', allow_remove=False)

    # Initialize group.
    rfile.create_group(sample_key)

    # Initialize criteria.
    # Note that the com_distance_sum function requires entity ids.
    # However, the add_structures_to_group and sample_structures functions
    # automatically include the entity_ids.
    criteria = Criteria(com_distance_sum, {}, 5.0)

    # Sample structures and automatically create the group.
    add_structures_to_group(
        rfile, group_key, rfile, sample_key, quantity,
        sample_comp_ids, structure_idxs=structure_idxs, criteria=criteria,
        center_structures=center_structures, sampling_updates=True,
        copy_EG=False, write=True
    )
