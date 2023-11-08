==
md
==

runtime_info
------------

.. glossary::

    algo_shake
        If the SHAKE algorithm is used for bond constrains.

        **Type:** :obj:`bool`

    md_duration
        TODO in ps

        **Type:** :obj:`float`

    md_restarted
        Whether this MD simulation was restarted from another.

        **Type:** :obj:`bool`

    t_step
        TODO

        **Type:** :obj:`float`

    thermostat_temp
        Thermostat temperature set point.

        **Type:** :obj:`float`

    thermostat_type
        Type of thermostat used.

        **Type:** :obj:`string`

    md_steps_dump_traj
        Steps between each trajectory dump in femtoseconds.

        **Type:** :obj:`int`

outputs
-------

.. glossary::

    avg_energy_tot
        Average total energy during a MD simulation in Hartrees.

        **Type:** :obj:`float`

    avg_temp
        Average temperature during a MD simulation.

        **Type:** :obj:`float`

    energy_ke
        Kinetic energy of the system in Hartrees.

        **Type:** :obj:`float`

    temp
        Instantaneous temperature in Kelvin.

        **Type:** :obj:`float`

    velcs
        Atomic velocities in Angstrom/fs.

        **Type:** :obj:`float`
