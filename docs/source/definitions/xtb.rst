===
xtb
===

runtime_info
------------

.. glossary::

    xtb_scc_accuracy
        SCC accuracy level in MD. Every 10th step the SCC is properly converged at sccconv=1.0. sccmd should be < 5 in critical cases, effects may show up as bad thermostating

        **Type:** :obj:`float`

    wall_potential
        Confining potential(s) for a MD simulation. Each potential is described with a dictionary (with its own definitions) in a list. ``sphere_radius`` is the radius of the confining potential in Angstroms.

        **Type:** :obj:`list`

    mass_h
        Mass of Hydrogen atoms in atomic mass units (if specified).

        **Type:** :obj:`float`
