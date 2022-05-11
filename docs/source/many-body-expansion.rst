===================
Many-body expansion
===================

The many-body expansion (MBE) is a computational framework that states a system's total energy (an its gradients) are equal to the sum of all :math:`n`-body contributions.

.. math::

    E = \sum_{i}^N E_i^{(1)} + \sum_{i < j}^N \Delta E_{ij}^{(2)} + \sum_{i < j < k}^N \Delta E_{ijk}^{(3)} + \cdots.

Here, :math:`N` is the number of monomers; :math:`i`, :math:`j`, :math:`k` are monomer indices; :math:`E_i^{(1)}` is the energy of an isolated monomer :math:`i`; and :math:`\Delta E_{\; i, \: j, \: \ldots}^{\; (n)}` represents the :math:`n`-body interaction energy contribution of the fragment containing monomers :math:`i`, :math:`j`, :math:`\ldots` with lower order (:math:`< n`) contributions removed.
For example, the 2-body contribution of the fragment containing monomers :math:`i` and :math:`j` is

.. math::

    \Delta E_{ij}^{(2)} = E_{ij}^{(2)} - E_{i}^{(1)} - E_{j}^{(1)}.

Reptar provides useful routines for computing :math:`n`-body properties.

*n*-body contributions
======================

.. autoclass:: reptar.mb.mb_contributions
    :noindex:
