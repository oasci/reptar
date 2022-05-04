# MIT License
# 
# Copyright (c) 2022, Alex M. Maldonado
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import numpy as np
from xtb.interface import Calculator, Param

def xtb_engrad(
    Z, R, R_idxs, charge, mult, calc_acc=0.1, params=Param.GFN2xTB
):
    """Ray remote function for computing total electronic energy and atomic
    gradients using xtb.

    Parameters
    ----------
    Z : ``ray.ObjectRef`` of :obj:`numpy.ndarray`
        Atomic numbers of the atoms with repsect to ``R``.
    R : ``ray.ObjectRef`` of :obj:`numpy.ndarray`
        Cartesian coordinates of all structures in group. This includes
        unused structures.
    R_idxs : :obj:`numpy.ndarray`
        Indices of the structures from ``R`` to compute energies and gradients
        for.
    charge : :obj:`int`
        Total molecular charge.
    mult : :obj:`int`
        Total multiplicity.
    calc_acc : :obj:`str`, optional
        Numerical accuracy for calculation. For more information, visit the
        `documentation <https://xtb-python.readthedocs.io/en/latest\
        /general-api.html#xtb.interface.Calculator.set_accuracy>`_.
        Defaults to ``0.1``.
    params, optional
        xTB parameters. Defaults to ``Param.GFN2xTB``.
    
    Returns
    -------
    :obj:`numpy.ndarray`
        ``R_idxs``
    :obj:`numpy.ndarray`
        Total electronic energy of computed structures in the same order as
        ``R_idxs``. Units of Hartree.
    :obj:`numpy.ndarray`
        Atomic gradients of computed structures in the same order as ``R_idxs``.
        Units of Hartree/Angstrom.
    """
    n_upair_ele = int(mult - 1)
    R = R[R_idxs]
    G = np.zeros(R.shape)
    E = np.zeros(R.shape[0])
    for i in range(len(R)):
        calc = Calculator(Param.GFN2xTB, Z, R[i], charge, n_upair_ele)
        calc.set_accuracy(calc_acc)
        res = calc.singlepoint()
        g = res.get_gradient()
        g /= 0.52917721067  # psi4.constants.bohr2angstroms
        G[i] = g
        E[i] = res.get_energy()
    return R_idxs, E, G
