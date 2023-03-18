# MIT License
#
# Copyright (c) 2022-2023, Alex M. Maldonado
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
from ..logger import ReptarLogger

log = ReptarLogger(__name__)

try:
    from xtb.interface import Calculator, Param

    _HAS_XTB = True
except ImportError:
    _HAS_XTB = False


def xtb_engrad(
    idxs, Z, R, charge=0, mult=1, calc_acc=0.1, max_iterations=300, params=None
):
    r"""Ray remote function for computing total electronic energy and atomic
    gradients using xtb.

    Parameters
    ----------
    idxs : :obj:`numpy.ndarray`, ndim: ``1``
        Indices of the structures from ``R`` to compute energies and gradients
        for.
    Z : :obj:`numpy.ndarray`, ndim: ``1``
        Atomic numbers of the atoms with respect to ``R``.
    R : :obj:`numpy.ndarray`, ndim: ``3``
        Cartesian coordinates of all structures in group. This includes
        unused structures.
    charge : :obj:`int`, default: ``0``
        Total molecular charge.
    mult : :obj:`int`, default: ``1``
        Total multiplicity.
    calc_acc : :obj:`int`, default: ``0.1``
        Numerical accuracy for calculation. For more information, visit the
        `documentation <https://xtb-python.readthedocs.io/en/latest/
        general-api.html#xtb.interface.Calculator.set_accuracy>`__.
    max_iterations : :obj:`int`, default: ``300``
        Maximum number of iterations for self-consistent charge methods. If the
        calculations fails to converge in a given number of cycles, no error
        is necessarily shown.
    params : default: ``None``
        xTB parameters. Defaults to ``xtb.interface.Param.GFN2xTB`` if ``None``.

    Returns
    -------
    :obj:`numpy.ndarray`
        ``idxs``
    :obj:`numpy.ndarray`
        Total electronic energy of computed structures in the same order as
        ``idxs``. Units of Hartree.
    :obj:`numpy.ndarray`
        Atomic gradients of computed structures in the same order as ``idxs``.
        Units of Hartree/Angstrom.
    """
    assert _HAS_XTB

    n_upair_ele = int(mult - 1)
    R = R[idxs]
    G = np.zeros(R.shape)
    E = np.zeros(R.shape[0])
    if params is None:
        params = Param.GFN2xTB
    for i, r in enumerate(R):
        calc = Calculator(params, Z, r, charge, n_upair_ele)
        calc.set_accuracy(calc_acc)
        calc.set_max_iterations(max_iterations)
        res = calc.singlepoint()
        g = res.get_gradient()
        g /= 0.52917721067  # psi4.constants.bohr2angstroms
        G[i] = g
        E[i] = res.get_energy()
    return idxs, E, G
