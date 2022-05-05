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
import psi4

def psi4_engrad(Z, R, R_idxs, charge, mult, method, options, threads, mem):
    """Worker function for computing total electronic energy and atomic
    gradients using Psi4.

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
        Total molecular multiplicity.
    method : :obj:`str`
        Specifies the Psi4 method used for the gradient. For more information,
        please see `the Psi4 documentation <https://psicode.org/psi4manual/\
        master/opt.html#geometry-optimization-w-w-optimize-and-gradient>`_ for
        your specific version.
    options : :obj:`dict`
        `Psi4 control keywords <https://psicode.org/psi4manual/master\
        /psithoninput.html#job-control-keywords>`_ using the PsiAPI format.
        Some common ones are ``basis``, ``e_convergence``, ``d_convergence``,
        and ``reference``.
    threads : :obj:`int`
        Number of threads for Psi4. This is almost always the number of cores
        being used for the worker. For more information, see the
        `documentation <https://psicode.org/psi4manual/master/api\
        /psi4.core.set_num_threads.html#psi4.core.set_num_threads>`_.
    mem : :obj:`int`, :obj:`float`, :obj:`str`
        The amount of memory available. For more information, see the
        `documentation <https://psicode.org/psi4manual/master/api\
        /psi4.driver.set_memory.html#psi4.driver.set_memory>`_.
    
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
    psi4.core.set_num_threads(threads)
    psi4.set_memory(mem)
    psi4.set_options(options)
    R = R[R_idxs]
    G = np.zeros(R.shape)
    E = np.zeros(R.shape[0])
    for i in range(len(R)):
        mol = psi4.core.Molecule.from_arrays(
            geom=R[i], elez=Z, molecular_charge=charge,
            molecular_multiplicity=mult
        )
        g, wfn_mp2 = psi4.gradient(
            method, molecule=mol, return_wfn=True
        )
        g = g.to_array()
        g /= psi4.constants.bohr2angstroms
        G[i] = g
        E[i] = wfn_mp2.energy()
    return R_idxs, E, G
