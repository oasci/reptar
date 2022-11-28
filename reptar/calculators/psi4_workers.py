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

try:
    import psi4
except ImportError:
    pass

# pylint: disable=invalid-name
def psi4_energy(
    idxs, Z, R, charge=0, mult=1, method="mp2", options=None, threads=1, mem="1 GB"
):
    r"""Worker function for computing total electronic energy using Psi4.

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
        Total molecular multiplicity.
    method : :obj:`str`, default: ``'b3lyp-d3bj'``
        Specifies the Psi4 method used for the gradient. For more information,
        please see `the Psi4 documentation <https://psicode.org/psi4manual/master/
        opt.html#geometry-optimization-w-w-optimize-and-gradient>`__
        for your specific version.
    options : :obj:`dict`
        `Psi4 control keywords <https://psicode.org/psi4manual/master/
        psithoninput.html#job-control-keywords>`__
        using the PsiAPI format. Some common ones are ``basis``,
        ``e_convergence``, ``d_convergence``, and ``reference``.
    threads : :obj:`int`, default: ``1``
        Number of threads for Psi4. This is almost always the number of cores
        being used for the worker. For more information, see the
        `documentation <https://psicode.org/psi4manual/master/api/
        psi4.core.set_num_threads.html#psi4.core.set_num_threads>`__.
    mem : :obj:`int`, :obj:`float`, :obj:`str`, default: ``'1 GB'``
        The amount of memory available. For more information, see the
        `documentation <https://psicode.org/psi4manual/master/api/
        psi4.driver.set_memory.html#psi4.driver.set_memory>`__.

    Returns
    -------
    :obj:`numpy.ndarray`
        ``idxs``
    :obj:`numpy.ndarray`
        Total electronic energy of computed structures in the same order as
        ``idxs``. Units of Hartree.

    Notes
    -----
    Psi4 uses QCElemental to build molecules from arrays.
    There is some postprocessing in the `from_arrays routine
    <http://docs.qcarchive.molssi.org/projects/QCElemental/en/latest/api/
    qcelemental.molparse.from_arrays.html#from-arrays>`__.
    QCElemental will translate and rotate molecules which affects computed gradients.
    Here, we set ``fix_com`` and ``fix_orientation`` to ``True`` to avoid this.
    """
    psi4.core.set_num_threads(threads)
    psi4.set_memory(mem)
    if options is not None:
        psi4.set_options(options)
    R = R[idxs]
    E = np.empty(R.shape)
    for i, r in enumerate(R):
        mol = psi4.core.Molecule.from_arrays(
            geom=r,
            elez=Z,
            molecular_charge=charge,
            molecular_multiplicity=mult,
            fix_com=True,
            fix_orientation=True,
        )
        E[i] = psi4.energy(method, molecule=mol, return_wfn=False)
    return idxs, E


def psi4_engrad(
    idxs, Z, R, charge=0, mult=1, method="mp2", options=None, threads=1, mem="1 GB"
):
    r"""Worker function for computing total electronic energy and atomic
    gradients using Psi4.

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
        Total molecular multiplicity.
    method : :obj:`str`, default: ``'b3lyp-d3bj'``
        Specifies the Psi4 method used for the gradient. For more information,
        please see `the Psi4 documentation <https://psicode.org/psi4manual/master/
        opt.html#geometry-optimization-w-w-optimize-and-gradient>`__
        for your specific version.
    options : :obj:`dict`
        `Psi4 control keywords <https://psicode.org/psi4manual/master/psithoninput.html
        #job-control-keywords>`__
        using the PsiAPI format. Some common ones are ``basis``,
        ``e_convergence``, ``d_convergence``, and ``reference``.
    threads : :obj:`int`, default: ``1``
        Number of threads for Psi4. This is almost always the number of cores
        being used for the worker. For more information, see the
        `documentation <https://psicode.org/psi4manual/master/api/
        psi4.core.set_num_threads.html#psi4.core.set_num_threads>`__.
    mem : :obj:`int`, :obj:`float`, :obj:`str`, default: ``'1 GB'``
        The amount of memory available. For more information, see the
        `documentation <https://psicode.org/psi4manual/master/api/
        psi4.driver.set_memory.html#psi4.driver.set_memory>`__.

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

    Notes
    -----
    Psi4 uses QCElemental to build molecules from arrays.
    There is some postprocessing in the `from_arrays routine
    <http://docs.qcarchive.molssi.org/projects/QCElemental/en/latest/api/
    qcelemental.molparse.from_arrays.html#from-arrays>`__.
    QCElemental will translate and rotate molecules which affects computed gradients.
    Here, we set ``fix_com`` and ``fix_orientation`` to ``True`` to avoid this.
    """
    psi4.core.set_num_threads(threads)
    psi4.set_memory(mem)
    if options is not None:
        psi4.set_options(options)
    R = R[idxs]
    G = np.zeros(R.shape)
    E = np.zeros(R.shape[0])
    for i, r in enumerate(R):
        mol = psi4.core.Molecule.from_arrays(
            geom=r,
            elez=Z,
            molecular_charge=charge,
            molecular_multiplicity=mult,
            fix_com=True,
            fix_orientation=True,
        )
        g, wfn_mp2 = psi4.gradient(method, molecule=mol, return_wfn=True)
        g = g.to_array()
        g /= psi4.constants.bohr2angstroms
        G[i] = g
        E[i] = wfn_mp2.energy()
    return idxs, E, G


def psi4_opt(
    idxs, Z, R, charge=0, mult=1, method="mp2", options=None, threads=1, mem="1 GB"
):
    r"""Worker function for optimizations using Psi4.

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
        Total molecular multiplicity.
    method : :obj:`str`, default: ``'b3lyp-d3bj'``
        Specifies the Psi4 method used for the gradient. For more information,
        please see `the Psi4 documentation <https://psicode.org/psi4manual/master/
        opt.html#geometry-optimization-w-w-optimize-and-gradient>`__
        for your specific version.
    options : :obj:`dict`
        `Psi4 control keywords <https://psicode.org/psi4manual/master/
        psithoninput.html#job-control-keywords>`__
        using the PsiAPI format. Some common ones are ``basis``,
        ``e_convergence``, ``d_convergence``, and ``reference``.
    threads : :obj:`int`, default: ``1``
        Number of threads for Psi4. This is almost always the number of cores
        being used for the worker. For more information, see the
        `documentation <https://psicode.org/psi4manual/master/api/
        psi4.core.set_num_threads.html#psi4.core.set_num_threads>`__.
    mem : :obj:`int`, :obj:`float`, :obj:`str`, default: ``'1 GB'``
        The amount of memory available. For more information, see the
        `documentation <https://psicode.org/psi4manual/master/api/
        psi4.driver.set_memory.html#psi4.driver.set_memory>`__.

    Returns
    -------
    :obj:`numpy.ndarray`
        ``idxs``
    :obj:`numpy.ndarray`
        If the optimizations converged or not.
    :obj:`numpy.ndarray`
        Optimized geometries.
    :obj:`numpy.ndarray`
        Total electronic energies of optimized structures. Units of Hartree.
    :obj:`numpy.ndarray`
        Atomic gradients of optimized structures. Units of Hartree/Angstrom.

    Notes
    -----
    Psi4 uses QCElemental to build molecules from arrays.
    There is some postprocessing in the `from_arrays routine
    <http://docs.qcarchive.molssi.org/projects/QCElemental/en/latest/api/
    qcelemental.molparse.from_arrays.html#from-arrays>`__.
    QCElemental will translate and rotate molecules which affects computed gradients.
    Here, we set ``fix_com`` and ``fix_orientation`` to ``True`` to avoid this.
    """
    psi4.core.set_num_threads(threads)
    psi4.set_memory(mem)
    if options is not None:
        psi4.set_options(options)
    R = R[idxs]
    opt_conv = np.full(R.shape[0], False, dtype=np.bool8)
    R_opt = np.empty(R.shape, dtype=np.float64)
    G = np.full(R.shape, np.nan, dtype=np.float64)
    E = np.full(R.shape[0], np.nan, dtype=np.float64)
    for i, r in enumerate(R):
        mol = psi4.core.Molecule.from_arrays(
            geom=r,
            elez=Z,
            molecular_charge=charge,
            molecular_multiplicity=mult,
            fix_com=True,
            fix_orientation=True,
        )
        try:
            e, wfn = psi4.opt(method, molecule=mol, return_wfn=True)
            r_opt = np.asarray(wfn.molecule().geometry())
            g = np.asarray(wfn.gradient())
            r_opt_conv = True
        except psi4.OptimizationConvergenceError as ex:
            r_opt_conv = False
            r_opt = np.asarray(ex.wfn.molecule().geometry())
            e = ex.wfn.energy()
            g = np.asarray(ex.wfn.gradient())
        finally:
            r_opt *= psi4.constants.bohr2angstroms
            g /= psi4.constants.bohr2angstroms

            opt_conv[i] = r_opt_conv  # pylint: disable=used-before-assignment
            R_opt[i] = r_opt
            E[i] = e  # pylint: disable=used-before-assignment
            G[i] = g

    return idxs, opt_conv, R_opt, E, G
