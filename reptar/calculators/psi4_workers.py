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

import os
from tempfile import TemporaryDirectory
import numpy as np
from .cube import initialize_grid_arrays
from ..logger import ReptarLogger
from ..parsers.gaussian_cube import parse_cube

log = ReptarLogger(__name__)

try:
    import psi4
    from optking.exceptions import AlgError
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
    E = np.empty(R.shape[0])
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
        g, wfn = psi4.gradient(method, molecule=mol, return_wfn=True)
        g = g.to_array()
        g /= psi4.constants.bohr2angstroms
        G[i] = g
        E[i] = wfn.energy()
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
    method : :obj:`str`, default: ``'mp2'``
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
            r_opt_conv = True
        except (psi4.OptimizationConvergenceError, AlgError) as ex:
            r_opt_conv = False
            r_opt = np.asarray(ex.wfn.molecule().geometry())
            e = ex.wfn.energy()
        finally:
            r_opt *= psi4.constants.bohr2angstroms

            opt_conv[i] = r_opt_conv  # pylint: disable=used-before-assignment
            R_opt[i] = r_opt
            E[i] = e  # pylint: disable=used-before-assignment

    return idxs, opt_conv, R_opt, E


def psi4_cube(
    idxs,
    Z,
    R,
    total_n_points,
    charge=0,
    mult=1,
    method="mp2",
    options=None,
    threads=1,
    mem="1 GB",
    cube_file_name=None,
):
    r"""Worker function for computing cube data using Psi4.

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
    total_n_points : :obj:`int`
        Number of points to initialize arrays with.
    charge : :obj:`int`, default: ``0``
        Total molecular charge.
    mult : :obj:`int`, default: ``1``
        Total molecular multiplicity.
    method : :obj:`str`, default: ``'b3lyp-d3bj'``
        Specifies the Psi4 method used for the calculation. For more information,
        please see `the Psi4 documentation <https://psicode.org/psi4manual/master/
        opt.html#geometry-optimization-w-w-optimize-and-gradient>`__
        for your specific version.
    options : :obj:`dict`
        `Psi4 control keywords <https://psicode.org/psi4manual/master/psithoninput.html
        #job-control-keywords>`__ using the PsiAPI format.

        You must specify ``cubeprop_tasks`` (the tasks are
        listed `here <https://psicode.org/psi4manual/master/cubeprop.html
        #cubeprop-tasks>`__). If ``cubeprop_filepath`` is not specified then no
        text-based cube file is saved.
    threads : :obj:`int`, default: ``1``
        Number of threads for Psi4. This is almost always the number of cores
        being used for the worker. For more information, see the
        `documentation <https://psicode.org/psi4manual/master/api/
        psi4.core.set_num_threads.html#psi4.core.set_num_threads>`__.
    mem : :obj:`int`, :obj:`float`, :obj:`str`, default: ``'1 GB'``
        The amount of memory available. For more information, see the
        `documentation <https://psicode.org/psi4manual/master/api/
        psi4.driver.set_memory.html#psi4.driver.set_memory>`__.
    cube_file_name : :obj:`str`, default: ``None``
        You can manually specify the name of the cube file to parse. Otherwise
        we will automatically choose. Typical options are ``ESP.cube``, ``Dt.cube``,


    Returns
    -------
    :obj:`numpy.ndarray`
        ``idxs``
    :obj:`numpy.ndarray`
        Cartesian coordinates of points where a property is probed.
    :obj:`numpy.ndarray`
        Property values at the same Cartesian coordinates.

    Notes
    -----
    Uses ```CubeProp()`` <https://psicode.org/psi4manual/master/cubeprop.html>`__ to
    generate a cube file, then parses back into an array.

    Psi4 uses QCElemental to build molecules from arrays.
    There is some postprocessing in the `from_arrays routine
    <http://docs.qcarchive.molssi.org/projects/QCElemental/en/latest/api/
    qcelemental.molparse.from_arrays.html#from-arrays>`__.
    QCElemental will translate and rotate molecules which affects computed gradients.
    Here, we set ``fix_com`` and ``fix_orientation`` to ``True`` to avoid this.
    """
    assert "cubeprop_tasks" in options
    if cube_file_name is None:
        if options["cubeprop_tasks"][0].lower() == "esp":
            cube_file_name = "ESP.cube"
        elif options["cubeprop_tasks"][0].lower() == "density":
            cube_file_name = "Dt.cube"

    with TemporaryDirectory() as temp_dir:
        options["cubeprop_filepath"] = temp_dir

        psi4.core.set_num_threads(threads)
        psi4.set_memory(mem)
        if options is not None:
            psi4.set_options(options)
        R = R[idxs]

        cube_R, cube_V = initialize_grid_arrays(R, max_points=total_n_points)

        for i, r in enumerate(R):
            mol = psi4.core.Molecule.from_arrays(
                geom=r,
                elez=Z,
                molecular_charge=charge,
                molecular_multiplicity=mult,
                fix_com=True,
                fix_orientation=True,
            )
            _, wfn = psi4.energy(method, molecule=mol, return_wfn=True)
            psi4.cubeprop(wfn)
            cube_r, cube_v = parse_cube(os.path.join(temp_dir, cube_file_name))

            cube_R[i, : cube_r.shape[0], :] = cube_r
            cube_V[i, : cube_v.shape[0]] = cube_v
    return idxs, cube_R, cube_V
