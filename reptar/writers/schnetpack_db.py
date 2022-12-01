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

from ase import Atoms
import numpy as np

try:
    from schnetpack import AtomsData
    from schnetpack.data.atoms import get_center_of_mass
except ImportError:
    _HAS_SCHNETPACK = False
else:
    _HAS_SCHNETPACK = True


def write_schnetpack_db(
    db_path,
    Z,
    R,
    energy=None,
    forces=None,
    e_units=1.0,
    f_units=1.0,
    centering_function=None,
):
    r"""Create a schnetpack database.

    Will add rows if the database already exists.

    Parameters
    ----------
    db_path : :obj:`str`
        Path to schnetpack database.
    Z : :obj:`numpy.ndarray`, ndim: ``1``
        Atomic numbers used to populate ``numbers`` in the ASE db.
    R : :obj:`numpy.ndarray`, ndim: ``3``
        Cartesian coordinates used to populate ``positions`` in the ASE db.
    energy : :obj:`numpy.ndarray`, ndim: ``1``, default: ``None``
        Energies to add to database in eV.
    forces : :obj:`numpy.ndarray`, ndim: ``3``, default: ``None``
        Atomic forces to add to database in eV/A.
    e_units : :obj:`float`, default: ``1.0``
        Units of energy using ``ase.units``. Essentially the conversion factor
        to eV.
    f_units : :obj:`float`, default: ``1.0``
        Units of forces using ``ase.units``. Essentially the conversion factor
        to eV/A.
    centering_function : ``callable``
        A function for centering the ``positions`` when querying the database.
        Defaults to centering the structure according to the center of mass.
        Set to ``None`` for no centering. Defaults to
        ``schnetpack.data.atoms.get_center_of_mass``

    Returns
    -------
    ``schnetpack.data.atoms.AtomsData``
        Schnetpack database.

    Notes
    -----
    When using ``db.get_properties(idx)`` to query the database, the
    ``ase.atoms`` object has the exact same positions as it was originally
    written. The ``_properties`` key in the dictionary could be centered
    (if ``centering_function`` is not ``None``) and has some precision issues.

    For example, if ``R[0][0][0] = -0.93996473515199`` is written to the
    database, ``db``, then we can get the positions using
    ``row_atom, row_prop = db.get_properties(0)``.
    ``row_atom.positions[0][0]`` would be exactly ``-0.93996473515199``.
    However, ``row_prop["_positions"][0][0]`` could be something like
    ``-0.93996471166611``.
    """
    assert _HAS_SCHNETPACK

    if centering_function is None:
        centering_function = get_center_of_mass

    avail_props = []
    unit_list = []
    if energy is not None:
        avail_props.append("energy")
        unit_list.append(e_units)
    if forces is not None:
        avail_props.append("forces")
        unit_list.append(f_units)

    db = AtomsData(
        db_path,
        available_properties=avail_props,
        units=unit_list,
        centering_function=centering_function,
    )

    all_atoms = []
    all_properties = []
    for i, r in enumerate(R):
        all_atoms.append(Atoms(Z, r))

        props = {}
        if energy is not None:
            props["energy"] = np.array([energy[i]])
        if forces is not None:
            props["forces"] = forces[i]
        all_properties.append(props)
    db.add_systems(all_atoms, all_properties)
    return db
