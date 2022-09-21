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

def write_ase_db(db_path, Z, R, energy=None, forces=None):
    """Create an atomic simulation environment database.

    Will add rows if the database already exists.

    Parameters
    ----------
    db_path : :obj:`str`
        Path to atomic simulation environment database.
    Z : :obj:`numpy.ndarray`, ndim: ``1``
        Atomic numbers used to populate ``numbers`` in the ASE db.
    R : :obj:`numpy.ndarray`, ndim: ``3``
        Cartesian coordinates used to populate ``positions`` in the ASE db.
    energy : :obj:`numpy.ndarray`, ndim: ``1``, default: ``None``
        Energies to population ``energy`` in the ASE db.
    forces : :obj:`numpy.ndarray`, ndim: ``3``, default: ``None``
        Atomic forces to populate ``forces`` in the ASE db.
    
    Returns
    -------
    ``ase.db.sqlite.SQLite3Database``
        Atomic simulation environment database.
    """
    from ase import Atoms
    from ase.db import connect
    from ase.calculators.singlepoint import SinglePointCalculator

    db = connect(db_path)
    for i in range(len(R)):
        if energy is not None:
            e = energy[i]
        else:
            e = None
        if forces is not None:
            f = forces[i]
        else:
            f = None
        # TODO: Other way to add energies and forces? This is just a workaround.
        spe = SinglePointCalculator(
            Atoms(Z, R[i]), **{'energy': e, 'forces': f}
        )
        atom = spe.get_atoms()
        db.write(atom)
    return db
