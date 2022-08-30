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
import os
from ..utils import z_to_element
    
def write_pdb(
    pdb_path, Z, R, entity_ids, comp_ids, atom_type='HETATM'
):
    """Write PDB file containing a single structure or trajectory.

    Parameters
    ----------
    pdb_path : :obj:`str`
        Path to PDB file to write.
    Z : :obj:`numpy.ndarray`, ndim: ``1``
        Atomic numbers of each structure in ``R``.
    R : :obj:`numpy.ndarray`, ndim: ``3``
        Cartesian coordinates.
    entity_ids : :obj:`numpy.ndarray`, ndim: ``1``
        A uniquely identifying integer specifying what atoms belong to
        which entities. Entities can be a related set of atoms, molecules,
        or functional group. For example, a water and methanol molecule
        would be ``[0, 0, 0, 1, 1, 1, 1, 1, 1]``.
    comp_ids : :obj:`numpy.ndarray`, ndim: ``1``
        Relates ``entity_id`` to a fragment label for chemical components
        or species. Labels could be ``WAT`` or ``h2o`` for water, ``MeOH``
        for methanol, ``bz`` for benzene, etc. There are no standardized
        labels for species. The index of the label is the respective
        ``entity_id``. For example, a water and methanol molecule could
        be ``['h2o', 'meoh']``.
    atom_type : :obj:`str`, default: ``'HETATM'``
        The PDB atom type to be used. Should almost always be ``HETATM``.
    """ 
    element_symbols = [z_to_element[z] for z in Z]

    num_structures = len(R)
    num_atoms = len(Z)

    file_name = os.path.splitext(os.path.basename(pdb_path))[0]
    # Trims component ids to the first three letters
    comp_ids = np.array([comp_id[:3] for comp_id in comp_ids])
    with open(pdb_path, 'w') as f:
        f.write(file_name+'\n')
        f.write(f'{num_atoms}\n')  
        for i_structure in range(num_structures):
            if i_structure != 0:
                f.write('MODEL\n')
            for i_atom in range(num_atoms):
                entity_id = entity_ids[i_atom]
                comp_id = comp_ids[entity_id]
                element_symbol = element_symbols[i_atom]
                # Determines the number of atoms with the same atomic number
                # in the entity before this atom.
                atom_type_count = np.count_nonzero(
                    Z[
                        np.intersect1d(
                            np.arange(i_atom), np.argwhere(entity_ids == entity_id).T[0]
                        )
                    ] == Z[i_atom]
                )
                atom_label = element_symbol + str(atom_type_count+1)

                coords = R[i_structure][i_atom]

                f.write(
                    '{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2s}          {:>2s}{:2s}\n'.format(
                        str(atom_type), i_atom+1, str(atom_label), '',
                        str(comp_id), 'A', entity_id+1, '',
                        coords[0], coords[1], coords[2],
                        1.00, '', str(atom_label), ''
                    )
                )
            if num_structures > 1:
                f.write('ENDMDL\n')
