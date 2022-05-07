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
from ..utils import z_to_element
    
def write_pdb(
    Z, R, entity_ids, comp_ids, atom_type='HETATM', file_name=None,
    save_dir=None,
):
    """Write PDB file.

    Parameters
    ----------
    Z : :obj:`numpy.ndarray`
        Atomic numbers used to populate ``numbers`` in the ASE db.
    R : :obj:`numpy.ndarray`
        Cartesian coordinates used to populate ``positions`` in the ASE db.
    entity_ids : :obj:`numpy.ndarray`
        A uniquely identifying integer specifying what atoms belong to
        which entities. Entities can be a related set of atoms, molecules,
        or functional group. For example, a water and methanol molecule
        could be ``[0, 0, 0, 1, 1, 1, 1, 1, 1]``.
    comp_ids : :obj:`numpy.ndarray`
        Relates ``entity_id`` to a fragment label for chemical components
        or species. Labels could be ``WAT`` or ``h2o`` for water, ``MeOH``
        for methanol, ``bz`` for benzene, etc. There are no standardized
        labels for species. The index of the label is the respective
        ``entity_id``. For example, a water and methanol molecule could
        be ``['h2o', 'meoh']``.
    atom_type : :obj:`str`, optional
        The PDB atom type to be used. Should almost always be ``HETATM``.
    file_name : :obj:`str`, optional
        Name of the PDB file. Defaults to ``'structure'``.
    save_dir : :obj:`str`, optional
        Directory to save the PDB file. Defaults to the current directory.
    
    Returns
    -------
    :obj:`str`
        Path to PDB file.
    """ 
    atom_labels = [z_to_element[z] for z in Z]
    
    if file_name is None:
        file_name = 'structure'
    if save_dir is None:
        save_dir = '.'
    else:
        if save_dir[-1] == '/':
            save_dir = save_dir[:-1]

    num_structures = len(R)
    num_atoms = len(Z)

    # Trims component ids to the first three letters
    comp_ids = np.array([comp_id[:3] for comp_id in comp_ids])
    file_path = f'{save_dir}/{file_name}.pdb'
    with open(file_path, 'w') as f:
        f.write(file_name+'\n')
        f.write(f'{num_atoms}\n')  
        for i_structure in range(num_structures):
            if i_structure != 0:
                f.write('MODEL\n')
            for i_atom in range(num_atoms):
                entity_id = entity_ids[i_atom]
                comp_id = comp_ids[entity_id]
                atom_label = atom_labels[i_atom]
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
    return file_path
            
    
    