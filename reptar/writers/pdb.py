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

from .writer import reptarWriter
from ..utils import z_to_element

class pdbWriter(reptarWriter):
    """Writes PDB file from reptar group.
    """

    def __init__(self):
        super().__init__()
    
    def write(
        self, data, group_key, atom_type='HETATM', file_name=None, save_dir=None,
        R_limits=None
    ):
        """Parses trajectory file and extracts information.

        Parameters
        ----------
        data : ``reptar.data``
            A reptar data object.
        group_key : :obj:`str`
            Key to the desired group.
        atom_type : :obj:`str`, optional
            The PDB atom type to be used. Should almost always be ``HETATM``.
        file_name : :obj:`str`, optional
            Name of the PDB file. Defaults to ``'structure'``.
        save_dir : :obj:`str`, optional
            Directory to save the PDB file. Defaults to the current directory.
        R_limits : :obj:`tuple` (:obj:`int`), optional
            Slicing limits specifying start and stop index. If ``None`` is used
            it does not change the beginning or end of the array. For example,
            ``(None, 100)`` would only print the first 100 structures.
        """
        R = data.get(f'{group_key}/geometry')
        if R.ndim == 2:
            R = np.array([R])
        if R_limits is not None:
            R = R[slice(*R_limits)]
        Z = data.get(f'{group_key}/atomic_numbers')
        atom_labels = [z_to_element[i] for i in Z]
        entity_ids = data.get(f'{group_key}/entity_ids')
        comp_ids = data.get(f'{group_key}/comp_ids')
        
        if file_name is None:
            file_name = 'structure'
        if save_dir is None:
            save_dir = '.'
        else:
            if save_dir[-1] == '/':
                save_dir = save_dir[:-1]

        num_structures = len(R)
        num_atoms = len(Z)
        with open(f'{save_dir}/{file_name}.pdb', 'w') as f:
            f.write(file_name+'\n')
            f.write(f'{num_atoms}\n')  
            for i_structure in range(num_structures):
                if i_structure != 0:
                    f.write('MODEL\n')
                for i_atom in range(num_atoms):
                    entity_id = entity_ids[i_atom]
                    comp_id = comp_ids[entity_id]  # TODO: Handle component ids that are longer than three letters.
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
                
    
    