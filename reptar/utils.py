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

import hashlib
import numpy as np
import collections.abc

element_to_z = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9,
    'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16,
    'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23,
    'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30,
    'Ga': 31,'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37,
    'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44,
    'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51,
    'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58,
    'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65,
    'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72,
    'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79,
    'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86,
    'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93,
    'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100,
    'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106,
    'Bh': 107, 'Hs': 108, 'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112,
    'Uuq': 114, 'Uuh': 116,
}
z_to_element = {v: k for k, v in element_to_z.items()}

# Standard atomic weight from version 4.1 of
# https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses
# If lower and upper bounds were provided the lower bound was selected.
z_to_mass = (
    None, 1.00784, 4.002602, 6.938, 9.0121831, 10.806, 12.0096, 14.00643,  # N
    15.99903, 18.998403163, 20.1797, 22.98976928, 24.304, 26.9815385, 28.084,  # Si
    30.973761998, 32.059, 35.446, 39.948, 39.0983, 40.078, 44.955908,  # Sc
    47.867, 50.9415, 51.9961, 54.938044, 55.845, 58.933194, 58.6934, 63.546,  # Cu
    65.38, 69.723, 72.630, 74.921595, 78.971, 79.901, 83.798, 85.4678, 87.62,  # Sr
    88.90584, 91.224, 92.90637, 95.95, 98.0, 101.07, 102.90550, 106.42, 107.8682,  # Ag
    112.414, 114.818, 118.710, 121.760, 127.60, 126.90447, 131.293,  # Xe
    132.90545196, 137.327, 138.90547, 140.116, 140.90766, 144.242, 145.0, 150.36,  # Sm
    151.964, 157.25, 158.92535, 162.500, 164.93033, 167.259, 168.93422, 173.054,  # Yb
    174.9668, 178.49, 180.94788, 183.84, 186.207, 190.23, 192.217, 195.084,  # Pt
    196.966569, 200.592, 204.382, 207.2, 208.98040, 209.0, 210.0, 222.0, 223.0, 226.0,  # Ra
    227.0, 232.0377, 231.03588, 238.02891, 237.0, 244.0  # Pu
)

def get_files(path, expression, recursive=True):
    """Returns paths to all files in a given directory that matches a provided
    expression in the file name.
    
    Commonly used to find all files of a certain type, e.g., output or xyz
    files.
    
    Parameters
    ----------
    path : :obj:`str`
        Specifies the directory to search.
    expression : :obj:`str`
        Expression to be tested against all file names in ``path``.
    recursive : :obj:`bool`, optional
        Recursively find all files in all subdirectories.
    
    Returns
    -------
    :obj:`list` [:obj:`str`]
        All absolute paths to files matching the provided expression.
    """
    if path[-1] != '/':
        path += '/'
    if recursive:
        all_files = []
        for (dirpath, _, filenames) in os.walk(path):
            index = 0
            while index < len(filenames):
                if dirpath[-1] != '/':
                    dirpath += '/'
                filenames[index] = dirpath + filenames[index]
                index += 1
            all_files.extend(filenames)
        files = []
        for f in all_files:
            if expression in f:
                files.append(f)
    else:
        files = []
        for f in os.listdir(path):
            filename = os.path.basename(f)
            if expression in filename:
                files.append(path + f)
    return files

def atoms_by_element(atom_list):
    """Converts a list of atoms identified by their atomic number to their
    elemental symbol in the same order.
    
    Parameters
    ----------
    atom_list : :obj:`list` [:obj:`int`]
        Atomic numbers of atoms within a structure.
    
    Returns
    -------
    :obj:`list` [:obj:`str`]
        Element symbols of atoms within a structure.
    """

    atom_list_elements = []
    for atom in atom_list:
        atom_list_elements.append(str(z_to_element[atom]))

    return atom_list_elements

def atoms_by_number(atom_list):
    """Converts a list of atoms identified by their elemental symbol to their
    atomic number.
    
    Parameters
    ----------
    atom_list : :obj:`list` [:obj:`str`]
        Element symbols of atoms within a structure.
    
    Returns
    -------
    :obj:`list` [:obj:`int`]
        Atomic numbers of atoms within a structure.
    """
    return [int(element_to_z[i]) for i in atom_list]

def parse_stringfile(stringfile_path):
    """Parses data from string file.

    A string file is data presented as consecutive xyz data. The data could be
    three Cartesian coordinates for each atom, three atomic force vector
    components, or both coordinates and atomic forces in one line (referred to
    as extended xyz).
    
    Parameters
    ----------
    stringfile_path : :obj:`str`
        Path to string file.
    
    Returns
    -------
    :obj:`tuple` [:obj:`list`]
        Parsed atoms (as element symbols :obj:`str`), comments, and data as
        :obj:`float` from string file.
    """
    Z, comments, data = [], [], []
    with open(stringfile_path, 'r') as f:
        for _, line in enumerate(f):
            line = line.strip()
            if not line:
                # Skips blank lines
                pass
            else:
                line_split = line.split()
                if len(line_split) == 1 \
                    and float(line_split[0]) % int(line_split[0]) == 0.0:
                    # Skips number of atoms line, adds comment line, and
                    # prepares next z and data item.
                    comment_line = next(f)
                    comments.append(comment_line.strip())
                    Z.append([])
                    data.append([])
                else:
                    # Grabs z and data information.
                    Z[-1].append(line_split[0])
                    data[-1].append([float(i) for i in line_split[1:]])
    return Z, comments, data

def get_md5(rfile, group_key, only_arrays=False, only_structures=False):
    """Creates MD5 hash for a group.

    Parameters
    ----------
    rfile : :obj:`reptar.File`
        A reptar file.
    group_key : :obj:`str`
        Key to the desired group.
    only_arrays : :obj:`bool`, optional
        Generate the MD5 hash using only arrays. This creates a data-centered
        MD5 that is not affected by data that are commonly added or
        changed. Defaults to ``False``.
    only_structures : :obj:`bool`, optional
        Generate the MD5 has with only ``atomic_numbers`` and ``geometry``
        if possible). This is more static than ``only_arrays`` and should be
        used to track sampling (i.e., ``r_prov_ids``).
    
    Returns
    -------
    :obj:`str`
        MD5 hash of a group.
    """
    md5_hash = hashlib.md5()
    # TODO: Figure out why different formats have different MD5s.

    if only_structures:
        try:
            Z = rfile.get(f'{group_key}/atomic_numbers')
            Z = Z.ravel()
            md5_hash.update(hashlib.md5(Z).digest())
        except Exception:
            pass
        try:
            R = rfile.get(f'{group_key}/geometry')
            R = R.ravel()
            md5_hash.update(hashlib.md5(R).digest())
        except Exception:
            pass
    else:
        keys = rfile.get_keys(group_key)
        for key in keys:
            if 'md5' in key:
                continue
            d = rfile.get(f'{group_key}/{key}')
            if isinstance(d, np.ndarray):
                d = d.ravel()
                md5_hash.update(hashlib.md5(d).digest())
            else:
                if only_arrays:
                    continue
                else:
                    md5_hash.update(repr(d).encode())
    
    return md5_hash.hexdigest()

def gen_entity_ids(
    atoms_per_mol, num_mol, starting_idx=0, add_to=None
):
    """Generates entity ids for a single species.

    Note that all of the atoms in each molecule must occur in the same order and
    be grouped together.

    Parameters
    ----------
    atoms_per_mol : :obj:`int`
        Number of atoms in the molecule.
    num_mol : :obj:`int`
        Number of molecules of this type in the system.
    starting_idx : :obj:`int`
        Number to start entity_id labels.
    add_to : :obj:`list`
        Entity ids to append new ids to.
    
    Returns
    -------
    :obj:`numpy.ndarray`
        Entity ids for a structure.
    """
    entity_ids = []
    for i in range(starting_idx, num_mol+starting_idx):
        entity_ids.extend([i for _ in range(0, atoms_per_mol)])
    
    if add_to is not None:
        if isinstance(add_to, np.ndarray):
            add_to = add_to.tolist()
        return np.array(add_to + entity_ids)
    else:
        return np.array(entity_ids)

def gen_comp_ids(label, num_mol, entity_ids, add_to=None):
    """Prepares the list of component ids for a system with only one species.

    Parameters
    ----------
    label : :obj:`int`
        Species label.
    num_mol : :obj:`int`
        Number of molecules of this type in the system.
    entity_ids : :obj:`int`
        Entity ids for these species.
    add_to : :obj:`list`
        Component ids to append new ids to.

    Returns
    -------
    :obj:`numpy.ndarray`
        Component ids for a structure.
    """
    comp_ids = [label for _ in range(0, num_mol)]
    if add_to is not None:
        if isinstance(add_to, np.ndarray):
            add_to = add_to.tolist()
        return np.array(add_to + comp_ids)
    else:
        return np.array(comp_ids)

def center_structures(Z, R):
    """Centers each structure's center of mass to the origin.

    Previously centered structures should not be affected by this technique.

    Parameters
    ----------
    Z : :obj:`numpy.ndarray`
        Atomic numbers of the atoms in every structure.
    R : :obj:`numpy.ndarray`
        Cartesian atomic coordinates of data set structures.
    
    Returns
    -------
    :obj:`numpy.ndarray`
        Centered Cartesian atomic coordinates.
    """
    # Masses of each atom in the same shape of R.
    if R.ndim == 2:
        R = np.array([R])
    
    masses = np.empty(R[0].shape)
    
    for i in range(len(masses)):
        masses[i,:] = z_to_mass[Z[i]]
    
    for i in range(len(R)):
        r = R[i]
        cm_r = np.average(r, axis=0, weights=masses)
        R[i] = r - cm_r
    
    if R.shape[0] == 1:
        return R[0]
    else:
        return R

def combine_dicts(dict1, dict2):
    """Combine two dictionaries.
    
    Parameters
    ----------
    dict1 : :obj:`dict`
    dict2 : :obj:`dict`

    Returns
    -------
    :obj:`dict`
    """
    for k, v in dict2.items():
        if isinstance(v, collections.abc.Mapping):
            dict1[k] = combine_dicts(dict1.get(k, {}), v)
        else:
            dict1[k] = v
    return dict1
