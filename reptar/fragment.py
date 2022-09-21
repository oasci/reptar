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

"""Groups atoms into fragments."""

import numpy as np
from scipy.spatial.distance import cdist

def pairwise_distance(R, cutoff=1.5, num_frags=None):
    """Group atoms together into fragments based on some cutoff distance.

    This is not very dependable.

    Parameters
    ---------
    R : :obj:`numpy.ndarray`, ndim: ``2``
        Cartesian coordinates of all atoms in the system in Angstroms.
    cutoff : :obj:`float`, default: ``1.5``
        Atoms within this distance is group together in a fragment.
    num_frags : :obj:`int`, default: ``None``
        Known number of fragments in structure. No check is performed if
        ``None``.
    
    Returns
    -------
    :obj:`numpy.ndarray`
        Computed ``entity_ids`` of the structure.
    """
    if R.ndim != 2:
        raise ValueError('Array must have two dimensions.')

    entity_ids = np.empty((R.shape[0],))
    entity_ids[:] = np.NaN

    all_distances = cdist(R, R, 'euclidean')

    frag_idx = 0
    for i in range(len(entity_ids)):
        frag_atoms = np.argwhere(all_distances[i]<cutoff).T[0]
        # If one of the atoms are already in fragment we include these atoms.
        atoms_entity_ids = entity_ids[frag_atoms]
        if np.isnan(atoms_entity_ids).all():
            entity_ids[frag_atoms] = frag_idx
            frag_idx += 1
        else:
            prev_frag_idx = np.unique(
                atoms_entity_ids[~np.isnan(atoms_entity_ids)]
            )
            assert len(prev_frag_idx) == 1
            prev_frag_idx = prev_frag_idx[0]
            entity_ids[frag_atoms] = prev_frag_idx
    
    if num_frags is not None:
        if num_frags != frag_idx-1:
            raise ValueError(
                f'Detected {frag_idx-1} fragments instead of {num_frags}.'
            )
    
    return entity_ids

