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

"""Implementations of structural criteria."""

import itertools
import numpy as np
from .utils import z_to_mass

def get_center_of_mass(Z, R):
    """Compute the center of mass.

    Parameters
    ----------
    Z : :obj:`numpy.ndarray`, ndim: ``1``
        Atomic numbers of all atoms in the system.
    R : :obj:`numpy.ndarray`, ndim: ``3``
        Cartesian coordinates.
    
    Returns
    -------
    :obj:`float`
        The center of mass cartesian coordinates.
    """
    if R.ndim == 2:
        R = np.array([R])
    masses = np.empty(R[0].shape)
    for i in range(len(masses)):
        masses[i,:] = z_to_mass[Z[i]]
    R_masses = np.full(R.shape, masses)
    cm_structure = np.average(R, axis=1, weights=R_masses)
    return cm_structure

def criteria(desc, desc_args, cutoff):
    """Evaluates some descriptor criteria.

    Parameters
    ----------
    desc : ``callable``
        A descriptor function.
    desc_args : :obj:`tuple`, ndim: ``1``
        All of the descriptor function arguments.
    cutoff : :obj:`float`
        Descriptor cutoff.
    
    Returns
    -------
    :obj:`numpy.ndarray`, ndim: ``1``
        Descriptor values.
    :obj:`numpy.ndarray`, ndim: ``1``
        If the descriptor value is less than or equal to the cutoff.
    """
    v_desc = desc(*desc_args)
    accept = v_desc <= cutoff
    return v_desc, accept

def max_atom_pair_dist(Z, R):
    """The largest atomic pairwise distance.

    Parameters
    ----------
    Z : :obj:`numpy.ndarray`, ndim: ``1``
        Atomic numbers of all atoms in the system.
    R : :obj:`numpy.ndarray`, ndim: ``3``
        Cartesian coordinates.

    Returns
    -------
    :obj:`float`
        The largest pairwise distance.
    """
    if R.ndim == 2:
        R = np.array([R])
    # Finds all atom pairs.
    all_pairs = np.array(list(itertools.combinations(range(len(Z)), 2)))
    # Creates arrays of all the points for all structures.
    pair0_points = np.array([[R[i,j] for j in all_pairs[:,0]] for i in range(len(R))])
    pair1_points = np.array([[R[i,j] for j in all_pairs[:,1]] for i in range(len(R))])
    # Computes the distance, then largest distances of all structures.
    distances = np.linalg.norm(pair0_points - pair1_points, axis=2)
    max_distances = np.amax(distances, axis=1)
    return max_distances

def com_distance_sum(Z, R, entity_ids):
    """The sum of pairwise distances from each entity's center of mass to
    the total structure center of mass.

    Parameters
    ----------
    Z : :obj:`numpy.ndarray`, ndim: ``1``
        Atomic numbers of all atoms in the system.
    R : :obj:`numpy.ndarray`, ndim: ``3``
        Cartesian coordinates.
    entity_ids : :obj:`numpy.ndarray`, ndim: ``1``
        A uniquely identifying integer specifying what atoms belong to
        which entities. Entities can be a related set of atoms, molecules,
        or functional group. For example, a water and methanol molecule
        could be ``[0, 0, 0, 1, 1, 1, 1, 1, 1]``.

    Returns
    -------
    :obj:`float`
        Calculated distance metric.
    """
    if R.ndim == 2:
        R = np.array([R])
    cm_structures = get_center_of_mass(Z, R)

    d_sum = np.zeros(R.shape[0])
    for entity_id in set(entity_ids):
        atom_idxs = np.where(entity_ids == entity_id)[0]
        cm_entity = get_center_of_mass(Z, R[:,atom_idxs])
        d_sum += np.linalg.norm(cm_structures - cm_entity, axis=1)
    
    return d_sum
