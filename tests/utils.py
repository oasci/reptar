# MIT License
#
# Copyright (c) 2022-2023, OASCI
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

"""Utilities for tests."""

# pylint: skip-file

import numpy as np
import scipy


def _pdist(r):
    r"""Compute pairwise Euclidean distance matrix between all atoms.

    Parameters
    ----------
    r : :obj:`numpy.ndarray`
        Array of size 3N containing the Cartesian coordinates of
        each atom.

    Returns
    -------
    :obj:`numpy.ndarray`
        Array of size N x N containing all pairwise distances between atoms.
    """

    r = r.reshape(-1, 3)
    n_atoms = r.shape[0]

    pdist = scipy.spatial.distance.pdist(r, "euclidean")
    tril_idxs = np.tril_indices(n_atoms, k=-1)
    return scipy.spatial.distance.squareform(pdist, checks=False)[tril_idxs]
