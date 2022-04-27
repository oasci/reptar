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

"""Computations using Psi4"""

import numpy as np
import psi4
import exdir

def psi4_pes(group, e_label, g_label, psi4_options, psi4_method):
    Z = grp['atomic_numbers']
    R = grp['geometry']
    E = np.array([np.nan for _ in range(len(R))])
    E = group.require_dataset(e_label, data=E)
    G = np.empty(R.shape)
    G[:] = np.nan
    G = group.require_dataset(g_label, data=G)

    psi4.set_options(psi4_options)

    for i in range(len(R)):
        mol = psi4.core.Molecule.from_arrays(geom=R[i], elez=Z)
        G, wfn_mp2 = psi4.gradient(
            psi4_method, molecule=mol, return_wfn=True
        )
        G[i] = G/psi4.constants.bohr2angstroms
        E[i] = wfn_mp2.energy()
