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

"""Analyze reptars from MD simulations."""
# TODO: Change this script into a function.

"""
import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rdf
import matplotlib.pyplot as plt



pdb_path = './test.pdb'
# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Approximate volume of sphere as box.
sphere_raidus = 12.5000030  # Angstroms
sphere_volume = (4/3)*np.pi*(sphere_raidus)**3
box_length = np.cbrt(sphere_volume)

u = mda.Universe(
    pdb_path, dt=0.005
)
u.dimensions=[box_length, box_length, box_length, 90, 90, 90]


water = u.select_atoms('resname WAT')
methanol = u.select_atoms('resname MET')

rdf_calc_all_sites = rdf.InterRDF(
    water, methanol,
    density=False,
    nbins=200,
    range=(0, 10),
    verbose=True
)
rdf_calc_all_sites.run()

plt.plot(rdf_calc_all_sites.results.bins, rdf_calc_all_sites.results.rdf)
plt.savefig('test.png', dpi=600)
plt.clf()


water_o = u.select_atoms('resname WAT and name O')
methanol_o = u.select_atoms('resname MET and name O')

print(np.sqrt(np.sum((water_o.positions-methanol_o.positions)**2,axis=1)))
exit()
rdf_calc = rdf.InterRDF_s(
    u,
    [[water_o, methanol_o]],
    density=False,
    nbins=100,
    range=(0, 10),
    verbose=True
)
rdf_calc.run()
print(rdf_calc.results.rdf[0].shape)
plt.plot(rdf_calc.results.bins, rdf_calc.results.rdf[0][0, 0])
plt.savefig('test-site.png', dpi=600)
"""