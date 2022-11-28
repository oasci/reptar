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

from ..extractors import extractorXTB
import numpy as np
from .parser import parser
from ..utils import atoms_by_number, parse_stringfile


class parserXTB(parser):
    """Custom parser for xtb output files."""

    def __init__(self, out_path=None, geom_path=None, traj_path=None, extractors=None):
        """
        Parameters
        ----------
        out_path : :obj:`str`
            Path to the main log file generated by the package.
        geom_path : :obj:`str`, default: ``None``
            Path to a file containing a single geometry.
        traj_path : :obj:`str`, default: ``None``
            Path to a trajectory file from a geometry optimization, MD
            simulation, etc.
        extractors : :obj:`list`, default: ``None``
            Additional extractors for the parser to use.
        """
        self.package = "xtb"
        if extractors is None:
            extractors = []
        extractors.insert(0, extractorXTB())
        super().__init__(out_path, extractors)

        self.geom_path = geom_path
        self.traj_path = traj_path
        if (traj_path is None) and (geom_path is None):
            raise ValueError("geom_path and traj_path cannot both be None")

        self.parsed_info["runtime_info"]["prov"] = "xTB"

    def parse(self):
        """Parses trajectory file and extracts information."""
        # Extract information.
        self.extract_data_out()

        # Adding any structure information.
        Z = []
        comments = []
        R = []
        if self.geom_path is not None:
            Z_geom, _, R_geom = parse_stringfile(self.geom_path)
            Z.extend(Z_geom)
            R.extend(R_geom)
        if self.traj_path is not None:
            Z_traj, comments, R_traj = parse_stringfile(self.traj_path)
            Z.extend(Z_traj)
            R.extend(R_traj)

            # Adds potential energies from trajectory file for MD simulations.
            if self.parsed_info["runtime_info"]["calc_driver"] == "molecular_dynamics":
                energy_pot = []
                # grad_norm = []
                for comment in comments:
                    comment = comment.split()
                    energy_pot.append(float(comment[1]))
                    # grad_norm.append(float(comment[3]))
                if "energy_pot" not in self.parsed_info["outputs"].keys():
                    self.parsed_info["outputs"]["energy_pot"] = []
                self.parsed_info["outputs"]["energy_pot"].extend(energy_pot)

        if len(set(tuple(i) for i in Z)) == 1:
            Z = Z[0]
        else:
            raise ValueError("Atomic numbers are not consistent.")
        Z = np.array(atoms_by_number(Z))
        R = np.array(R)
        if R.ndim == 2:
            R = np.array([R])
        assert R.ndim == 3

        self.parsed_info["system_info"]["atomic_numbers"] = Z
        self.parsed_info["system_info"]["geometry"] = R

        self.after_parse()
        return self.parsed_info

    def after_parse(self):
        """Checks to perform after parsing output file."""
        # xtb prints the last energy twice during optimizations.
        # The last printed energy has more significant figures, so we will
        # get rid of the second to last one.
        # It is also unclear why the structure considered "CYCLE   1" is
        # missing from the geometry log. It goes from the provided structure
        # to the "CYCLE   2" structure. So we will also remove this energy.
        # Note that we overwrite this data later, but it is better to ensure
        # consistency than punting to later in the code.
        if self.parsed_info["runtime_info"]["calc_driver"] == "optimization":
            del self.parsed_info["outputs"]["energy_scf"][1]
            del self.parsed_info["outputs"]["energy_scf"][-2]

        if "success" not in self.parsed_info["runtime_info"].keys():
            self.parsed_info["runtime_info"]["success"] = False
