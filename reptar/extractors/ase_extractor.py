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

from math import sqrt
import numpy as np
import qcelemental as qcel
from .extractor import Extractor

try:
    from ase.calculators.calculator import PropertyNotImplementedError
except ImportError:
    pass


class ExtractorASE(Extractor):
    r"""ASE extractor for the Atoms object."""
    # We always pass the file object into methods here.
    # pylint: disable=unused-argument

    def __init__(self):
        super().__init__()
        self._e_conv = qcel.constants.conversion_factor("eV", "hartree")
        self._velcs_conv = (
            1e-15
            * 1e10
            * sqrt(
                qcel.constants.elementary_charge  # pylint: disable=no-member
                / qcel.constants.unified_atomic_mass_unit  # pylint: disable=no-member
            )
        )
        # pylint: disable-next=invalid-name
        self._kB = qcel.constants.kb * qcel.constants.conversion_factor("joule", "eV")

    @property
    def triggers(self):
        r"""Activates for every Atoms object."""
        trig = ((lambda line: True, "atoms"),)
        return trig

    def atoms(self, traj, atoms):
        r"""Parse properties from Atoms object.

        Parameters
        ----------
        traj : ``ase.io.trajectory.Trajectory``
            Loaded ASE trajectory file.
        atoms : ``ase.Atoms``
            ASE Atoms object with at least atomic numbers and positions.
        """
        if "prov_version" not in self.parsed_info["runtime_info"]:
            self.parsed_info["runtime_info"]["prov_version"] = traj.ase_version

        self._atomic_numbers(atoms)
        self._geometry(atoms)
        self._periodic(atoms)
        self._energy_pot(atoms)
        self._energy_ke(atoms)
        self._forces(atoms)
        self._velcs(atoms)

    def _atomic_numbers(self, atoms):
        r"""Parse atomic numbers using ``get_atomic_numbers()``.

        Parameters
        ----------
        atoms : ``ase.Atoms``
            ASE Atoms object.
        """
        atomic_numbers = atoms.get_atomic_numbers()
        if atomic_numbers.shape == (0,):
            return
        if "atomic_numbers" in self.parsed_info["system_info"]:
            assert np.all(
                self.parsed_info["system_info"]["atomic_numbers"] == atomic_numbers
            )
        else:
            self.parsed_info["system_info"]["atomic_numbers"] = atomic_numbers

    def _geometry(self, atoms):
        r"""Parse geometry using ``get_positions()``.

        Parameters
        ----------
        atoms : ``ase.Atoms``
            ASE Atoms object.
        """
        geometry = atoms.get_positions()
        if geometry.shape == (0, 3):
            return
        if "geometry" in self.parsed_info["system_info"]:
            self.parsed_info["system_info"]["geometry"].append(geometry)
        else:
            self.parsed_info["system_info"]["geometry"] = [geometry]

    def _periodic(self, atoms):
        r"""Parse periodic cell information using ``pbc`` and ``get_cell()``.

        Parameters
        ----------
        atoms : ``ase.Atoms``
            ASE Atoms object.
        """
        periodic = atoms.pbc
        if np.any(periodic):
            periodic_cell = atoms.get_cell()[:]
            if "periodic" in self.parsed_info["system_info"]:
                assert np.all(self.parsed_info["system_info"]["periodic"] == periodic)
                assert np.all(
                    self.parsed_info["system_info"]["periodic_cell"] == periodic_cell
                )
            else:
                self.parsed_info["system_info"]["periodic"] = periodic
                self.parsed_info["system_info"]["periodic_cell"] = periodic_cell

    def _energy_pot(self, atoms):
        r"""Parse total potential energy using ``get_potential_energy()``.

        Parameters
        ----------
        atoms : ``ase.Atoms``
            ASE Atoms object.
        """
        try:
            energy_pot = atoms.get_potential_energy()  # eV
            energy_pot *= self._e_conv  # Eh
            if "energy_pot" in self.parsed_info["outputs"]:
                self.parsed_info["outputs"]["energy_pot"].append(energy_pot)
            else:
                self.parsed_info["outputs"]["energy_pot"] = [energy_pot]
        except PropertyNotImplementedError:
            pass
        except RuntimeError:
            pass

    def _energy_ke(self, atoms):
        r"""Parse total kinetic energy using ``get_kinetic_energy()``.

        The temperature is calculated using the expression
        :math:`T = \frac{2}{3} \frac{E_{kinetic}}{N_{atoms}}`.

        Parameters
        ----------
        atoms : ``ase.Atoms``
            ASE Atoms object.
        """
        energy_ke = atoms.get_kinetic_energy()  # eV

        if energy_ke == 0.0:
            return

        # Get temperature before converting.
        n_atoms = len(self.parsed_info["system_info"]["atomic_numbers"])
        energy_ke_per_atom = energy_ke / n_atoms
        temp = energy_ke_per_atom / (1.5 * self._kB)

        energy_ke *= self._e_conv  # Eh
        if "energy_ke" in self.parsed_info["outputs"]:
            self.parsed_info["outputs"]["energy_ke"].append(energy_ke)
            self.parsed_info["outputs"]["temp"].append(temp)
        else:
            self.parsed_info["outputs"]["energy_ke"] = [energy_ke]
            self.parsed_info["outputs"]["temp"] = [temp]

    def _forces(self, atoms):
        r"""Parse atomic forces using ``get_forces()``.

        Parameters
        ----------
        atoms : ``ase.Atoms``
            ASE Atoms object.
        """
        try:
            forces = atoms.get_forces()  # eV/Ang
            forces *= self._e_conv  # Eh/Ang
            if "forces" in self.parsed_info["outputs"]:
                self.parsed_info["outputs"]["forces"].append(forces)
            else:
                self.parsed_info["outputs"]["forces"] = [forces]
        except PropertyNotImplementedError:
            pass
        except RuntimeError:
            pass

    def _velcs(self, atoms):
        r"""Parse velocities using ``get_velocities()``.

        Parameters
        ----------
        atoms : ``ase.Atoms``
            ASE Atoms object.
        """
        velcs = atoms.get_velocities()  # Ang / (Ang * sqrt(amu / eV))
        if velcs.shape == (0, 3):
            return
        velcs *= self._velcs_conv  # Ang / fs
        if "velcs" in self.parsed_info["outputs"]:
            self.parsed_info["outputs"]["velcs"].append(velcs)
        else:
            self.parsed_info["outputs"]["velcs"] = [velcs]
