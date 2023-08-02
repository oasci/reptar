# MIT License
#
# Copyright (c) 2022-2023, Alex M. Maldonado
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

from __future__ import annotations
import os
from collections.abc import Iterable
import numpy as np
import qcelemental as qcel
from .cube import get_R_span, get_max_grid_points
from ..utils import common_elements
from ..logger import ReptarLogger

log = ReptarLogger(__name__)


class MissingDataException(Exception):
    r"""Data that should be provided for specific task."""

    def __init__(self, message):
        super().__init__(message)
        self.message = message


# pylint: disable=invalid-name
class Data:
    r"""Automates handling data with reptar calculations. This can be used in a couple
    of ways.

    - **Write data.** Providing ``rfile`` and keys allows you to use the
      :meth:`~reptar.calculators.Data.save` method to write stored arrays to the reptar
      file.
    - **Worker data.** Used to store data computed by a worker that is then used with
      :meth:`~reptar.calculators.Data.update` and ``idxs_source`` to update another
      :class:`~reptar.calculators.Data` object.

    We store data in this class with properties (e.g., ``E``, ``R_opt``, ``cube_V``)
    that have a companion property with a ``_key`` suffix. These ``_key``
    properties are used to write data to ``rfile``.
    """

    def __init__(
        self,
        rfile: "File" | None = None,
        idxs_source: np.ndarray | None = None,
    ) -> None:
        r"""
        Parameters
        ----------
        rfile
            A file to save data to if desired.
        idxs_source
            Specifies slice indices of sources where values will eventually be updated.
            This is normally specified in calculation workers and used in
            :meth:`~reptar.calculators.Data.update`.
        """
        self.rfile = rfile

        self.idxs_source = idxs_source
        self.data_labels = ["Z", "R", "E", "G", "conv_opt", "R_opt", "cube_R", "cube_V"]
        self.data_key_labels = [key + "_key" for key in self.data_labels]

        self.task_data_map = {
            "E": ["E"],
            "G": ["E", "G"],
            "opt": ["E", "conv_opt", "R_opt"],
            "cube": ["cube_R", "cube_V"],
        }
        self._require_init_cube = False

    def update(self, data: Data) -> None:
        r"""Update self with another data object. ``rfile`` must be specified.

        Parameters
        ----------
        data
            A subset of data that will be used to update ``self``.
        """
        if data.idxs_source is None:
            raise ValueError("data must specify `idxs_source`")
        for data_label in self.data_labels:
            data_attr = getattr(data, data_label)
            if data_attr is not None:
                self.add_subset(data_label, data.idxs_source, data_attr)

    def add_subset(self, prop: str, idxs: np.ndarray, values: np.ndarray) -> None:
        r"""Add a subset of values to a property"""
        data_original = getattr(self, prop)
        data_original[idxs] = values
        setattr(self, prop, data_original)

    def get_idxs_todo(
        self,
        tasks: Iterable[str],
        start_slice: int | None = None,
        end_slice: int | None = None,
    ) -> np.ndarray:
        r"""Determines indices of missing data depending on what is loaded and desired
        tasks.

        Parameters
        ----------
        tasks
            Desired computations to perform.
        start_slice
            Slice arrays in ``args`` starting at this index.
        end_slice
            Slice arrays` in ``args`` stopping at this index.
        """
        todo = {}
        if isinstance(self.E, np.ndarray):
            todo["E"] = np.argwhere(np.isnan(self.E[start_slice:end_slice]))[:, 0]
        if isinstance(self.G, np.ndarray):
            todo["G"] = np.argwhere(np.isnan(self.G[start_slice:end_slice][:, 0]))[:, 0]
        if isinstance(self.conv_opt, np.ndarray):
            todo["opt"] = np.where(~self.conv_opt[start_slice:end_slice])[0]
        if isinstance(self.cube_V, np.ndarray):
            todo["cube"] = np.argwhere(
                np.isnan(self.cube_V[start_slice:end_slice][:, 0])
            )[:, 0]

        # We make a list of indices that are missing all requested calculations.
        idxs_todo = todo[tasks[0]]
        if len(tasks) > 1:
            for i in range(1, len(tasks)):
                idxs_todo = common_elements(idxs_todo, todo[tasks[i]])
        return idxs_todo

    def validate(self, tasks: Iterable[str] | str | None) -> None:
        r"""Will check if required data is loaded to complete tasks.

        Parameters
        ----------
        tasks
            Calculations that will be performed with these data. If ``None``, then this
            only checks for ``Z`` and ``R``.

        Raises
        ------
        ``MissingDataException``
            If data that should be present for specific task is missing.
        """
        if self.Z is None:
            raise MissingDataException("Z is required for all tasks but is missing")
        if self.R is None:
            raise MissingDataException("R is required for all tasks but is missing")

        if tasks is None:
            return None

        if isinstance(tasks, str):
            tasks = (tasks,)

        for task in tasks:
            data_attrs = self.task_data_map[task]
            for data_attr in data_attrs:
                if getattr(self, data_attr) is None:
                    raise MissingDataException(
                        f"{data_attr} is required for {task} but is missing"
                    )
        return None

    def initialize_array(self, label: str) -> np.ndarray:
        r"""Initialize array given ``Z`` and ``R``.

        Properties
        ----------
        label
            Name of the attribute.
        """
        if "cube" in label:
            self._require_init_cube = True
        setattr(self, label, np.full(**self.array_init_specs[label]))

    def save(self) -> None:
        r"""Write all data to the reptar file."""
        for data_key_label in self.data_key_labels:
            data_key = getattr(self, data_key_label)
            data = getattr(self, data_key_label[:-4])
            if (data_key is not None) and (data is not None):
                self.rfile.put(data_key, data)

    @property
    def array_init_specs(self) -> dict[str, "Any"]:
        r"""Specifications for initializing arrays."""
        self.validate(None)
        array_specs = {
            "E": {
                "shape": self.R.shape[0],
                "fill_value": np.nan,
                "dtype": np.float64,
            },
            "G": {"shape": self.R.shape, "fill_value": np.nan, "dtype": np.float64},
            "conv_opt": {"shape": self.R.shape[0], "fill_value": False, "dtype": bool},
            "R_opt": {"shape": self.R.shape, "fill_value": np.nan, "dtype": np.float64},
            "cube_R": {
                "shape": self.R.shape,
                "fill_value": np.nan,
                "dtype": np.float64,
            },
            "cube_V": {
                "shape": self.R.shape,
                "fill_value": np.nan,
                "dtype": np.float64,
            },
        }
        if self._require_init_cube:
            if self._max_grid_points is None:
                self._max_grid_points = get_max_grid_points(
                    get_R_span(self.R), self.cube_grid_overage, self.cube_grid_spacing
                )
            array_specs["cube_R"] = {
                "shape": (self.R.shape[0], self._max_grid_points, 3),
                "fill_value": np.nan,
                "dtype": np.float64,
            }
            array_specs["cube_V"] = {
                "shape": (self.R.shape[0], self._max_grid_points),
                "fill_value": np.nan,
                "dtype": np.float64,
            }

        return array_specs

    def prepare_tasks(
        self,
        tasks: Iterable[str],
        source_key: str,
        source_labels: dict[str, str],
        dest_key: str,
        dest_labels: dict[str, str],
    ) -> None:
        r"""Prepare group and data for reptar calculations

        Parameters
        ----------
        tasks
            Reptar calculations that will be ran with
            :class:`~reptar.calculators.Driver`.
        source_key
            Key to group containing data sources.
        source_labels
            Labels of data from group ``source_key`` to populate the data object. No
            ``_key`` properties in :class:`~reptar.calculators.Data` will be added from
            this ``dict``.
        dest_key
            Key to store calculation results.
        dest_labels
            Labels of data in group ``dest_key`` to populate the data object. Any data
            that is not provided will be initialized here based on ``tasks``.
            ``_key`` properties in :class:`~reptar.calculators.Data` will be added from
            this ``dict``.

        Returns
        -------
        :obj:`reptar.calculators.Data`
            Data for calculations.
        """
        log.debug("Retrieving source data")
        for data_attr, label in source_labels.items():
            data_key = os.path.join(source_key, label)
            value = self.rfile.get(data_key)
            setattr(self, data_attr, value[:])
        self.validate(None)  # Checks for Z and R

        log.debug("Checking calculation data")
        try:
            self.rfile.get(dest_key)
        except RuntimeError as e:
            if " does not exist" in str(e):
                self.rfile.create_group(dest_key)
            else:
                raise RuntimeError from e

        for data_attr, label in dest_labels.items():
            data_key = os.path.join(dest_key, label)
            setattr(self, data_attr + "_key", data_key)

            # If included from source, then we may have previous data already.
            # Check before initializing.
            if getattr(self, data_attr) is None:
                self.initialize_array(data_attr)

        self.validate(tasks)
        self.save()

    @property
    def Z(self) -> np.ndarray | None:
        r"""Atomic numbers"""
        if hasattr(self, "_Z"):
            return self._Z
        return None

    @Z.setter
    def Z(self, value: np.ndarray | None) -> None:
        # We directly store the value if it is None or an array.
        if isinstance(value, type(None)):
            self._Z = value
        elif isinstance(value, np.ndarray):
            assert value.ndim == 1
            self._Z = value

    @property
    def Z_key(self) -> str | None:
        r"""Atomic number key to save in reptar file."""
        if hasattr(self, "_Z_key"):
            return self._Z_key
        return None

    @Z_key.setter
    def Z_key(self, value: str):
        self._Z_key = value

    @property
    def R(self) -> np.ndarray | None:
        r"""Last geometry during an optimization."""
        if hasattr(self, "_R"):
            return self._R
        return None

    @R.setter
    def R(self, value: Iterable[np.ndarray, np.ndarray] | np.ndarray | None) -> None:
        # We directly store the value if it is None or an array.
        if isinstance(value, type(None)):
            self._R = value
        elif isinstance(value, np.ndarray):
            assert value.ndim == 3
            self._R = value
        self._max_grid_points = None  # Force recalculate if requested.

    @property
    def R_key(self) -> str | None:
        r"""Atomic number key to save in reptar file."""
        if hasattr(self, "_R_key"):
            return self._R_key
        return None

    @R_key.setter
    def R_key(self, value: str):
        self._R_key = value

    @property
    def E(self) -> np.ndarray | None:
        r"""Energy."""
        if hasattr(self, "_E"):
            return self._E
        return None

    @E.setter
    def E(self, value: Iterable[np.ndarray, np.ndarray] | np.ndarray | None) -> None:
        # We directly store the value if it is None or an array.
        if isinstance(value, type(None)):
            self._E = value
        elif isinstance(value, np.ndarray):
            assert value.ndim == 1
            self._E = value
        # Next case is if we provided a tuple of (idxs, data). This happens in
        # self.update().
        elif isinstance(value, (tuple, list)):
            if len(value) != 2:
                raise RuntimeError(
                    "If updating this data object, you must provide idxs and values"
                )
            idxs, e_values = value
            self._E[idxs] = e_values

    @property
    def E_key(self) -> str | None:
        r"""Energy key to save in reptar file."""
        if hasattr(self, "_E_key"):
            return self._E_key
        return None

    @E_key.setter
    def E_key(self, value: str):
        self._E_key = value

    @property
    def G(self) -> np.ndarray | None:
        r"""Atomic gradients."""
        if hasattr(self, "_G"):
            return self._G
        return None

    @G.setter
    def G(self, value: Iterable[np.ndarray, np.ndarray] | np.ndarray | None):
        # We directly store the value if it is None or an array.
        if isinstance(value, type(None)):
            self._G = value
        elif isinstance(value, np.ndarray):
            assert value.ndim == 3
            self._G = value

    @property
    def G_key(self) -> str | None:
        r"""Gradient key to save in reptar file."""
        if hasattr(self, "_G_key"):
            return self._G_key
        return None

    @G_key.setter
    def G_key(self, value: str):
        self._G_key = value

    @property
    def conv_opt(self) -> np.ndarray | None:
        r"""Convergence of geometry optimizations."""
        if hasattr(self, "_conv_opt"):
            return self._conv_opt
        return None

    @conv_opt.setter
    def conv_opt(
        self, value: Iterable[np.ndarray, np.ndarray] | np.ndarray | None
    ) -> None:
        # We directly store the value if it is None or an array.
        if isinstance(value, type(None)):
            self._conv_opt = value
        elif isinstance(value, np.ndarray):
            assert value.ndim == 1
            self._conv_opt = value

    @property
    def conv_opt_key(self) -> str | None:
        r"""Optimization convergence key to save in reptar file."""
        if hasattr(self, "_conv_opt_key"):
            return self._conv_opt_key
        return None

    @conv_opt_key.setter
    def conv_opt_key(self, value: str):
        self._conv_opt_key = value

    @property
    def R_opt(self) -> np.ndarray | None:
        r"""Last geometry during an optimization."""
        if hasattr(self, "_R_opt"):
            return self._R_opt
        return None

    @R_opt.setter
    def R_opt(
        self, value: Iterable[np.ndarray, np.ndarray] | np.ndarray | None
    ) -> None:
        # We directly store the value if it is None or an array.
        if isinstance(value, type(None)):
            self._R_opt = value
        elif isinstance(value, np.ndarray):
            assert value.ndim == 3
            self._R_opt = value

    @property
    def R_opt_key(self) -> str | None:
        r"""Key for last geometry optimization step to save in reptar file."""
        if hasattr(self, "_R_opt_key"):
            return self._R_opt_key
        return None

    @R_opt_key.setter
    def R_opt_key(self, value: str):
        self._R_opt_key = value

    @property
    def cube_R(self) -> np.ndarray | None:
        r"""Cartesian coordinates of points where a property is probed on a grid."""
        if hasattr(self, "_cube_R"):
            return self._cube_R
        return None

    @cube_R.setter
    def cube_R(self, value: Iterable[np.ndarray, np.ndarray] | np.ndarray | None):
        # We directly store the value if it is None or an array.
        if isinstance(value, type(None)):
            self._cube_R = value
        elif isinstance(value, np.ndarray):
            assert value.ndim == 3
            self._cube_R = value

    @property
    def cube_R_key(self) -> str | None:
        r"""Grid/cube coordinate key to save in reptar file."""
        if hasattr(self, "_cube_R_key"):
            return self._cube_R_key
        return None

    @cube_R_key.setter
    def cube_R_key(self, value: str):
        self._cube_R_key = value

    @property
    def cube_V(self) -> np.ndarray | None:
        r"""Property values on a grid specified by ``cube_R``."""
        if hasattr(self, "_cube_V"):
            return self._cube_V
        return None

    @cube_V.setter
    def cube_V(self, value: Iterable[np.ndarray, np.ndarray] | np.ndarray | None):
        # We directly store the value if it is None or an array.
        if isinstance(value, type(None)):
            self._cube_V = value
        elif isinstance(value, np.ndarray):
            assert value.ndim == 2
            self._cube_V = value

    @property
    def cube_V_key(self) -> str | None:
        r"""Grid/cube coordinate key to save in reptar file."""
        if hasattr(self, "_cube_V_key"):
            return self._cube_V_key
        return None

    @cube_V_key.setter
    def cube_V_key(self, value: str):
        self._cube_V_key = value

    @property
    def cube_grid_overage(self) -> np.ndarray | None:
        r"""Overage for cube in Angstroms.

        `Defaults <https://psicode.org/psi4manual/master/
        autodoc_glossary_options_c.html#term-CUBIC_GRID_OVERAGE-GLOBALS>`__ to ``4.0``
        Bohr (approximately ``2.117`` Angstroms) in ``x``, ``y``, and ``z``.
        """
        if hasattr(self, "_cube_grid_overage"):
            return self._cube_grid_overage
        default_overage = np.array([4.0, 4.0, 4.0], dtype=np.float64)
        default_overage *= qcel.constants.bohr2angstroms
        return default_overage

    @cube_grid_overage.setter
    def cube_grid_overage(self, value: np.ndarray | float) -> None:
        if isinstance(value, float):
            value = np.array([value, value, value], dtype=np.float64)
        self._cube_grid_overage = value

    @property
    def cube_grid_spacing(self) -> np.ndarray | None:
        r"""Spacing for cube grid in Angstroms.

        `Defaults <https://psicode.org/psi4manual/master/
        autodoc_glossary_options_c.html#term-CUBIC_GRID_SPACING-GLOBALS>`__ to ``0.2``
        Bohr (approximately ``0.108`` Angstroms) in ``x``, ``y``, and ``z``.
        """
        if hasattr(self, "_cube_grid_spacing"):
            return self._cube_grid_spacing
        default_spacing = np.array([0.2, 0.2, 0.2], dtype=np.float64)
        default_spacing *= qcel.constants.bohr2angstroms
        return default_spacing

    @cube_grid_spacing.setter
    def cube_grid_spacing(self, value: np.ndarray | float) -> None:
        if isinstance(value, float):
            value = np.array([value, value, value], dtype=np.float64)
        self._cube_grid_spacing = value
