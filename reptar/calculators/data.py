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
from typing import Iterable
import numpy as np
from .logger import ReptarLogger

log = ReptarLogger(__name__)

# pylint: disable=invalid-name
class Data:
    r"""Automates handling data with reptar calculations. This can be used in a couple
    of ways.

    - **Write data.** Providing ``rfile`` and keys allows you to use the
      :meth:`~reptar.calculators.Data.save` method to write stored arrays to the reptar file.
    - **Worker data.** Used to store data computed by a worker that is then used with
      :meth:`~reptar.calculators.Data.update` and ``idxs_source`` to update another
      :class:`~reptar.calculators.Data` object.

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

    def update(self, data: Data) -> None:
        r"""Update self with another data object. ``rfile`` must be specified.

        Parameters
        ----------
        data
            A subset of data that will be used to update ``self``.
        """
        if data.idxs_source is None:
            raise ValueError("data must specify `idxs_source`")
        for data_attr_label in ["E", "G", "cube_R", "cube_V"]:
            data_attr = getattr(data, data_attr_label)
            if data_attr is not None:
                setattr(self, data_attr_label, (data.idxs_source, data_attr))

    def save(self) -> None:
        r"""Will write all data to the reptar file."""
        for data_key_label in ["E_key", "G_key", "cube_R_key", "cube_V_key"]:
            if data_key_label is not None:
                self.rfile.put(data_key_label, getattr(self, data_key_label[:-4]))

    @property
    def E(self) -> np.ndarray | None:
        r"""Energy."""
        if hasattr(self, "_E"):
            return self._E
        return None

    @E.setter
    def E(self, value: Iterable[str, np.ndarray] | np.ndarray | None) -> None:
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
    def G(self, value: Iterable[np.ndarray] | np.ndarray | None):
        # We directly store the value if it is None or an array.
        if isinstance(value, type(None)):
            self._G = value
        elif isinstance(value, np.ndarray):
            assert value.ndim == 3
            self._G = value
        # Next case is if we provided a tuple of (idxs, data). This happens in
        # self.update().
        elif isinstance(value, (tuple, list)):
            if len(value) != 2:
                raise RuntimeError(
                    "If updating this data object, you must provide idxs and values"
                )
            idxs, g_values = value
            self._G[idxs] = g_values

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
    def cube_R(self) -> np.ndarray | None:
        r"""Cartesian coordinates of points where a property is probed on a grid."""
        if hasattr(self, "_cube_R"):
            return self._cube_R
        return None

    @cube_R.setter
    def cube_R(self, value: Iterable[np.ndarray] | np.ndarray | None):
        # We directly store the value if it is None or an array.
        if isinstance(value, type(None)):
            self._cube_R = value
        elif isinstance(value, np.ndarray):
            assert value.ndim == 3
            self._cube_R = value
        # Next case is if we provided a tuple of (idxs, data). This happens in
        # self.update().
        elif isinstance(value, (tuple, list)):
            if len(value) != 2:
                raise RuntimeError(
                    "If updating this data object, you must provide idxs and values"
                )
            idxs, cube_r_values = value
            self._cube_R[idxs] = cube_r_values

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
    def cube_V(self, value: Iterable[np.ndarray] | np.ndarray | None):
        # We directly store the value if it is None or an array.
        if isinstance(value, type(None)):
            self._cube_V = value
        elif isinstance(value, np.ndarray):
            assert value.ndim == 3
            self._cube_V = value
        # Next case is if we provided a tuple of (idxs, data). This happens in
        # self.update().
        elif isinstance(value, (tuple, list)):
            if len(value) != 2:
                raise RuntimeError(
                    "If updating this data object, you must provide idxs and values"
                )
            idxs, cube_r_values = value
            self._cube_V[idxs] = cube_r_values

    @property
    def cube_V_key(self) -> str | None:
        r"""Grid/cube coordinate key to save in reptar file."""
        if hasattr(self, "_cube_V_key"):
            return self._cube_V_key
        return None

    @cube_V_key.setter
    def cube_V_key(self, value: str):
        self._cube_V_key = value
