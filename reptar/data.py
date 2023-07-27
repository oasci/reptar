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

import numpy as np
from .logger import ReptarLogger

log = ReptarLogger(__name__)

# pylint: disable=invalid-name
class Data:
    r"""Automates handling data with reptar calculations."""

    def __init__(
        self,
        rfile,
        idxs_parent=None,
        E=None,
        E_key=None,
        G=None,
        G_key=None,
        cube_R=None,
        cube_R_key=None,
        cube_V=None,
        cube_V_key=None,
    ):
        r""" """
        self.rfile = rfile

        self.idxs_parent = idxs_parent

        self.E = E
        self.E_key = E_key
        self.G = G
        self.G_key = G_key
        self.cube_R = cube_R
        self.cube_R_key = cube_R_key
        self.cube_V = cube_V
        self.cube_V_key = cube_V_key

    @property
    def E(self):
        r"""Energy.

        :type:`numpy.ndarray`
        """
        return self._E

    @E.setter
    def E(self, value):
        # We directly store the value if it is None or an array.
        if isinstance(value, type(None)):
            self._E = value
        elif isinstance(value, np.ndarray):
            assert value.ndim == 1
            self._E = value
        # Next case is if we provided a tuple of (data, idxs). This happens in
        # self.update().
        elif isinstance(value, (tuple, list)):
            if len(value) != 2:
                raise RuntimeError(
                    "If updating this data object, you must provide idxs and values"
                )
            idxs, e_values = value
            self._E[idxs] = e_values

    @property
    def E_key(self):
        r"""Energy.

        :type:`numpy.ndarray`
        """
        return self._E_key

    @E_key.setter
    def E_key(self, value):
        assert isinstance(value, str)
        self._E_key = value

    def update(self, data_obj):
        r"""Update self with another data object.

        Parameters
        ----------
        data_obj : :obj:`reptar.Data`
        """
        # TODO

    def save(self):
        r"""Will write all data to the reptar file."""
        # TODO
