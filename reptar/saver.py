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

from .reptar_file import File

class Saver:
    r"""Saves data to File."""

    def __init__(self, rfile_path, data_keys):
        """
        Parameters
        ----------
        rfile_path : :obj:`reptar.File`
            Path to File.
        data_keys : :obj:`tuple`
            Keys to save data for the particular worker
        """
        self.rfile_path = rfile_path
        self.data_keys = data_keys

    def save(self, *data):
        r"""Put data to rfile.

        Parameters
        ----------
        *data
            Data to add to file in the same order as ``data_keys``.

        Returns
        -------
        :obj:`reptar.File`
            File after putting data.

        Examples
        --------
        Suppose you have a script that calculates electronic energies of water
        molecules. For some reason, you want to continuously update the atomic
        numbers and the electronic energy under the keys
        ``/1h2o/atomic_numbers`` and ``/1h2o/energy_ele``, respectively.

        We first initialize the ``Saver`` class with the path to the file and
        the keys where we want to store the data.

        >>> saver = Saver(
        ...    'path/to/file.exdir', ('/1h2o/atomic_numbers', '/1h2o/energy_ele')
        ... )

        Then in our script we provide the data in the same order as we specified
        the keys when initializing ``saver``.

        >>> saver.save(*(np.array([8, 1, 1]), -76.35015973429272))

        .. note::

            Note we use the ``*`` operator to unpack the data; however,
            it is fine if you forget it.
        """
        rfile = File(self.rfile_path, mode="a")

        # Handles data where * is dropped.
        if isinstance(data, tuple):
            if len(data) == 1 and isinstance(data[0], tuple):
                data = data[0]

        for key, d in zip(self.data_keys, data):
            rfile.put(key, d)
        return rfile
