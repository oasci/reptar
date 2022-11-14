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

from .. import File

class Saver:
    """Saves data to File.

    Useful for long driver calculations that should save periodically. This
    :obj:`reptar.calculators.save.Saver` uses provided ``data_keys`` that
    specify where to store data in ``rfile`` in the same order as the driver.
    """

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
        """Put data to rfile.

        Returns
        -------
        :obj:`reptar.File`
            File after putting data.
        """
        rfile = File(self.rfile_path, mode='a', allow_remove=False)
        for key,data in zip(self.data_keys, data):
            rfile.put(key, data)
        return rfile
