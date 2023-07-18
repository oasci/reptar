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


def write_qdata(file_path, R, energy=None, forces=None):
    """Write a ``qdata.txt`` file in the `ForceBalance format 
    <http://leeping.github.io/forcebalance/doc/html/usage.html>`__.

    Parameters
    ----------
    file_path : :obj:`str`
        Path to write ``qdata.txt`` file.
    R : :obj:`numpy.ndarray`, ndim: ``3``
        Cartesian coordinates for the ``coord`` field in ``qdata.txt``.
    energy : :obj:`numpy.ndarray`, ndim: ``1``, default: ``None``
        Energies used for force field fitting.
    forces : :obj:`numpy.ndarray`, ndim: ``3``, default: ``None``
        Forces used for force field fitting.
    """
    assert R.ndim == 3
    with open(file_path, "w", encoding="utf-8") as f:
        for i in range(R.shape[0]):
            f.write(f"JOB {i}\n")
