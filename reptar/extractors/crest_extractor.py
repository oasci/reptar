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

# pylint: disable=line-too-long

import numpy as np
from .extractor import Extractor


class ExtractorCREST(Extractor):
    r"""CREST extractor"""
    # We always pass the file object into methods here.
    # pylint: disable=unused-argument

    def __init__(self, xyz_type):
        r"""
        Parameters
        ----------
        xyz_type : :obj:`str`
            Specifies the type of CREST xyz file: ``'conformer'`` or
            ``'rotamer'``.
        """
        super().__init__()
        self.xyz_type = xyz_type

    @property
    def triggers(self):
        trig = (
            (lambda line: bool("   Version " in line), "crest_version"),
            (
                lambda line: bool(
                    "total number unique points considered further :" in line
                ),
                "ensemble_info",
            ),
        )
        return trig

    def crest_version(self, f, line):
        r"""Version of CREST.

        Parameters
        ----------
        f : :obj:`io.TextIOWrapper`
            Buffered text stream of the file.
        line : :obj:`str`
            Parsed line from ``f``.

        Notes
        -----
        Example trigger text for this extractor.

        .. code-block:: text

            Version 2.12,   Thu 19. Mai 16:32:32 CEST 2022
        """
        line_split = line.split()
        version = line_split[1][:-1]
        self.parsed_info["runtime_info"]["prov_version"] = version
        next(f)

    def ensemble_info(self, f, line):
        r"""Ensemble information for conformers and rotamers.

        Parameters
        ----------
        f : :obj:`io.TextIOWrapper`
            Buffered text stream of the file.
        line : :obj:`str`
            Parsed line from ``f``.

        Notes
        -----
        Example trigger text for this extractor.

        .. code-block:: text

             total number unique points considered further :          16
                   Erel/kcal        Etot weight/tot  conformer     set   degen     origin
                   1   0.000  -254.47460    0.26798    0.26798       1       1     mtd4
                   2   0.007  -254.47459    0.26473    0.52938       2       2     mtd3
                   3   0.007  -254.47459    0.26465                                mtd6
                   4   0.893  -254.47318    0.05942    0.05942       3       1     mtd1

                   ...
        """
        if self.xyz_type == "conformer":
            # The ensemble weight is not printed in xyz, so we parse it here.
            line = self.skip_lines(f, 2)
            ensemble_weights = []
            while "T /K  " not in line:
                line_split = line.split()
                if len(line_split) == 8:
                    ensemble_weights.append(float(line_split[-4]))
                line = next(f)
            self.parsed_info["outputs"]["crest_ensemble_weights"] = np.array(
                ensemble_weights, dtype=np.float64
            )
        elif self.xyz_type == "rotamer":
            # The ensemble weight is printed in xyz with more precision.
            next(f)
