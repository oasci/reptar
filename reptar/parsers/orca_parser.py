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

from .parser import Parser  # pylint: disable=no-name-in-module
from ..extractors import ExtractorORCA  # pylint: disable=no-name-in-module


class ParserORCA(Parser):
    """Custom parser for ORCA output files."""

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
        self.package = "orca"
        if extractors is None:
            extractors = []
        extractors.insert(0, ExtractorORCA())
        super().__init__(out_path, extractors)
        # TODO: Handle geom and traj paths

    def parse(self):
        """Parses output file and extracts information."""
        # cclib parsed information.
        try:
            import cclib  # pylint: disable=import-outside-toplevel

            self.cclib_data = cclib.io.ccread(self.out_path)
            self.map_cclib_data()
        except Exception as e:
            raise e

        # Extract information.
        self.extract_data_out()

        # Postprocessing
        self.after_parse()
        # TODO: Figure out how to get driver.

        return self.parsed_info

    def after_parse(self):
        """Checks to perform after parsing output file."""
        if "scf_grid_level_final" not in self.parsed_info["runtime_info"].keys():
            self.parsed_info["runtime_info"]["scf_grid_level_final"] = self.parsed_info[
                "runtime_info"
            ]["scf_grid_level"]

        if "conv_val_geo_energy" in self.parsed_info["runtime_info"].keys():
            self.parsed_info["runtime_info"]["conv_val_geo_energy"].insert(0, 0.0)
