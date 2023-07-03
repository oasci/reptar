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

from abc import ABC, abstractmethod


class Parser(ABC):
    r"""Base class for parsing output files.

    Attributes
    ----------
    out_path : :obj:`str`
        Path to output file.
    file_name : :obj:`str`
        The name of the file without extension.

    """

    def __init__(self, out_path, extractors):
        """
        Parameters
        ----------
        out_path : :obj:`str`
            Path to output file.
        extractors : :obj:`list`, ndim: ``1``
            Additional extractors for the parser to use.
        """
        self.out_path = out_path
        self.file_name = ".".join(self.out_path.split("/")[-1].split(".")[:-1])
        self.parsed_info = {"system_info": {}, "runtime_info": {}, "outputs": {}}
        self.extractors = extractors

    @abstractmethod
    def parse(self):
        r"""Drive the parsing and postprocessing of parsed data.

        In general, there are usually three stages.

        :**Preprocessing**:
            Anything that is done before reptar does its parsing.
            Reptar does not have extractors for all information, so sometimes
            `cclib <https://cclib.github.io/>`__ is used to extract additional
            information.

        :**Extracting**:
            Iterates through all lines in :term:`out_path` using all extractors with
            :meth:`~reptar.parsers.Parser.extract_data_out`.
            Each extractor has its own ``parsed_info`` that needs to be added to
            :attr:`~reptar.parsers.Parser.parsed_info`.

        :**Postprocessing**:
            Any additional tasks to complete parsing.
            There is usually a ``after_parse`` method that adds or changes data
            based on the final :attr:`~reptar.parsers.Parser.parsed_info` attribute.
            Handling of other files such as :term:`geom_path` and :term:`traj_path`.
        """
        return NotImplemented

    @property
    def parsed_info(self):
        r"""Information parsed from files. Contains the following keys.

        ``system_info``
            Information specifying the system prior to any computation. Such
            as the initial cartesian coordinates, total system charge and
            multiplicity, etc.

        ``runtime_info``
            Contains information about setting up the job/calculation or running
            the job. Defining convergence criteria, parameters, etc.

        ``outputs``
            Results, requested or not, of the job. For example, SCF
            cycle values, optimized coordinates, trajectory, number of
            electrons, generated structures, etc.

        :type: :obj:`dict`
        """
        return self._parsed_info

    @parsed_info.setter
    def parsed_info(self, value):
        self._parsed_info = value

    @parsed_info.deleter
    def parsed_info(self):
        del self._parsed_info

    def after_parse(self):
        r"""Replace if desired"""

    def extract_data_out(self):
        r"""Extract data from ``out_path`` using all extractors"""
        with open(self.out_path, mode="r", encoding="utf-8") as f:
            for line in f:
                for extractor in self.extractors:
                    for trigger in extractor.triggers:
                        if trigger[0](line):
                            getattr(extractor, trigger[1])(f, line)
                            break
                    else:
                        continue
                    break
        self.combine_extracted()

    def combine_extracted(self):
        r"""Combines all parsed_info from extractors into the parsed_info in the
        parser object.
        """
        for extractor in self.extractors:
            extracted_info = extractor.parsed_info
            for key_cat in extracted_info.keys():
                for key_def in extracted_info[key_cat].keys():
                    if key_def not in self.parsed_info[key_cat].keys():
                        self.parsed_info[key_cat][key_def] = extracted_info[key_cat][
                            key_def
                        ]
                    else:
                        pass

    def map_cclib_data(self):
        r"""Assign cclib-parsed data to our ``parsed_info`` dictionary."""
        # pylint: disable=too-many-branches
        parsed_info = self.parsed_info
        cclib_data = self.cclib_data  # pylint: disable=no-member

        # Loop through every attribute in cclib_data and adds data to parsed data.
        for attr in dir(cclib_data):
            if not attr.startswith("__"):
                if attr == "atomcharges":
                    parsed_info["outputs"]["charges_lowdin"] = cclib_data.atomcharges[
                        "lowdin"
                    ]
                    parsed_info["outputs"]["charges_mulliken"] = cclib_data.atomcharges[
                        "mulliken"
                    ]

                if attr == "atomcoords":
                    parsed_info["system_info"]["geometry"] = cclib_data.atomcoords

                if attr == "atomnos":
                    parsed_info["system_info"]["atomic_numbers"] = cclib_data.atomnos

                if attr == "charge":
                    parsed_info["system_info"]["charge"] = int(cclib_data.charge)

                if attr == "mult":
                    parsed_info["system_info"]["mult"] = int(cclib_data.mult)

                if attr == "grads":
                    parsed_info["outputs"]["grads"] = cclib_data.grads

                if attr == "moenergies":
                    parsed_info["outputs"]["energy_mos"] = cclib_data.moenergies

                if attr == "moments":
                    parsed_info["outputs"]["dipole_moment"] = cclib_data.moments[1]

                if attr == "scfenergies":
                    parsed_info["outputs"]["energy_scf"] = (
                        cclib_data.scfenergies / 27.21138505
                    )  # eV -> Eh

                if attr == "nbasis":
                    parsed_info["runtime_info"]["basis_n_func"] = int(cclib_data.nbasis)

                if attr == "nelectrons":
                    parsed_info["system_info"]["n_ele"] = int(cclib_data.nelectrons)

        self.parsed_info = parsed_info
