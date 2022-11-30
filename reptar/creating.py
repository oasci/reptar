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

import os
import shutil
import ase
from ase.io.ulm import InvalidULMFileError
import numpy as np
from pkg_resources import resource_stream
import yaml
from . import _version
from .reptar_file import File

# pylint: disable=no-name-in-module
from .parsers import ParserORCA, ParserXTB, ParserASE, ParserCREST

__version__ = _version.get_versions()["version"]

defs_reserved = ["base", "md", "molprop", "pes", "qc", "sampl", "solv", "xtb"]


def identify_parser(out_path):
    r"""Identifies the correct parser depending on some trigger.
    Only supported packaged should be included in ``triggers``.

    Parameters
    ----------
    out_path : :obj:`str`
        Path to output file.

    Returns
    -------
    :obj:`reptar.Creator`
        One of the supported creator classes.
    """
    with open(out_path, "r", encoding="utf-8") as f:
        for line in f:
            for parser, phrases, do_break in triggers:
                if all(line.lower().find(p.lower()) >= 0 for p in phrases):
                    filetype = parser
                    if do_break:
                        return filetype
    return None


# Triggers to identify output files.
triggers = [
    (ParserORCA, ["O   R   C   A"], True),
    (ParserXTB, ["x T B"], True),
    (ParserCREST, ["C R E S T"], True),
]


def identify_trajectory(traj_path):
    r"""Identifies the type of trajectory depending on a series of tests."""
    # ASE trajectory
    try:
        ase.io.trajectory.Trajectory(traj_path)
        return ParserASE
    except InvalidULMFileError:
        # Does have ase installed but is not a trajectory file.
        pass
    return None


class Creator:
    r"""Create groups from computational chemistry data."""

    def __init__(self, rfile=None):
        r"""
        Parameters
        ----------
        rfile : :obj:`reptar.File`, default: ``None``
            An initialized reptar File.
        """
        if rfile is not None:
            self.rfile = rfile

    def load(self, file_path, mode="r", allow_remove=False, plugins=None):
        r"""Load a reptar file for creating/adding information.

        Parameters
        ----------
        file_path : :obj:`str`
            Path to a file supported by reptar. If it does not exist, then one
            will be created if possible.
        mode : :obj:`str`, default: ``'r'``
            A file mode string that defines the read/write behavior.
        allow_remove : :obj:`bool`, default: ``False``
            Allow the removal of exdir groups in ``w`` operation.
        plugins : :obj:`list`, default: ``None``
            A list of instantiated exdir plugins.
        """
        self.rfile = File(
            file_path, mode=mode, allow_remove=allow_remove, plugins=plugins
        )

    @property
    def rfile(self):
        r"""The reptar file to manage.

        :obj:`reptar.File`
        """
        return self._rfile

    @rfile.setter
    def rfile(self, value):
        self._rfile = value

    @rfile.deleter
    def rfile(self):
        del self._rfile

    def parse_output(
        self, out_path, geom_path=None, traj_path=None, extractors=None, **kwargs
    ):
        r"""Parse output file using cclib and custom parser.

        Sets the ``reptar.creator.parser`` and ``reptar.creator.parsed_info``
        attributes.

        Parameters
        ----------
        out_path : :obj:`str`
            Path to computational chemistry output file.
        geom_path : :obj:`str`, default: ``None``
            Path to a file containing a single geometry.
        traj_path : :obj:`str`, default: ``None``
            Path to a trajectory file from a geometry optimization, MD
            simulation, etc.
        extractors : :obj:`list`, default: ``None``
            Additional extractors for the parser to use.
        **kwargs
            Additional paths to output files. Passed to the package's parser.
        """
        # Handle paths.
        self.out_path = os.path.abspath(out_path)
        if geom_path is not None:
            self.geom_path = os.path.abspath(geom_path)
        if traj_path is not None:
            self.traj_path = os.path.abspath(traj_path)

        # Identify and run the parser.
        package_parser = identify_parser(self.out_path)
        self.parser = package_parser(
            self.out_path,
            geom_path=geom_path,
            traj_path=traj_path,
            extractors=extractors,
            **kwargs,
        )
        self.parsed_info = self.parser.parse()

    def parse_traj(self, traj_path, extractors=None):
        r"""Parse a trajectory file.

        Parameters
        ----------
        traj_path : :obj:`str`
            Path to a trajectory file from a geometry optimization, MD
            simulation, etc.
        extractors : :obj:`list`, default: ``None``
            Additional extractors for the parser to use.
        """
        self.traj_path = os.path.abspath(traj_path)
        package_parser = identify_trajectory(self.traj_path)
        self.parser = package_parser(
            out_path=None, geom_path=None, traj_path=traj_path, extractors=extractors
        )
        self.parsed_info = self.parser.parse()

    def _create_extras_xtb(self, group_key):
        r"""Extra actions for groups created from xtb calculations.

        Parameters
        ----------
        group_key : :obj:`str`
            Key to the desired new group (including parent).
        """
        out_dir = os.path.dirname(self.out_path)
        group = self.rfile.get(group_key)

        if self.rfile.ftype == "exdir":
            # Copy MD and xTB restart files if possible.
            md_restart_path = f"{out_dir}/mdrestart"
            if os.path.exists(md_restart_path):
                raw = group.require_raw("restart_files")
                raw_path = os.path.join(raw.root_directory, raw.relative_path)
                shutil.copy(md_restart_path, os.path.join(raw_path, "mdrestart"))
            xtb_restart_path = f"{out_dir}/xtbrestart"
            if os.path.exists(xtb_restart_path):
                raw = group.require_raw("restart_files")
                raw_path = os.path.join(raw.root_directory, raw.relative_path)
                shutil.copy(xtb_restart_path, os.path.join(raw_path, "xtbrestart"))

    def from_calc(
        self,
        group_key,
        out_path=None,
        geom_path=None,
        traj_path=None,
        extractors=None,
        **kwargs,
    ):
        r"""Create a group from a supported calculation.

        Parameters
        ----------
        group_key : :obj:`str`
            Key to the desired new group (including parent).
        out_path : :obj:`str`, default: ``None``
            Path to the main log file generated by the package.
        geom_path : :obj:`str`, default: ``None``
            Path to a file containing a single geometry.
        traj_path : :obj:`str`, default: ``None``
            Path to a trajectory file from a geometry optimization, MD
            simulation, etc.
        extractors : :obj:`list`, default: ``None``
            Additional extractors for the parser to use.
        **kwargs
            Additional paths to output files. Passed into
            ``creator.parse_output()``.

        Returns
        -------
        :obj:`reptar.File`
            The entire reptar file after creating the new group.

        Notes
        -----
        If both ``geom_path`` and ``traj_path`` are provided, it is assumed that
        ``geom_path`` provides an initial geometry not included in
        ``traj_path``.
        """
        assert hasattr(self, "rfile")

        if self.rfile.ftype == "exdir":
            self.rfile.create_group(group_key)

        # Parsable calculations using an output file.
        if out_path is not None:
            self.parse_output(
                out_path,
                geom_path=geom_path,
                traj_path=traj_path,
                extractors=extractors,
                **kwargs,
            )
            parsed_info = self.parsed_info

        # Only a trajectory is provided. Likely coordinates or package-specific
        # trajectory file like ASE.
        elif out_path is None and traj_path is not None:
            self.parse_traj(traj_path, extractors=extractors)
            parsed_info = self.parsed_info
        else:
            return None

        # Add all data to group.
        self.rfile.put_all(group_key, parsed_info, nested=True)

        # Extra stuff to do depending on package.
        if self.parser.package == "xtb":
            self._create_extras_xtb(group_key)

        # MD5 stuff
        self.rfile.update_md5(group_key)

        # Adding version
        self.rfile.put(f"{group_key}/reptar_version", __version__)

        return self.rfile

    def ids(self, group_key, entity_ids, comp_ids):
        r"""Add ``entity_ids`` and ``comp_ids`` to a group.

        Parameters
        ----------
        group_key : :obj:`str`
            Key to the group to add IDs to (including parent).
        entity_ids : :obj:`numpy.ndarray`, ndim: ``1``
            A uniquely identifying integer specifying what atoms belong to which
            entities. Entities can be a related set of atoms, molecules, or
            functional group.
        comp_ids : :obj:`numpy.ndarray`, ndim: ``1``
            Relates ``entity_id`` to a fragment label for chemical components
            or species. Labels could be ``WAT`` or ``h2o`` for water, ``MeOH``
            for methanol, ``bz`` for benzene, etc. There are no standardized
            labels for species. The index of the label is the respective
            ``entity_id``.

        Returns
        -------
        ``obj``
            The reptar file group.
        """
        if isinstance(comp_ids, list):
            comp_ids = np.array(comp_ids)

        self.rfile.put(f"{group_key}/entity_ids", entity_ids)
        self.rfile.put(f"{group_key}/comp_ids", comp_ids)

        comp_ids_num = {}
        unique_comp_ids, comp_ids_freq = np.unique(comp_ids, return_counts=True)
        for comp_id, num in zip(unique_comp_ids, comp_ids_freq):
            comp_ids_num[str(comp_id)] = int(num)

        self.rfile.put(f"{group_key}/comp_ids_num", comp_ids_num)
        return self.rfile.get(group_key)

    # TODO: Update ways we get data here.
    def definitions(self, definitions=None):
        r"""Add definitions of data to the file.

        The base definitions are included by default. There should be no overlap
        of any definitions. If there is, precedence is given in order of
        ``definitions`` with base always being first.

        Parameters
        ----------
        definitions : :obj:`list` [:obj:`str`], default: ``None``
            Paths to data definition YAML files to congregate. Only the name
            (not a path) is required for ones provided by reptar.

        Notes
        -----
        Reserved definition YAML files: ``base``, ``md``, ``molprop``, ``pes``,
        ``qc``, ``sampl``, ``solv``, ``xtb``.
        """
        defs = {}

        if definitions is None:
            definitions = ["base"]
        else:
            definitions.insert(0, "base")

        for def_path in definitions:
            if def_path in defs_reserved:
                stream = resource_stream("reptar.definitions", f"{def_path}.yaml")
            else:
                # pylint: disable-next=consider-using-with
                stream = open(def_path, encoding="utf-8")
            def_add = yaml.safe_load(stream)

            # Add all definitions.
            for key_cat in def_add.keys():
                if key_cat not in defs:
                    defs[key_cat] = {}
                for key_def in def_add[key_cat].keys():
                    if key_def not in defs[key_cat].keys():
                        defs[key_cat][key_def] = def_add[key_cat][key_def]
                    else:
                        pass

        self.rfile.put("definitions", defs)
