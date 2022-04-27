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

import numpy as np
from .creator import creator
from .data import data
from pkg_resources import resource_stream
import yaml

# TODO: Find a way to handle and generate MD5 hash upon request and when saving.

defs_reserved = [
    'base', 'md', 'molprop', 'pes', 'qc', 'sampl', 'solv', 'xtb'
]

class manager:
    """Manages reptar files.
    """

    def __init__(self):
        pass
    
    def load(
        self, file_path, mode='r', allow_remove=False, plugins=None
    ):
        """Creates a data object.

        Parameters
        ----------
        file_path : :obj:`str`
            Path to a file supported by reptar. If it does not exist, then one will
            be created if possible.
        mode : :obj:`str`, optional
            A file mode string that defines the read/write behavior. Defaults to
            ``'r'``.
        allow_remove : :obj:`bool`, optional
            Allow the removal of exdir groups in ``w`` operation. Defaults to
            ``False``.
        plugins : :obj:`list`, optional
            A list of instantiated exdir plugins. Defaults to ``None``.
        """
        self.data = data(
            file_path, mode=mode, allow_remove=allow_remove, plugins=plugins
        )
    
    def create_group(
        self, group_key, out_path=None, geom_path=None, traj_path=None,
        extractors=None
    ):
        """Creates a group using any data reptar can find from provided
        paths.

        Parameters
        ----------
        group_key : :obj:`str`, optional
            The key to the desired new group (including parent if applicable).
        out_path : :obj:`str`, optional
            Path to the main log file generated by the package.
        geom_path : :obj:`str`, optional
            Path to a file containing a single geometry.
        traj_path : :obj:`str`, optional
            Path to a trajectory file from a geometry optimization, MD
            simulation, etc.
        extractors : :obj:`list`, optional
            Additional extractors for the parser to use.
        Notes
        -----
        If both ``geom_path`` and ``traj_path`` are provided, it is assumed that
        ``geom_path`` provides an initial geometry not included in
        ``traj_path``.
        """
        reptar_creator = creator(data=self.data)
        reptar_creator.parse_output(
            out_path, geom_path=geom_path, traj_path=traj_path,
            extractors=extractors
        )
        grp = reptar_creator.create_group(group_key)
        self.data = reptar_creator.data
        return grp
    
    def add_ids(self, group_key, entity_ids, comp_ids):
        """Add ``entity_ids`` and ``comp_ids`` to a group.

        Parameters
        ----------
        group_key : :obj:`str`
            Key to the group to add IDs to (including parent).
        entity_ids : :obj:`numpy.ndarray`
            A uniquely identifying integer specifying what atoms belong to which
            entities. Entities can be a related set of atoms, molecules, or
            functional group.
        comp_ids : :obj:`numpy.ndarray`
            Relates ``entity_id`` to a fragment label for chemical components
            or species. Labels could be ``WAT`` or ``h2o`` for water, ``MeOH``
            for methanol, ``bz`` for benzene, etc. There are no standardized
            labels for species. The index of the label is the respective
            ``entity_id``.
        """
        if isinstance(comp_ids, list):
            comp_ids = np.array(comp_ids)
        
        self.data.add(f'{group_key}/entity_ids', entity_ids)
        self.data.add(f'{group_key}/comp_ids', comp_ids)

        comp_ids_num = {}
        unique_comp_ids, comp_ids_freq = np.unique(comp_ids, return_counts=True)
        for comp_id, num in zip(unique_comp_ids, comp_ids_freq):
            comp_ids_num[str(comp_id)] = int(num)
        
        self.data.add(f'{group_key}/comp_ids_num', comp_ids_num)
        return self.data.get(group_key)
    
    # TODO: Update ways we get data here.
    def add_definitions(self, definitions=None):
        """Add definitions of data to the file.

        The base definitions are included by default. There should be no overlap
        of any definitions. If there is, precedence is given in order of
        ``definitions`` with base always being first.

        Parameters
        ----------
        definitions : :obj:`list` [:obj:`str`], optional
            Paths to data definition YAML files to congregate. Only the name
            (not a path) is required for ones provided by reptar.
        
        Notes
        -----
        Reserved definition YAML files: ``base``, ``md``, ``molprop``, ``pes``,
        ``qc``, ``sampl``, ``solv``, ``xtb``.
        """        
        defs = {}

        if definitions is None:
            definitions = ['base']
        else:
            definitions.insert(0, 'base')
        
        for def_path in definitions:
            if def_path in defs_reserved:
                stream = resource_stream(
                    'reptar.definitions', f'{def_path}.yaml'
                )
            else:
                stream = open(def_path)
            def_add = yaml.safe_load(stream)
            
            # Add all definitions.
            for key_cat in def_add.keys():
                if key_cat not in defs.keys():
                    defs[key_cat] = {}
                for key_def in def_add[key_cat].keys():
                    if key_def not in defs[key_cat].keys():
                        defs[key_cat][key_def] = def_add[key_cat][key_def]
                    else:
                        pass
        
        self.data.add('definitions', defs)
