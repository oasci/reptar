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
import exdir
from .creator import creator
from pkg_resources import resource_stream
import yaml

# TODO: Find a way to handle and generate MD5 hash upon request and when saving.

class manager:
    """Manages reptar functionality with exdir.
    """

    def __init__(self):
        pass
    
    ###   Properties   ###
    @property
    def File(self):
        """The exdir file to manage. To load the exdir file, set this property
        to the path.
        """
        if hasattr(self, '_File'):
            return self._File
        else:
            return None

    def get_file(
        self, directory, mode=None, allow_remove=False, name_validation=None,
        plugins=None
    ):
        """Creates or loads an exdir file.

        Parameters
        ----------
        directory : :obj:`str`

        mode : :obj:`str`

        allow_remove : :obj:`bool`

        name_validation : :obj:`str`

        plugins : :obj:`str`

        """
        self._File = exdir.File(
            directory, mode, allow_remove, name_validation, plugins
        )
    
    def close(self):
        self.File.close()
        delattr(self, '_File')
    
    def create_group(
        self, name, parent=None, out_path=None, geom_path=None, traj_path=None,
        extractors=None
    ):
        """Creates a group using any data Reptar can find..

        Parameters
        ----------
        out_path : :obj:`str`
            Path to the quantum chemistry outfile.
        parent : ``exdir.Object``, optional
            The parent exdir object (``File`` or ``Group``). Defaults to
            loaded ``File``.
        """
        assert self.File is not None
        if parent is None:
            parent = self.File
        
        if out_path is None:
            grp = parent.create_group(name)
        else:
            reptar_creator = creator()
            reptar_creator.parse_output(
                out_path, geom_path=geom_path, traj_path=traj_path
            )
            grp = reptar_creator.create_group(name, parent)
        return grp
    
    def add_ids(self, group, entity_ids, comp_ids):
        """

        Parameters
        ----------
        entity_ids : :obj:`numpy.ndarray`

        comp_ids : :obj:`numpy.ndarray` or :obj:`list`
        """
        try:
            assert hasattr(self, '_File')
        except AssertionError:
            e = 'Must load exdir file first.'
            raise(AssertionError(e))
        
        if isinstance(comp_ids, list):
            comp_ids = np.array(comp_ids)
        
        group.create_dataset('entity_ids', data=entity_ids)
        group.create_dataset('comp_ids', data=comp_ids)
        comp_ids_num = {}
        unique_comp_ids, comp_ids_freq = np.unique(comp_ids, return_counts=True)
        for comp_id, num in zip(unique_comp_ids, comp_ids_freq):
            comp_ids_num[str(comp_id)] = int(num)
        group.attrs['comp_ids_num'] = comp_ids_num
        return group
    
    def add_definitions(self, definitions=None):
        """Add YAML file containing definitions of data included in the exdir
        file.

        The base definitions are included by default. There should be no overlap
        of any definitions. If there is, precedence is given in order of
        ``definitions`` with base always being first.

        Parameters
        ----------
        definitions : :obj:`list` (:obj:`str`), optional
            Paths to data definition YAML files to congregate. Only the name
            (not a path) is required for ones provided by reptar.
        
        Notes
        -----
        Reserved definition YAML files: ``base``, ``md``, ``molprop``, ``pes``,
        ``qc``, ``sampl``, ``solv``, ``xtb``.
        """
        try:
            assert hasattr(self, '_File')
        except AssertionError:
            e = 'Must load exdir file first.'
            raise(AssertionError(e))
        
        try:
            defs = dict(self.File.attrs['r_prov_ids'])
        except Exception as e:
            defs = {}
        
        defs_reserved = [
            'base', 'md', 'molprop', 'pes', 'qc', 'sampl', 'solv', 'xtb'
        ]
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
        self.File.attrs['definitions'] = defs
