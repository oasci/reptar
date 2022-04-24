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
import numpy as np
import json
import exdir
import cclib
from collections import defaultdict
from .utils import combine_dicts

class data:
    """A reptar data object that creates, stores, retrieves, and manages data.

    Parameters
    ----------
    file_path : :obj`str`
        Path to a file supported by reptar. If it does not exist, then one will
        be created if possible.
    mode : :obj:`str`, optional
        A file mode string that defines the read/write behavior. Defaults to
        ``'r'``.
    allow_remove : :obj:`bool`, optional
        Allow the removal of exdir groups in ``w`` operation. Defaults to
        ``False``.
    name_validation : :obj:`str`, optional
        Validation mode for exdir group names. Defaults to ``'thorough'``.
    plugins : :obj:`list`, optional
        A list of instantiated exdir plugins. Defaults to ``None``.
    """

    def __init__(self, file_path, mode='r', allow_remove=False,
        name_validation='thorough', plugins=None
    ):
        exists = os.path.exists(file_path)
        _, f_ext = os.path.splitext(file_path)

        if f_ext == '.exdir':
            File = exdir.File(
                file_path, mode, allow_remove, name_validation, plugins
            )
        elif f_ext == '.json':
            if exists:
                with open(file_path, 'r') as f:
                    File = json.load(f)
            else:
                File = defaultdict(lambda: defaultdict(dict))
        else:
            raise TypeError(f'{f_ext} is not supported.')
        
        self.fpath = file_path
        self.ftype = f_ext[1:]
        self.File = File
    
    def clean_key(self, key):
        """Clean key and remove any common mistakes in keys.
        
        Parameters
        ----------
        key : :obj:`str`
            Key of the desired data.
        
        Returns
        -------
        :obj:`str`
            Cleaned key.
        """
        if '//' in key:
            key = key.replace('//', '/')
        return key
    
    def split_key(self, key):
        """Split the key into a parent and data key.

        Parameters
        ----------
        key : :obj:`str`
            Key of the desired data. Nested keys should be separated by ``/``.
        
        Returns
        -------
        :obj:`str`
            Parent key.
        :obj:`str`
            Data key.
        """
        if key[0] != '/':
            key = '/' + key
        key_split = key.rsplit('/', 1)
        if key_split[0] == '':
            key_split[0] = '/'
        return key_split
    
    def _get_from_dict(self, key):
        """Get data from dictionary-like file.
        
        Parameters
        ----------
        key : :obj:`str`
            Key of the desired data. Nested keys should be separated by ``/``.
        """
        if key == '/':
            return self.File
        
        keys = key.split('/')
        if keys[0] == '':
            del keys[0]
        
        data = self.File[keys[0]]
        for k in keys[1:]:
            data = data[k]
        if isinstance(data, list):
            data = np.array(data)
        return data
    
    def _get_from_exdir(self, key):
        """Get data from exdir file source.
        
        Parameters
        ----------
        key : :obj:`str`
            Key of the desired data. Nested keys should be separated by ``/``.
        """
        if key == '/':
            return self.File
        
        key_parent, key_data = self.split_key(key)
        parent = self.File[key_parent]
        if key_data in list(parent):
            data = parent[key_data]
            if isinstance(data, exdir.core.dataset.Dataset):
                data = data.data
        elif key_data in list(parent.attrs):
            data = parent.attrs[key_data]
            if isinstance(data, exdir.core.attribute.Attribute):
                data = dict(data)
        
        return data
    
    def get_keys(self, group_key):
        """A list of keys in a group (not included nested keys).
        
        Parameters
        ----------
        group_key : :obj:`str`
            Key specifying which group to get the data keys from.
        
        Returns
        -------
        :obj:`list`
            A list of all available keys not leading to another group.
        """
        group_key = self.clean_key(group_key)
        keys = []
        group = self.get(group_key)
        
        if self.ftype == 'exdir':
            keys.extend(list(sorted(group.attrs.keys())))
            keys.extend(list(sorted(group.keys())))
        elif self.ftype == 'json':
            keys.extend(list(group.keys()))
        return keys

    def get(self, key):
        """Retrives data from a loaded file.

        Provides a uniform method of query data from different sources.

        Parameters
        ----------
        key : :obj:`str`
            Key of the desired data. Nested keys should be separated by ``/``.
        
        Examples
        --------
        >>> data.get('/prod_1/charge')
        0
        >>> data.get('energy_scf')
        -12419.360138637763
        """
        key = self.clean_key(key)
        if self.ftype == 'exdir':
            data = self._get_from_exdir(key)
        elif self.ftype == 'json':
            data = self._get_from_dict(key)
        return data
    
    def _is_iter(self, data):
        if isinstance(data, np.ndarray) or isinstance(data, list) or isinstance(data, tuple):
            return True
        else:
            return False
    
    def simplify_iter_data(self, data):
        """Simplify iterative data if possible.
        
        Parameters
        ----------
        data : :obj:`numpy.ndarray`, :obj:`list`, :obj:`tuple`
            An iterative data object.
        """
        # If data is a list we check the types of data it contains.
        # If all of them are strings, we do not convert to array.
        if not isinstance(data, np.ndarray):
            # If there is only one item in the list we can likely
            # just store it as an attribute. However, we have
            # to be careful when the single item is an array.
            if len(data) == 1:
                data = data[0]
            # If all the items are not all strings then we make
            # the dataset.
            elif not all(isinstance(i, str) for i in data):
                data = np.array(data)
            # At this point only data that contains all strings
            # should be left. We put these as attributes (i.e., 
            # comp_ids)
            else:
                pass
        # Handling arrays.
        else:
            # If only one item we store it as a value.
            if data.shape == (1,):
                data = data[0].item()
        
        return data
    
    def _add_to_dict(self, key, data):
        """Add data to dictionary-like file.

        Parameters
        ----------
        key : :obj:`str`
            Where to add the data. Can be a nested key.
        data : various
            Data to add to file.
        """
        if isinstance(data, np.ndarray):
            if 'U' in str(data.dtype):
                data = data.tolist()

        keys = key.split('/')
        if keys[0] == '':
            del keys[0]
        
        is_iter = self._is_iter(data)
        if is_iter:
            data = self.simplify_iter_data(data)

        add_dic = {keys[-1]: data}
        for key in reversed(keys[:-1]):
            add_dic = {key: add_dic}

        self.File = combine_dicts(self.File, add_dic)
    
    def _add_to_exdir(self, key, data):
        """Add data to an exdir group.

        Parameters
        ----------
        key : :obj:`str`
            Where to add the data. Can be a nested key.
        data : various
            Data to add to exdir file.
        """
        parent_key, data_key = self.split_key(key)
        group = self.get(parent_key)
        assert isinstance(group, exdir.core.exdir_file.File) \
            or isinstance(group, exdir.core.group.Group)
        
        # Custom handling of parsed data to ensure logical behavior.
        if data_key == 'dipole_moment':
            if isinstance(data, np.ndarray):
                data = data.tolist()
            group.attrs[data_key] = data
            return None

        # If data is not array/list or array/list of one dimension and
        # length we make it an exdir attribute.
        # We make anything else an exdir dataset.
        is_iter = self._is_iter(data)
        if is_iter:
            data = self.simplify_iter_data(data)
            if isinstance(data, np.ndarray):
                store_as = 'dset'
            else:
                store_as = 'attr'
        # Handles all non iterable data. 
        else:
            store_as = 'attr'
        
        if store_as == 'attr':
            group.attrs[data_key] = data
        elif store_as == 'dset':
            try:
                group.create_dataset(data_key, data=data)
            except RuntimeError:
                # If data set already exists, exdir will throw a RuntimeError.
                # We remove this dataset and try again.
                group.__delitem__(data_key)
                group.create_dataset(data_key, data=data)

    def add(self, key, data):
        """Add information to a loaded file.

        Parameters
        ----------
        key : :obj:`str`
            Key of the desired data (including parent). Nested keys should be
            separated by ``/``.
        data
            Data to add to file.
        """
        key = self.clean_key(key)
        if self.ftype == 'exdir':
            self._add_to_exdir(key, data)
        elif self.ftype == 'json':
            self._add_to_dict(key, data)
    
    def init_group(self, key):
        """Create a new group with the specified key.

        Parameters
        ----------
        key : :obj:`str`
            Key of the desired data (including parent). Nested keys should be
            separated by ``/``.
        """
        parent_key, group_key = self.split_key(key)
        if self.ftype == 'exdir':
            parent = self.get(parent_key)
            group = parent.create_group(group_key)
        return group
    
    def save(self, json_prettify=True):
        """Saves non-exdir files.

        Parameters
        ----------
        json_prettify : :obj:`bool`
            Indents JSON objects if True. If false the JSON file is only one
            line.
        """
        if self.ftype == 'json':
            json_dict = self.File

            if json_prettify:
                json_string = json.dumps(
                    json_dict, cls=cclib.io.cjsonwriter.JSONIndentEncoder,
                    sort_keys=True, indent=4
                )
            else:
                json_string = json.dumps(
                    json_dict, cls=cclib.io.cjsonwriter.NumpyAwareJSONEncoder,
                    sort_keys=True
                )
            with open(self.fpath, 'w') as f:
                f.write(json_string)
    