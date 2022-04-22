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
        else:
            raise TypeError(f'{f_ext} is not supported.')
        self.ftype = f_ext[1:]
        self.File = File
    
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
        if '//' in key:
            key = key.replace('//', '/')
        if key[0] != '/':
            key = '/' + key
        key_split = key.rsplit('/', 1)
        if key_split[0] == '':
            key_split[0] = '/'
        return key_split

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
        key_parent, key_data = self.split_key(key)
        
        if self.ftype == 'exdir':
            if key == '/':
                return self.File
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
    
    def _add_to_exdir(self, parent_key, data_key, data):
        """Add data to an exdir group.

        Parameters
        ----------
        parent_key : :obj:`str`
        data_key : :obj:`str`
        data : :obj:`dict`
        
        """
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
        if isinstance(data, np.ndarray) or isinstance(data, list) or isinstance(data, tuple):
            # If data is a list we check the types of data it contains.
            # If all of them are strings, we do not convert to array.
            if not isinstance(data, np.ndarray):
                # If there is only one item in the list we can likely
                # just store it as an attribute. However, we have
                # to be careful when the single item is an array.
                if len(data) == 1:
                    data = data[0]
                    if isinstance(data, np.ndarray):
                        store_as = 'dset'
                    else:
                        store_as = 'attr'
                # If all the items are not all strings then we make
                # the dataset.
                elif not all(isinstance(i, str) for i in data):
                    data = np.array(data)
                    store_as = 'dset'
                # At this point only data that contains all strings
                # should be left. We put these as attributes (i.e., 
                # comp_ids)
                else:
                    store_as = 'attr'
            # Handling arrays.
            else:
                # If only one item we store it as an attribute.
                # Otherwise we create the dataset.
                if data.shape == (1,):
                    data = data[0].item()
                    store_as = 'attr'
                else:
                    store_as = 'dset'
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
        parent_key, data_key = self.split_key(key)
        if self.ftype == 'exdir':
            self._add_to_exdir(parent_key, data_key, data)
    
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
    