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

import exdir
import json
import numpy as np
import os
from .utils import combine_dicts, get_md5

class File:
    """Create, store, and access data from a variety of formats.
    """

    def __init__(self, file_path, mode='r', allow_remove=False,
        plugins=None, from_dict=None
    ):
        """
        Parameters
        ----------
        file_path : :obj:`str`
            Path to a file supported by reptar. If it does not exist, then one will
            be created if possible.
        mode : :obj:`str`, default: ``'r'``
            A file mode string that defines the read/write behavior.
        allow_remove : :obj:`bool`, default: ``False``
            Allow the removal of exdir groups in ``'w'`` operation.
        plugins : :obj:`list`, default: ``None``
            A list of instantiated exdir plugins.
        from_dict : :obj:`dict`, default: ``None``
            Load data from a dictionary.
        """
        if from_dict is None:
            self._from_path(
                file_path, mode, allow_remove, plugins
            )
        else:
            self._from_dict(
                file_path, from_dict, mode, allow_remove, plugins
            )
    
    def _from_path(
        self, file_path, mode, allow_remove, plugins
    ):
        """Populates the File object from a file path.

        Parameters
        ----------
        file_path : :obj:`str`
            Path to a file supported by reptar. If it does not exist, then one will
            be created if possible.
        mode : :obj:`str`, default: ``'r'``
            A file mode string that defines the read/write behavior.
        allow_remove : :obj:`bool`, default: ``False``
            Allow the removal of exdir groups in ``'w'`` operation.
        plugins : :obj:`list`, default: ``None``
            A list of instantiated exdir plugins.
        """
        exists = os.path.exists(file_path)
        _, f_ext = os.path.splitext(file_path)

        if f_ext == '.exdir':
            File_ = exdir.File(
                file_path, mode, allow_remove, plugins
            )
        elif f_ext == '.json':
            if exists:
                with open(file_path, 'r') as f:
                    File_ = json.load(f)
            else:
                File_ = {}
        elif f_ext == '.npz':
            if exists:
                File_ = dict(np.load(file_path, allow_pickle=True))
                # Since everything is stored in arrays, we clean up array data.
                for k,v in File_.items():
                    if self._is_iter(v):
                        v = self.simplify_iter_data(v)
                        File_[k] = v
            else:
                File_ = {}
        else:
            raise TypeError(f'{f_ext} is not supported.')
        
        self.fpath = file_path
        self.ftype = f_ext[1:]
        self.fmode = mode
        self.File_ = File_
    
    def _from_dict(
        self, file_path, group_dict, mode, allow_remove, plugins
    ):
        """Populates the File object from a dictionary.

        Parameters
        ----------
        file_path : :obj:`str`
            Path to a file supported by reptar. If it does not exist, then one
            will be created if possible.
        group_dict : :obj:`dict`
            Dictionary to populate the File object with.
        mode : :obj:`str`, default: ``'r'``
            A file mode string that defines the read/write behavior.
        allow_remove : :obj:`bool`, default: ``False``
            Allow the removal of exdir groups in ``'w'`` operation.
        plugins : :obj:`list`, default: ``None``
            A list of instantiated exdir plugins.
        """
        _, f_ext = os.path.splitext(file_path)
        
        if f_ext == '.exdir':
            self.File_ = exdir.File(
                file_path, mode, allow_remove, plugins
            )
        elif f_ext == '.json' or f_ext == '.npz':
            self.File_ = {}
        self.fpath = file_path
        self.ftype = f_ext[1:]
        self.fmode = mode

        for pair in self._iter_dict(group_dict):
            key = '/'.join(pair[:-1])
            data = pair[-1]
            self.add(key, data)
    
    def clean_key(self, key):
        """Clean key and remove any common mistakes.

        Currently this only corrects instances of ``//``.
        
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
        """Get data from dictionary-like file (e.g., json, npz).
        
        Parameters
        ----------
        key : :obj:`str`
            Key of the desired data. Nested keys should be separated by ``/``.
        
        Returns
        -------
        Requested data from a dictionary source.
        """
        if key == '/':
            return self.File_
        
        keys = key.split('/')
        if keys[0] == '':
            del keys[0]
        
        data = self.File_[keys[0]]
        for k in keys[1:]:
            data = data[k]
        if isinstance(data, list):
            data_array = np.array(data)
            if data_array.dtype != 'O':
                data = data_array
        return data
    
    def _get_from_exdir(self, key):
        """Get data from exdir file.
        
        Parameters
        ----------
        key : :obj:`str`
            Key of the desired data. Nested keys should be separated by ``/``.
        
        Returns
        -------
        Requested data from a exdir source.
        """
        if key == '/':
            return self.File_
        
        key_parent, key_data = self.split_key(key)
        parent = self.File_[key_parent]
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
        """A list of keys in a group.

        Does not include keys of nested groups.
        
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
        group = self.get(group_key)

        keys = []
        if self.ftype == 'exdir':
            keys.extend(list(sorted(group.attrs.keys())))
            keys.extend(
                list(sorted(
                    key for key in group.keys() if not \
                    isinstance(group[key], exdir.core.raw.Raw)
                ))
            )
        elif self.ftype == 'json' or self.ftype == 'npz':
            keys.extend(list(group.keys()))
        return keys

    def get(self, key):
        """Retrieve data.

        Parameters
        ----------
        key : :obj:`str`
            Key of the desired data. Nested keys should be separated by ``/``.
        
        Examples
        --------
        >>> rfile.get('/prod_1/charge')
        0
        >>> rfile.get('energy_scf')
        -12419.360138637763
        """
        key = self.clean_key(key)
        if self.ftype == 'exdir':
            data = self._get_from_exdir(key)
        elif self.ftype == 'json' or self.ftype == 'npz':
            data = self._get_from_dict(key)
        return data
    
    def _is_iter(self, data):
        """If data is iterative (i.e., array, list, or tuple).
        
        Returns
        -------
        :obj:`bool`
            If the data is iterative.
        """
        if isinstance(data, np.ndarray) or isinstance(data, list) or isinstance(data, tuple):
            return True
        else:
            return False
    
    def simplify_iter_data(self, data):
        """Simplify iterative data if possible.

        Attempts to simplify iterative objects by the type and/or number of
        elements. For example, iterative data with only one element will return
        just that element.

        This is mainly used when adding data to a file. When parsing data we
        will not know if there will be multiple values of certain data (e.g.,
        multiple single point energies) until we finish parsing. Thus, we assume
        there will be multiple values while parsing then simplify any iterative
        objects at the end.
        
        Parameters
        ----------
        data : :obj:`numpy.ndarray`, :obj:`list`, :obj:`tuple`
            An iterative data object.
        
        Returns
        -------
        ``obj``
            Simplified data (if possible).
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
                data_array = np.array(data)
                # Sometimes we have a list of objects that cannot be easily
                # converted to numpy data types. Numpy will use an object
                # dtype, so we only convert to array if its not an 'O' dtype.
                if data_array.dtype != 'O':
                    data = data_array
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

        self.File_ = combine_dicts(self.File_, add_dic)
    
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
        elif data_key == 'wall_potential':
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
        """Add data to file.

        Note that there is some data postprocessing using
        :meth:`~reptar.reptar_file.File.simplify_iter_data`.

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
        elif self.ftype == 'json' or self.ftype == 'npz':
            self._add_to_dict(key, data)
    
    def _iter_dict(self, dic):
        """Iterate over nested dictionary.

        Parameters
        ----------
        dic : :obj:`dict`
            An arbitrarily nested dictionary.
        
        Yields
        ------
        A :obj:`tuple` of keys where the last element is the value.
        """
        for k, v in dic.items():
            if isinstance(v, dict):
                # If value is dict then iterate over all its values
                for pair in self._iter_dict(v):
                    yield (k, *pair)
            else:
                # If value is not dict type then yield the value
                yield (k, v)
    
    def init_group(self, key):
        """Initialize an exdir group with the specified key.

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
    
    def as_dict(self, group_key):
        """Get a group as a dictionary.

        Parameters
        ----------
        group_key : :obj:`str`
            Desired group.
        
        Returns
        -------
        :obj:`dict`
            The desired group as a dictionary.
        """
        group = self.get(group_key)
        if self.ftype == 'json' or self.ftype == 'npz':
            return group
        elif self.ftype == 'exdir':
            group_dict = {}
            data_keys = self.get_keys(group_key)
            get_keys = [f'{group_key}/{d_key}' for d_key in data_keys]
            for d_key,g_key in zip(data_keys, get_keys):
                group_dict[d_key] = self.get(g_key)
            return group_dict
    
    def update_md5(self, group_key):
        """Update all possible MD5 hashes of a specific group.

        Parameters
        ----------
        group_key : :obj:`str`
            Desired group.
        """
        md5 = get_md5(self, group_key)
        self.add(f'{group_key}/md5', md5)
        
        try:
            md5_arrays = get_md5(self, group_key, only_arrays=True)
            self.add(f'{group_key}/md5_arrays', md5_arrays)
        except Exception:
            pass

        try:
            md5_structures = get_md5(self, group_key, only_structures=True)
            self.add(f'{group_key}/md5_structures', md5_structures)
        except Exception:
            pass
    
    def save(self, json_prettify=True):
        """Saves non-exdir files.

        Parameters
        ----------
        json_prettify : :obj:`bool`
            Indents JSON objects if True. If false the JSON file is only one
            line.
        """
        assert self.fmode == 'w'
        from cclib.io.cjsonwriter import JSONIndentEncoder, NumpyAwareJSONEncoder
        if self.ftype == 'json':
            json_dict = self.File_

            if json_prettify:
                json_string = json.dumps(
                    json_dict, cls=JSONIndentEncoder,
                    sort_keys=True, indent=4
                )
            else:
                json_string = json.dumps(
                    json_dict, cls=NumpyAwareJSONEncoder,
                    sort_keys=True
                )
            with open(self.fpath, 'w') as f:
                f.write(json_string)
        elif self.ftype == 'npz':
            npz_dict = self.File_
            np.savez_compressed(self.fpath, **npz_dict)
