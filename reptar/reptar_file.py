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
                        v = self.simplify_iter_data(v, k.split('/')[-1])
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
            self.put(key, data)
    
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
    
    def _get_from_exdir(self, key, as_memmap=False):
        """Get data from exdir file.
        
        Parameters
        ----------
        key : :obj:`str`
            Key of the desired data. Nested keys should be separated by ``/``.
        as_memmap : :obj:`bool`, default: ``False``
            Keep NumPy memmap instead of converting to arrays.
        
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
                if as_memmap:
                    data = data.data
                else:
                    data = np.array(data.data)
        elif key_data in list(parent.attrs):
            data = parent.attrs[key_data]
            if isinstance(data, exdir.core.attribute.Attribute):
                data = dict(data)
        
        if 'data' not in locals():
            raise RuntimeError(f'{key} does not exist')
        
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

    def get(self, key, as_memmap=False):
        """Retrieve data.

        Parameters
        ----------
        key : :obj:`str`
            Key of the desired data. Nested keys should be separated by ``/``.
        as_memmap : :obj:`bool`, default: ``False``
            Keep NumPy memmap (from exdir files) instead of converting to arrays.
        
        Examples
        --------
        >>> rfile.get('/prod_1/charge')
        0
        >>> rfile.get('energy_scf')
        -12419.360138637763
        """
        key = self.clean_key(key)
        if self.ftype == 'exdir':
            data = self._get_from_exdir(key, as_memmap=as_memmap)
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
    
    def simplify_iter_data(self, data, data_key):
        """Simplify iterative data if possible.

        Checks contents of lists, tuples, and arrays to see if we can simplify.
        For example, if a list has only one element this method will return
        just the element.

        Some exceptions are used to maintain consistency. For example, atomic
        numbers will always remain a 1D array even if there is only one
        atom.

        This becomes important when adding parsed data to a file. During
        parsing, there is a possibility we could have one or more of the same
        property. We instead assume there will be multiple. Then as a
        postprocessing step we simplify cases where only one value was parsed.
        
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
                if data_key in ['atomic_numbers', 'entity_ids', 'comp_ids']:
                    return data
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
                if data_key in ['atomic_numbers', 'entity_ids', 'comp_ids']:
                    return data
                data = data[0].item()
        
        return data
    
    def _put_to_dict(self, key, data):
        """Add data to dictionary-like file.

        Parameters
        ----------
        key : :obj:`str`
            Where to put the data. Can be a nested key.
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
            data = self.simplify_iter_data(data, keys[-1])

        add_dic = {keys[-1]: data}
        for key in reversed(keys[:-1]):
            add_dic = {key: add_dic}

        self.File_ = combine_dicts(self.File_, add_dic)
    
    def _put_to_exdir(self, key, data):
        """Add data to an exdir group.

        Parameters
        ----------
        key : :obj:`str`
            Where to put the data. Can be a nested key.
        data : various
            Data to add to exdir file.
        """
        parent_key, data_key = self.split_key(key)
        group = self.get(parent_key)
        assert isinstance(group, exdir.core.exdir_file.File) \
            or isinstance(group, exdir.core.group.Group)
        
        # Custom handling of parsed data to ensure logical behavior.
        ndarray_to_list_keys = [
            'dipole_moment', 'periodic', 'periodic_cell'
        ]
        if data_key in ndarray_to_list_keys:
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
            data = self.simplify_iter_data(data, data_key)
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

    def put(self, key, data):
        """Put data to file in a specific location.

        Note that there is some data postprocessing using
        :meth:`~reptar.File.simplify_iter_data`.

        Parameters
        ----------
        key : :obj:`str`
            Key of the desired data (including parent). Nested keys should be
            separated by ``/``.
        data : ``obj``
            Data to add to file.
        """
        key = self.clean_key(key)
        if self.ftype == 'exdir':
            self._put_to_exdir(key, data)
        elif self.ftype == 'json' or self.ftype == 'npz':
            self._put_to_dict(key, data)
    
    def put_all(self, group_key, data, nested=False):
        """Adds all data from :obj:`dict` to group.

        This is just a loop over :meth:`~reptar.File.put` for each key-value pair.

        Parameters
        ----------
        group_key : :obj:`str`
            Key to the desired new group.
        data : :obj:`dict`
            Key-value pairs of data to add to group. For example, the
            ``parsed_info`` attribute after parsing a calculation.
        nested : :obj:`bool`, default: ``True``
            If``data`` contains one level of nested dictionaries. This is the
            case for ``parsed_info``.
        
        Returns
        -------
        """
        if nested:
            for cat_key in data.keys():
                for data_key in data[cat_key].keys():
                    value = data[cat_key][data_key]
                    self.put(f'{group_key}/{data_key}', value)
        else:
            for data_key in data.keys():
                value = data[data_key]
                self.put(f'{group_key}/{data_key}', value)
        
        return self.File_
    
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
    
    def create_group(self, key):
        """Initialize/create an exdir group with the specified key.

        Parameters
        ----------
        key : :obj:`str`
            Key of the desired data (including parent). Nested keys should be
            separated by ``/``.
        """
        return self.File_.create_group(key)
    
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
        self.put(f'{group_key}/md5', md5)
        
        try:
            md5_arrays = get_md5(self, group_key, only_arrays=True)
            self.put(f'{group_key}/md5_arrays', md5_arrays)
        except Exception:
            pass

        try:
            md5_structures = get_md5(self, group_key, only_structures=True)
            self.put(f'{group_key}/md5_structures', md5_structures)
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
