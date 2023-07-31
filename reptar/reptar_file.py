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

import os
import json
import exdir
import numpy as np
import yaml
import zarr
from .utils import combine_dicts, dict_iterator, get_md5, remove_nested_key
from .logger import ReptarLogger

log = ReptarLogger(__name__)


class File:
    r"""Create, store, and access data from a variety of formats."""

    # pylint: disable=unnecessary-dunder-call

    def __init__(self, file_path, mode="r", plugins=None, from_dict=None):
        """
        Parameters
        ----------
        file_path : :obj:`str`
            Path to a file supported by reptar. If it does not exist, then one will
            be created if possible.
        mode : :obj:`str`, default: ``'r'``
            A file mode string that defines the read/write behavior.
        plugins : :obj:`list`, default: ``None``
            A list of instantiated exdir plugins.
        from_dict : :obj:`dict`, default: ``None``
            Load data from a dictionary.
        """
        if "w" in mode:
            allow_remove = True
        elif mode in ["a", "r"]:
            allow_remove = False
        self.allow_remove = allow_remove
        self.plugins = plugins
        if from_dict is None:
            self._from_path(file_path, mode, allow_remove, plugins)
        else:
            self._from_dict(file_path, from_dict, mode, allow_remove, plugins)

    def _from_path(self, file_path, mode, allow_remove, plugins):
        r"""Populates the File object from a file path.

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
        log.debug("Opening file from path")
        exists = os.path.exists(file_path)
        _, f_ext = os.path.splitext(file_path)

        if f_ext == ".exdir":
            File_ = exdir.File(file_path, mode, allow_remove, plugins)
        elif f_ext == ".json":
            if exists:
                with open(file_path, "r", encoding="utf-8") as f:
                    File_ = json.load(f)
            else:
                File_ = {}
        elif f_ext == ".npz":
            if exists:
                File_ = dict(np.load(file_path, allow_pickle=True))
                # Since everything is stored in arrays, we clean up array data.
                for k, v in File_.items():
                    if self._is_iter(v):
                        v = self.simplify_iter_data(v, k.split("/")[-1])
                        File_[k] = v
            else:
                File_ = {}
        elif f_ext == ".zarr":
            File_ = zarr.open(file_path, mode)
        else:
            raise TypeError(f"{f_ext} is not supported.")

        self.fpath = os.path.abspath(file_path)
        self.ftype = f_ext[1:]
        self.fmode = mode
        self.File_ = File_

    def __getstate__(self):
        r"""Having File_ makes this not serializable."""
        state = self.__dict__.copy()
        del state["File_"]
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._from_path(self.fpath, self.fmode, self.allow_remove, self.plugins)

    def _from_dict(self, file_path, group_dict, mode, allow_remove, plugins):
        r"""Populates the File object from a dictionary.

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
        log.debug("Opening file from dictionary")
        _, f_ext = os.path.splitext(file_path)

        if f_ext == ".exdir":
            self.File_ = exdir.File(file_path, mode, allow_remove, plugins)
        elif f_ext in (".json", ".npz"):
            self.File_ = {}
        elif f_ext == ".zarr":
            self.File_ = zarr.open(file_path, mode)
        self.fpath = os.path.abspath(file_path)
        self.ftype = f_ext[1:]
        self.fmode = mode

        for item in dict_iterator(group_dict):
            key = "/".join(item[:-1])
            data = item[-1]
            self.put(key, data)

    @staticmethod
    def clean_key(key):
        r"""Clean key and remove any common mistakes.

        This corrects instances of ``//`` and enforces instances where the start of the
        key contains ``.``.

        Parameters
        ----------
        key : :obj:`str`
            Key of the desired data.

        Returns
        -------
        :obj:`str`
            Cleaned key.
        """
        log.debug("Cleaning key")
        log.debug("Original key: %s", key)
        if "//" in key:
            key = key.replace("//", "/")
        while key[0] == ".":
            log.debug("Removing . at beginning")
            key = key[1:]
        log.debug("Cleaned key: %s", key)
        return key

    @staticmethod
    def split_key(key):
        r"""Split the key into a parent and data key.

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
        log.debug("Splitting key")
        if key[0] != "/":
            key = "/" + key
        key_split = key.rsplit("/", 1)
        if key_split[0] == "":
            key_split[0] = "/"
        log.debug("Parent key: %s", key_split[0])
        log.debug("Data key: %s", key_split[1])
        return key_split

    def _get_from_dict(self, key):
        r"""Get data from dictionary-like file (e.g., json, npz).

        Parameters
        ----------
        key : :obj:`str`
            Key of the desired data. Nested keys should be separated by ``/``.

        Returns
        -------
        ``various``
            Requested data from a dictionary source.
        """
        log.debug("Getting data from dictionary")

        keys = key.split("/")
        if keys[0] == "":
            del keys[0]

        data = self.File_[keys[0]]
        for k in keys[1:]:
            data = data[k]
        if isinstance(data, list):
            data_array = np.array(data)
            if data_array.dtype != "O":
                data = data_array
        return data

    def _get_from_exdir(self, key, as_memmap=False):
        r"""Get data from exdir file.

        Parameters
        ----------
        key : :obj:`str`
            Key of the desired data. Nested keys should be separated by ``/``.
        as_memmap : :obj:`bool`, default: ``False``
            Keep NumPy memmap instead of converting to arrays.

        Returns
        -------
        ``various``
            Requested data from a exdir source.
        """
        log.debug("Getting data from exdir")

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

        if "data" not in locals():
            raise RuntimeError(f"{key} does not exist")

        return data

    def _get_from_zarr(self, key):
        r"""Get data from zarr file.

        Parameters
        ----------
        key : :obj:`str`
            Key of the desired data. Nested keys should be separated by ``/``.

        Returns
        -------
        ``various``
            Requested data from a zarr source.
        """
        log.debug("Getting data from zarr")

        # Check if key is attribute
        key_parent, key_data = self.split_key(key)
        attr_keys = list(self.File_[key_parent].attrs.keys())
        if key_data in attr_keys:
            data = self.File_[key_parent].attrs[key_data]
        else:
            # Should be array.
            try:
                data = self.File_[key]
            except KeyError as e:
                raise RuntimeError(f"{key} does not exist") from e

        return data

    def get_keys(self, group_key):
        r"""A list of keys in a group.

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
        if self.ftype == "exdir":
            keys.extend(list(sorted(group.attrs.keys())))
            keys.extend(
                list(
                    sorted(
                        key
                        for key in group.keys()
                        if not isinstance(group[key], exdir.core.raw.Raw)
                    )
                )
            )
        elif self.ftype in ("json", "npz"):
            keys.extend(list(group.keys()))
        return keys

    def get(self, key, as_memmap=False, missing_is_none=False):
        r"""Retrieve data.

        Parameters
        ----------
        key : :obj:`str`
            Key of the desired data. Nested keys should be separated by ``/``.
        as_memmap : :obj:`bool`, default: ``False``
            Keep NumPy memmap (from exdir files) instead of converting to arrays.
        missing_is_none : :obj:`bool`, default: ``False``
            Catch the ``RuntimeError`` and return ``None`` if the key does not
            exits.

        Examples
        --------
        >>> rfile.get('/prod_1/charge')
        0
        >>> rfile.get('energy_scf')
        -12419.360138637763
        """
        key = self.clean_key(key)

        if key == "/":
            log.debug("Requested root (i.e., '/'). Returning whole file.")
            return self.File_

        try:
            if self.ftype == "exdir":
                data = self._get_from_exdir(key, as_memmap=as_memmap)
            elif self.ftype in ("json", "npz"):
                data = self._get_from_dict(key)
            elif self.ftype == "zarr":
                data = self._get_from_zarr(key)
        except RuntimeError as e:
            if "does not exist" in str(e) and missing_is_none:
                data = None
            else:
                raise
        return data

    @staticmethod
    def _is_iter(data):
        r"""If data is iterative (i.e., array, list, or tuple).

        Returns
        -------
        :obj:`bool`
            If the data is iterative.
        """
        if isinstance(data, (list, tuple, np.ndarray, zarr.core.Array)):
            return True
        return False

    @staticmethod
    def simplify_iter_data(data, data_key):
        r"""Checks contents of lists, tuples, and arrays to see if we can simplify.

        Parameters
        ----------
        data : :obj:`numpy.ndarray`, :obj:`list`, :obj:`tuple`
            An iterative data object.
        data_key : :obj:`str`, default: ``None``
            Key for the data. There are some hard-coded checks and behaviors.

        Returns
        -------
        ``obj``
            Simplified data (if possible).

        Notes
        -----
        There are some circumstances where an iterable is not simplified.
        For example, if there is only one atom and we are trying to simplify
        ``atomic_numbers`` (i.e., ``[7]``) it is kept as an iterable.

        This routine is important after sampling. As the default behavior is to assume
        multiple values will be parsed. This routine is used in the end to clean up
        those data.

        The following cases -> behaviors are enforced on ``data_keys``.

        - Length of iterable is one -> Keep as iterable.
            - ``atomic_numbers``, ``entity_ids``, ``comp_ids``
        """
        log.debug("Checking if data can be simplified")
        # If data is a list we check the types of data it contains.
        # If all of them are strings, we do not convert to array.
        if not isinstance(data, np.ndarray):
            # If there is only one item in the list we can likely
            # just store it as an attribute. However, we have
            # to be careful when the single item is an array.
            if len(data) == 1:
                if data_key in ["atomic_numbers", "entity_ids", "comp_ids"]:
                    log.debug("Key is reserved for arrays")
                    return data
                log.debug("Returning the single item")
                data = data[0]
            # If all the items are not all strings then we make
            # the dataset.
            elif not all(isinstance(i, str) for i in data):
                data_array = np.array(data)
                # Sometimes we have a list of objects that cannot be easily
                # converted to numpy data types. Numpy will use an object
                # dtype, so we only convert to array if its not an 'O' dtype.
                if data_array.dtype != "O":
                    data = data_array
            # At this point only data that contains all strings
            # should be left. We put these as attributes (i.e.,
            # comp_ids)
            else:
                pass
        # Handling arrays.
        else:
            log.debug("Data type is array")
            # If only one item we store it as a value.
            if data.shape == (1,):
                log.debug("Shape is (1, )")
                if data_key in ["atomic_numbers", "entity_ids", "comp_ids"]:
                    log.debug("Key is reserved for arrays")
                    return data
                log.debug("Returning the single item")
                data = data[0].item()

        return data

    def _put_to_dict(self, key, data):
        r"""Add data to dictionary-like file.

        Parameters
        ----------
        key : :obj:`str`
            Where to put the data. Can be a nested key.
        data : ``various``
            Data to add to file.
        """
        if isinstance(data, np.ndarray):
            if "U" in str(data.dtype):
                data = data.tolist()

        keys = key.split("/")
        if keys[0] == "":
            del keys[0]

        is_iter = self._is_iter(data)
        if is_iter:
            data = self.simplify_iter_data(data, keys[-1])

        add_dic = {keys[-1]: data}
        for key in reversed(keys[:-1]):  # pylint: disable=redefined-argument-from-local
            add_dic = {key: add_dic}

        self.File_ = combine_dicts(self.File_, add_dic)

    # pylint: disable-next=too-many-branches
    def _put_to_exdir(self, key, data):
        r"""Add data to an exdir group.

        Parameters
        ----------
        key : :obj:`str`
            Where to put the data. Can be a nested key.
        data : ``various``
            Data to add to exdir file.
        """
        parent_key, data_key = self.split_key(key)

        # Handle nested creation of groups.
        # Suppose we want to put "/group_1/nested_group/data_key" before group_1 or
        # nested_group exists. We have to create these groups first.
        # create_group
        try:
            group = self.get(parent_key)
        except KeyError as e:
            if "No such object:" not in str(e):
                raise
            group = self.create_group(parent_key)
        except RuntimeError as e:
            if " does not exist" not in str(e):
                raise
            group = self.create_group(parent_key)
        assert isinstance(group, (exdir.core.exdir_file.File, exdir.core.group.Group))

        # Custom handling of parsed data to ensure logical behavior.
        ndarray_to_list_keys = ["dipole_moment", "periodic", "periodic_cell"]
        # pylint: disable-next=no-else-return
        if data_key in ndarray_to_list_keys:
            if isinstance(data, np.ndarray):
                data = data.tolist()
            group.attrs[data_key] = data
            return None
        elif data_key == "wall_potential":
            group.attrs[data_key] = data
            return None

        # If data is not array/list or array/list of one dimension and
        # length we make it an exdir attribute.
        # We make anything else an exdir dataset.
        is_iter = self._is_iter(data)
        if is_iter:
            data = self.simplify_iter_data(data, data_key)
            if isinstance(data, np.ndarray):
                store_as = "dset"
            else:
                store_as = "attr"
        # Handles all non iterable data.
        else:
            store_as = "attr"

        if store_as == "attr":
            group.attrs[data_key] = data
        elif store_as == "dset":
            try:
                group.create_dataset(data_key, data=data)
            except RuntimeError:
                # If data set already exists, exdir will throw a RuntimeError.
                # We remove this dataset and try again.
                group.__delitem__(data_key)
                group.create_dataset(data_key, data=data)
        return None

    def _put_to_zarr(self, key, data):
        r"""Add data to an zarr group.

        Parameters
        ----------
        key : :obj:`str`
            Where to put the data. Can be a nested key.
        data : ``various``
            Data to add to zarr file.
        """
        parent_key, data_key = self.split_key(key)

        # Handle nested creation of groups.
        # Suppose we want to put "/group_1/nested_group/data_key" before group_1 or
        # nested_group exists. We have to create these groups first.
        try:
            log.debug("Getting group: %s", parent_key)
            group = self.File_[parent_key]
        except KeyError:
            log.debug("Group not found")
            group = self.create_group(parent_key)
        assert isinstance(group, (zarr.hierarchy.Group))

        # Custom handling of parsed data to ensure logical behavior.
        if data_key == "wall_potential":
            group.attrs[data_key] = data
            return None

        # If data is not array/list or array/list of one dimension and
        # length we make it an exdir attribute.
        # We make anything else an exdir dataset.
        is_iter = self._is_iter(data)
        if is_iter:
            data = self.simplify_iter_data(data, data_key)
            if isinstance(data, (np.ndarray, zarr.core.Array)):
                store_as = "array"
            else:
                store_as = "attr"
        # Handles all non iterable data.
        else:
            store_as = "attr"

        if store_as == "attr":
            group.attrs[data_key] = data
        elif store_as == "array":
            group[data_key] = data
        return None

    def put(self, key, data, with_md5_update=False):
        r"""Put data to file in a specific location.

        Note that there is some data postprocessing using
        :meth:`~reptar.File.simplify_iter_data`.

        Parameters
        ----------
        key : :obj:`str`
            Key of the desired data (including parent). Nested keys should be
            separated by ``/``.
        data : ``obj``
            Data to add to file.
        with_md5_update : :obj:`bool`, default: ``False``
            Update MD5 hashes after putting data.
        """
        log.debug("Putting data with key %s", key)
        key = self.clean_key(key)
        if self.ftype == "exdir":
            self._put_to_exdir(key, data)
        elif self.ftype in ("json", "npz"):
            self._put_to_dict(key, data)
        elif self.ftype == "zarr":
            self._put_to_zarr(key, data)

        # Update MD5 of group
        if with_md5_update:
            group_key, _ = self.split_key(key)
            self.update_md5(group_key)

    def put_all(self, group_key, data, nested=False):
        r"""Adds all data from :obj:`dict` to group.

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
        """
        if nested:
            for cat_key in data.keys():
                for data_key in data[cat_key].keys():
                    value = data[cat_key][data_key]
                    self.put(f"{group_key}/{data_key}", value)
        else:
            for data_key in data.keys():
                value = data[data_key]
                self.put(f"{group_key}/{data_key}", value)

        return self.File_

    def _remove_dict(self, key):
        r"""Delete dictionary data under ``key``.

        Parameters
        ----------
        key : :obj:`str`
            Key of the desired data (including parent). Nested keys should be
            separated by ``/``.
        """
        # We know that npz cannot be nested; so we directly delete it.
        if self.ftype == "npz":
            del self.File_[key.replace("/", "")]

        # Any other file type we assume could be nested.
        remove_nested_key(self.File_, [k for k in key.split("/") if k != ""])

    def _remove_exdir(self, key):
        r"""Delete exdir data under ``key``.

        Parameters
        ----------
        key : :obj:`str`
            Key of the desired data (including parent). Nested keys should be
            separated by ``/``.
        """
        try:
            # Handles datasets
            self.File_.__delitem__(key)
        except KeyError as e:
            if "No such object: " not in str(e):
                raise

            # Handles attributes
            group_key, attr_key = self.split_key(key)
            if group_key[0] == "/":
                group_key = group_key[1:]
            yaml_path = os.path.join(self.fpath, group_key, "attributes.yaml")

            with open(yaml_path, "r", encoding="utf-8") as f:
                attrs = yaml.safe_load(f)

            del attrs[attr_key]

            with open(yaml_path, "w", encoding="utf-8") as f:
                yaml.dump(attrs, f)

    def _remove_zarr(self, key):
        r"""Delete zarr data under ``key``.

        Parameters
        ----------
        key : :obj:`str`
            Key of the desired data (including parent). Nested keys should be
            separated by ``/``.
        """
        try:
            # Handles datasets
            self.File_.__delitem__(key)
            log.debug("Deleted array key: %s", key)
        except KeyError:
            log.debug("Key is not an array. Must be an attribute.")

            # Handles attributes
            group_key, attr_key = self.split_key(key)
            if group_key[0] == "/":
                group_key = group_key[1:]
            json_path = os.path.join(self.fpath, group_key, ".zattrs")
            log.debug("Loading attribute file at %s", json_path)

            with open(json_path, "r", encoding="utf-8") as f:
                attrs = json.load(f)

            del attrs[attr_key]

            json_string = json.dumps(attrs)

            with open(json_path, "w", encoding="utf-8") as f:
                f.write(json_string)

    def remove(self, key, with_md5_update=False):
        r"""Delete data under ``key``.

        Parameters
        ----------
        key : :obj:`str`
            Key of the desired data (including parent). Nested keys should be
            separated by ``/``.
        with_md5_update : :obj:`bool`, default: ``False``
            Update MD5 hashes after putting data.
        """
        key = self.clean_key(key)

        if self.ftype == "exdir":
            self._remove_exdir(key)
        elif self.ftype in ("json", "npz"):
            self._remove_dict(key)
        elif self.ftype == "zarr":
            self._remove_zarr(key)

        # Update MD5 of group
        if with_md5_update:
            group_key, _ = self.split_key(key)
            self.update_md5(group_key)

    def copy(self, source_key, dest_key, with_md5_update=False):
        r"""Copy data from a source to a destination.

        Parameters
        ----------
        source_key : :obj:`str`
            Key of the data to copy.
        dest_key : :obj:`str`
            Where to copy the data do.
        """
        self.put(dest_key, self.get(source_key), with_md5_update)

    def create_group(self, key):
        r"""Initialize/create a group in heigherarcal files with the specified key.
        This can handle creating nested groups.

        Parameters
        ----------
        key : :obj:`str`
            Key of the desired data (including parent). Nested keys should be
            separated by ``/``.

        Returns
        -------
        :obj:`exdir.core.Group`
            The newly created group.
        """
        log.debug("Creating %s group with key %s", self.ftype, key)
        if self.ftype == "exdir":
            self.File_.create_group(key)
        elif self.ftype == "zarr":
            self.File_.create_group(key)
        return self.get(key)

    def as_dict(self, group_key):
        r"""Get a group as a dictionary.

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
        if self.ftype == "exdir":
            group_dict = {}
            data_keys = self.get_keys(group_key)
            get_keys = [f"{group_key}/{d_key}" for d_key in data_keys]
            for d_key, g_key in zip(data_keys, get_keys):
                group_dict[d_key] = self.get(g_key)
            return group_dict

        # json or npz
        return group

    def update_md5(self, group_key):
        r"""Update all possible MD5 hashes of a specific group.

        Parameters
        ----------
        group_key : :obj:`str`
            Desired group.
        """
        md5 = get_md5(self, group_key)
        self.put(f"{group_key}/md5", md5, with_md5_update=False)

        try:
            md5_arrays = get_md5(self, group_key, only_arrays=True)
            self.put(f"{group_key}/md5_arrays", md5_arrays, with_md5_update=False)
        except Exception:
            pass

        try:
            md5_structures = get_md5(self, group_key, only_structures=True)
            self.put(
                f"{group_key}/md5_structures", md5_structures, with_md5_update=False
            )
        except Exception:
            pass

    def save(self, json_prettify=True):
        r"""Saves non-exdir files.

        Parameters
        ----------
        json_prettify : :obj:`bool`
            Indents JSON objects if True. If false the JSON file is only one
            line.
        """
        assert self.fmode == "w"
        # pylint: disable=import-outside-toplevel
        from cclib.io.cjsonwriter import JSONIndentEncoder, NumpyAwareJSONEncoder

        if self.ftype == "json":
            json_dict = self.File_

            if json_prettify:
                json_string = json.dumps(
                    json_dict, cls=JSONIndentEncoder, sort_keys=True, indent=4
                )
            else:
                json_string = json.dumps(
                    json_dict, cls=NumpyAwareJSONEncoder, sort_keys=True
                )
            with open(self.fpath, "w", encoding="utf-8") as f:
                f.write(json_string)
        elif self.ftype == "npz":
            npz_dict = self.File_
            np.savez_compressed(self.fpath, **npz_dict)
