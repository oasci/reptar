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

import collections
from functools import reduce
import operator
import importlib
import itertools
import hashlib
import os
import numpy as np
from qcelemental import periodictable as ptable
from qcelemental.molparse.from_arrays import validate_and_fill_geometry
from qcelemental.exceptions import ValidationError
from .descriptors import get_center_of_mass
from .logger import ReptarLogger

log = ReptarLogger(__name__)


def get_files(path, expression, recursive=True):
    r"""Returns paths to all files in a given directory that matches a provided
    expression in the file name.

    Commonly used to find all files of a certain type, e.g., output or xyz
    files.

    Parameters
    ----------
    path : :obj:`str`
        Specifies the directory to search.
    expression : :obj:`str`
        Expression to be tested against all file names in ``path``.
    recursive : :obj:`bool`, default: ``True``
        Recursively find all files in all subdirectories.

    Returns
    -------
    :obj:`list` [:obj:`str`]
        All absolute paths to files matching the provided expression.
    """
    if path[-1] != "/":
        path += "/"
    if recursive:
        all_files = []
        for (dirpath, _, filenames) in os.walk(path):
            index = 0
            while index < len(filenames):
                if dirpath[-1] != "/":
                    dirpath += "/"
                filenames[index] = dirpath + filenames[index]
                index += 1
            all_files.extend(filenames)
        files = []
        for f in all_files:
            if expression in f:
                files.append(f)
    else:
        files = []
        for f in os.listdir(path):
            filename = os.path.basename(f)
            if expression in filename:
                files.append(path + f)
    return files


def atoms_by_element(atom_list):
    r"""Converts a list of atoms identified by their atomic number to their
    elemental symbol in the same order.

    Parameters
    ----------
    atom_list : :obj:`list` [:obj:`int`]
        Atomic numbers of atoms within a structure.

    Returns
    -------
    :obj:`list` [:obj:`str`]
        Element symbols of atoms within a structure.
    """
    return [ptable.to_symbol(z) for z in atom_list]


def atoms_by_number(atom_list):
    r"""Converts a list of atoms identified by their elemental symbol to their
    atomic number.

    Parameters
    ----------
    atom_list : :obj:`list` [:obj:`str`]
        Element symbols of atoms within a structure.

    Returns
    -------
    :obj:`list` [:obj:`int`]
        Atomic numbers of atoms within a structure.
    """
    return [ptable.to_atomic_number(symbol) for symbol in atom_list]


def parse_xyz(xyz_path):
    r"""Parses data from xyz file.

    An xyz file is data presented as consecutive xyz structures. The data could be
    three Cartesian coordinates for each atom, three atomic force vector
    components, or both coordinates and atomic forces in one line (referred to
    as extended xyz).

    Parameters
    ----------
    xyz_path : :obj:`str`
        Path to xyz file.

    Returns
    -------
    :obj:`tuple` [:obj:`list`]
        Parsed atoms (as element symbols :obj:`str`), comments, and data as
        :obj:`float` from string file.
    """
    Z, comments, data = [], [], []
    with open(xyz_path, "r", encoding="utf-8") as f:
        for _, line in enumerate(f):
            line = line.strip()
            if not line:
                # Skips blank lines
                pass
            else:
                line_split = line.split()
                if (
                    len(line_split) == 1
                    and float(line_split[0]) % int(line_split[0]) == 0.0
                ):
                    # Skips number of atoms line, adds comment line, and
                    # prepares next z and data item.
                    comment_line = next(f)
                    comments.append(comment_line.strip())
                    Z.append([])
                    data.append([])
                else:
                    # Grabs z and data information.
                    Z[-1].append(line_split[0])
                    data[-1].append([float(i) for i in line_split[1:]])
    return Z, comments, data


def get_md5(rfile, group_key, only_arrays=False, only_structures=False):
    r"""Creates MD5 hash for a group.

    Parameters
    ----------
    rfile : :obj:`reptar.File`
        A reptar file.
    group_key : :obj:`str`
        Key to the desired group.
    only_arrays : :obj:`bool`, default: ``False``
        Generate the MD5 hash using only arrays. This creates a data-centered
        MD5 that is not affected by data that are commonly added or
        changed. Defaults to ``False``.
    only_structures : :obj:`bool`, default: ``False``
        Generate the MD5 has with only ``atomic_numbers`` and ``geometry``
        if possible). This is more static than ``only_arrays`` and should be
        used to track sampling (i.e., ``r_prov_ids``).

    Returns
    -------
    :obj:`str`
        MD5 hash of a group.
    """
    md5_hash = hashlib.md5()
    # TODO: Figure out why different formats have different MD5s.

    if only_structures:
        try:
            Z = rfile.get(f"{group_key}/atomic_numbers")
            Z = Z.ravel()
            md5_hash.update(hashlib.md5(Z).digest())
        except Exception:
            pass
        try:
            R = rfile.get(f"{group_key}/geometry")
            R = R.ravel()
            md5_hash.update(hashlib.md5(R).digest())
        except Exception:
            pass
    else:
        keys = rfile.get_keys(group_key)
        for key in keys:
            if "md5" in key:
                continue
            d = rfile.get(f"{group_key}/{key}")
            if isinstance(d, np.ndarray):
                d = d.ravel()
                md5_hash.update(hashlib.md5(d).digest())
            else:
                if not only_arrays:
                    md5_hash.update(repr(d).encode())

    return md5_hash.hexdigest()


def gen_entity_ids(atoms_per_mol, num_mol, starting_idx=0, add_to=None):
    r"""Generates entity ids for a single species.

    Note that all of the atoms in each molecule must occur in the same order and
    be grouped together.

    Parameters
    ----------
    atoms_per_mol : :obj:`int`
        Number of atoms in the molecule.
    num_mol : :obj:`int`
        Number of molecules of this type in the system.
    starting_idx : :obj:`int`
        Number to start entity_id labels.
    add_to : :obj:`list`
        Entity ids to append new ids to.

    Returns
    -------
    :obj:`numpy.ndarray`
        Entity ids for a structure.
    """
    entity_ids = []
    for i in range(starting_idx, num_mol + starting_idx):
        entity_ids.extend([i for _ in range(0, atoms_per_mol)])

    if add_to is not None:
        if isinstance(add_to, np.ndarray):
            add_to = add_to.tolist()
        return np.array(add_to + entity_ids)

    return np.array(entity_ids)


def gen_comp_ids(label, num_mol, add_to=None):
    r"""Prepares the list of component ids for a system with only one species.

    Parameters
    ----------
    label : :obj:`int`
        Species label.
    num_mol : :obj:`int`
        Number of molecules of this type in the system.
    add_to : :obj:`list`
        Component ids to append new ids to.

    Returns
    -------
    :obj:`numpy.ndarray`
        Component ids for a structure.
    """
    comp_ids = [label for _ in range(0, num_mol)]
    if add_to is not None:
        if isinstance(add_to, np.ndarray):
            add_to = add_to.tolist()
        return np.array(add_to + comp_ids)
    return np.array(comp_ids)


def center_structures(Z, R):
    r"""Centers each structure's center of mass to the origin.

    Previously centered structures should not be affected by this technique.

    Parameters
    ----------
    Z : :obj:`numpy.ndarray`, ndim: ``1``
        Atomic numbers of the atoms in every structure.
    R : :obj:`numpy.ndarray`, ndim: ``3``
        Cartesian atomic coordinates of data set structures.

    Returns
    -------
    :obj:`numpy.ndarray`
        Centered Cartesian atomic coordinates.
    """
    # Masses of each atom in the same shape of R.
    if R.ndim == 2:
        R = np.array([R])
    n_atoms = R.shape[1]

    R -= np.repeat(get_center_of_mass(Z, R), n_atoms, axis=0).reshape(R.shape)

    if R.shape[0] == 1:
        R = R[0]

    return R


def combine_dicts(dict1, dict2):
    r"""Combine two dictionaries.

    Parameters
    ----------
    dict1 : :obj:`dict`
    dict2 : :obj:`dict`

    Returns
    -------
    :obj:`dict`
    """
    for k, v in dict2.items():
        if isinstance(v, collections.abc.Mapping):
            dict1[k] = combine_dicts(dict1.get(k, {}), v)
        else:
            dict1[k] = v
    return dict1


def find_parent_r_idxs(r_prov_specs, r_prov_specs_subset):
    r"""Find the structure indices of a parent r_prov_specs using a subset of
    the specifications.

    Useful for identifying structure indices when ``r_prov_specs_subset`` is in
    a different order.

    Parameters
    ----------
    r_prov_specs : :obj:`numpy.ndarray`, ndim: ``2``
        Structure provenance specifications.
    r_prov_specs_subset : :obj:`numpy.ndarray`, ndim: ``2``
        Specifications in no particular order.

    Returns
    -------
    :obj:`numpy.ndarray`, ndim: ``1``
        Structure indices of ``r_prov_specs_subset`` with respect to
        ``r_prov_specs``. Specifications that cannot be found will have
        ``NaN`` as their element.

    Example
    -------
    >>> r_prov_specs = np.array([[0, 0, 1], [0, 0, 2], [1, 1, 1]])
    >>> r_prov_specs_subset = np.array([[0, 0, 1], [1, 1, 1]])
    >>> find_parent_r_idxs(r_prov_specs, r_prov_specs_subset)
    [0, 2]

    """
    r_idxs = np.empty(r_prov_specs_subset.shape[0])
    r_idxs[:] = np.NaN
    # pylint: disable=consider-using-enumerate
    for i in range(len(r_idxs)):
        r_prov_spec = r_prov_specs_subset[i]
        r_idx = np.where((r_prov_spec == r_prov_specs).all(1))[0]
        if len(r_idx) == 1:
            r_idxs[i] = r_idx[0]
    return r_idxs.astype("int")


def gen_combs(sets, replacement=False):
    r"""Generate combinations from multiple sets.

    Parameters
    ----------
    sets : :obj:`list` or :obj:`tuple`, ndim: ``2``
        An iterable that contains multiple sets.
    replacement : :obj:`bool`, default: ``False``
        Allows repeated combinations in different order. If ``False``,
        ``(0, 1)`` and ``(1, 0)`` could be possible if there is overlap
        in the sets.

    Yields
    ------
    :obj:`tuple`
        Combination of one element per set in ``sets``.

    Examples
    --------
    >>> sets = ((0,) (1, 2), (1, 2, 3))
    >>> combs = gen_combs(sets)
    >>> for comb in combs:
    ...     print(comb)
    ...
    (0, 1, 2)
    (0, 1, 3)
    (0, 2, 3)
    """
    combs = itertools.product(*sets)
    # Excludes combinations that have repeats (e.g., (0, 0) and (1, 1. 2)).
    combs = itertools.filterfalse(lambda x: len(set(x)) < len(x), combs)
    # At this point, there are still duplicates in this iterator.
    # For example, (0, 1) and (1, 0) are still included.
    for comb in combs:
        # Sorts options is to avoid duplicate structures.
        # For example, if combination is (1, 0) the sorted version is not
        # equal and will not be included.
        if not replacement:
            if sorted(comb) != list(comb):
                continue
        yield comb


def chunk_iterable(iterable, n):
    r"""Chunk an iterable into ``n`` objects.

    Parameters
    ----------
    iterable : ``iterable``
        Iterable to chunk.
    n : :obj:`int`
        Size of each chunk.

    Yields
    ------
    :obj:`tuple`
        ``n`` objects.
    """
    iterator = iter(iterable)
    for first in iterator:
        yield tuple(itertools.chain([first], itertools.islice(iterator, n - 1)))


def exists_in_array(a_slice, array):
    r"""Check if ``a_slice`` exists in an ``array``.

    Parameters
    ----------
    a_slice : :obj:`numpy.ndarray`
        An example slice of ``array``'s first dimension to check.
    array : :obj:`numpy.ndarray`
        Array to check.

    Returns
    -------
    :obj:`bool`
        If ``row`` is present in ``array``.
    """
    ndim = int(array.ndim)
    exists = array == a_slice
    # For each additional dimension after 1 we check the last dimension
    # if they are all true.
    for _ in range(1, ndim):
        exists = exists.all(axis=-1)
    # At the very end, we will have a 1D array. If any are True, then a_slice
    # exists in array
    return exists.any()


def dict_iterator(dictionary):
    r"""Iterate over nested dictionary.

    Parameters
    ----------
    dictionary : :obj:`dict`
        An arbitrarily nested dictionary.

    Yields
    ------
    :obj:`tuple`
        Keys specifying the walk through a dictionary and the last item as the value.
    """
    for k, v in dictionary.items():
        if isinstance(v, dict):
            # If value is dict then iterate over all its values
            for pair in dict_iterator(v):
                yield k, *pair
        else:
            # If value is not dict type then yield the value
            yield k, v


def get_nested_key(dictionary, keys):
    r"""Access a nested object in a dictionary by key sequence.

    Parameters
    ----------
    dictionary : :obj:`dict`
        Dictionary to get data from.
    keys : :obj:`list`
        Keys that lead to data in ``dictionary``.

    Notes
    -----
    Some code here is from `this Stack Overflow answer
    <https://stackoverflow.com/a/14692747>`__ by
    `Martijn Pieters <https://stackoverflow.com/users/100297/martijn-pieters>`__,
    licensed under `CC BY-SA 4.0 <https://creativecommons.org/licenses/by-sa/4.0/>`__.
    """
    return reduce(operator.getitem, keys, dictionary)


def add_nested_key(dictionary, keys, data):
    r"""Set a value in a nested dictionary by key sequence.

    Parameters
    ----------
    dictionary : :obj:`dict`
        Dictionary to add data to.
    keys : :obj:`list`
        Keys that lead to data in ``dictionary``.
    data : ``obj``
        Data to add under ``keys`` in ``dictionary``.

    Returns
    -------
    :obj:`dict`
        ``dictionary`` with the added data.

    Notes
    -----
    Some code here is from `this Stack Overflow answer
    <https://stackoverflow.com/a/14692747>`__ by
    `Martijn Pieters <https://stackoverflow.com/users/100297/martijn-pieters>`__,
    licensed under `CC BY-SA 4.0 <https://creativecommons.org/licenses/by-sa/4.0/>`__.
    """
    get_nested_key(dictionary, keys[:-1])[keys[-1]] = data
    return dictionary


def remove_nested_key(dictionary, keys):
    r"""Delete data in a nested dictionary.

    Parameters
    ----------
    dictionary : :obj:`dict`
        Dictionary to delete data from.
    keys : :obj:`list`
        Keys that lead to data in ``dictionary``.

    Returns
    -------
    :obj:`dict`
        ``dictionary`` with data removed.

    Notes
    -----
    Some code here is from `this Stack Overflow answer
    <https://stackoverflow.com/a/14692747>`__ by
    `Martijn Pieters <https://stackoverflow.com/users/100297/martijn-pieters>`__,
    licensed under `CC BY-SA 4.0 <https://creativecommons.org/licenses/by-sa/4.0/>`__.
    """
    del get_nested_key(dictionary, keys[:-1])[keys[-1]]
    return dictionary


def validate_geometry(R, tooclose=0.1):
    """Checks if any atoms would clash.

    Parameters
    ----------
    R : :obj:`numpy.ndarray`, ndim: ``3``
        Atomic coordinates.
    too_close : :obj:`float`, default: ``0.1``
        Threshold for atom distances for validated geometries. All atoms must be at
        least this far away.

    Returns
    -------
    :obj:`int`
        Number of invalid structures.
    :obj:`numpy.ndarray`
        (shape: ``R.shape[0]``) If structures in R are valid.
    """
    if R.ndim == 2:
        R = R[None, ...]
    is_valid = np.full((R.shape[0],), True, dtype=np.bool)
    for i in range(R.shape[0]):
        try:
            validate_and_fill_geometry(geom=R[i], tooclose=tooclose)
        except ValidationError:
            is_valid[i] = False

    n_invalid = (~is_valid).sum()

    return n_invalid, is_valid


def prep_array(rfile, key, shape, dtype=None, fill_with=np.nan):
    """Prepares array in files by first trying to ``get()`` data, but if it does not
    exist then initializes with ``put()``.

    Parameters
    ----------
    rfile : :obj:`reptar.File`
        File object to prepare data in.
    key : :obj:`str`
        Key to desired data to prepare in ``rfile``.
    shape : :obj:`tuple`
        Desired shape of the array if it does not already exist.
    dtype : :obj:`str`, default: ``None``
        Desired NumPy data type. If ``None`` then NumPy will try to use a default
        data type that can represent the values.
    fill_with : ``various``, default: ``NaN``
        Fill the array with ``fill_with`` if we are initializing it.

    Returns
    -------
    ``various``
        Requested data.
    """
    try:
        data = rfile.get(key)
    except RuntimeError as e:
        if "does not exist" in str(e):
            data = np.full(shape, fill_with, dtype=dtype)
            rfile.put(key, data)
            data = rfile.get(key)
    return data


# pylint: disable=invalid-name
def prep_group_opt(
    rfile,
    Z_key,
    R_key,
    dest_key,
    Z_opt_label="atomic_numbers",
    conv_opt_label="conv_opt",
    R_opt_label="geometry",
    E_opt_label="energy_ele",
):
    """Prepare group for geometry optimization calculations in a reptar file.

    Parameters
    ----------
    rfile : :obj:`reptar.File`
        File to prepare a group for optimization-like data.
    Z_key : :obj:`str`
        Source key for atomic numbers of structures to optimize.
    R_key : :obj:`str`
        Source key for Cartesian coordinates of structures to optimize.
    dest_key : :obj:`str`
        Key to create new group for optimization. All data will be stored here.
    Z_opt_label : :obj:`str`, default: ``atomic_numbers``
        Label to store atomic numbers of new optimization group.
    conv_opt_label : :obj:`str`, default: ``conv_opt``
        Label to store boolean flag if the geometry optimization was successful.
    R_opt_label : :obj:`str`, default: ``geometry``
        Label to store the last structure from the optimization. If ``converged_opt``
        is ``True``, this is the optimized geometry. If ``False``, then it is the last
        structure before the optimization terminated.
    E_opt_key : :obj:`str`, default: ``energy_ele``
        Label to store the last electronic energy after the optimization.

    Returns
    -------
    :obj:`tuple`
        ``Z_opt_key``, ``conv_opt_key``, ``R_opt_key``, and ``E_opt_key``.
    :obj:`tuple`
        ``Z_opt``, ``conv_opt``, ``R_opt``, and ``E_opt`` arrays.
    """
    log.debug("Retrieving initial configurations")
    Z = rfile.get(Z_key)
    R = rfile.get(R_key)

    log.debug("Initializing optimization data")
    try:
        rfile.get(dest_key)
    except RuntimeError as e:
        if " does not exist" in str(e):
            rfile.create_group(dest_key)
        else:
            raise RuntimeError from e

    Z_opt_key = os.path.join(dest_key, Z_opt_label)
    rfile.put(Z_opt_key, Z)

    R_opt_key = os.path.join(dest_key, R_opt_label)
    R_opt = prep_array(rfile, R_opt_key, R.shape, dtype="float64")

    conv_opt_key = os.path.join(dest_key, conv_opt_label)
    conv_opt = prep_array(
        rfile, conv_opt_key, R.shape[0], dtype="bool", fill_with=False
    )

    E_opt_key = os.path.join(dest_key, E_opt_label)
    E_opt = prep_array(rfile, E_opt_key, R.shape[0], dtype="float64")

    return (Z_opt_key, conv_opt_key, R_opt_key, E_opt_key), (Z, conv_opt, R_opt, E_opt)


# pylint: disable=invalid-name
def prep_group_engrad(
    rfile,
    R_key,
    dest_key,
    E_label="energy_ele",
    G_label="grads",
):
    """Prepare group for geometry optimization calculations in a reptar file.

    Parameters
    ----------
    rfile : :obj:`reptar.File`
        File to prepare a group for optimization-like data.
    R_key : :obj:`str`
        Source key for Cartesian coordinates of structures to optimize.
    dest_key : :obj:`str`
        Key to create new group for optimization. All data will be stored here.
    E_key : :obj:`str`, default: ``energy_ele``
        Label to store the electronic energy.
    G_key : :obj:`str`, default: ``grads``
        Label to store the atomic gradient.

    Returns
    -------
    :obj:`tuple`
        ``E_key`` and ``G_key``.
    :obj:`tuple`
        ``E``, and ``G`` arrays.
    """
    R = rfile.get(R_key)

    log.debug("Initializing energy and gradient data")
    try:
        rfile.get(dest_key)
    except RuntimeError as e:
        if " does not exist" in str(e):
            rfile.create_group(dest_key)
        else:
            raise RuntimeError from e

    E_key = os.path.join(dest_key, E_label)
    E = prep_array(rfile, E_key, R.shape[0], dtype="float64", fill_with=np.nan)

    G_key = os.path.join(dest_key, G_label)
    G = prep_array(rfile, G_key, R.shape, dtype="float64", fill_with=np.nan)

    return (E_key, G_key), (E, G)


def get_obj_from_string(import_string):
    """Retrieves a function object based on an import string and object name.

    Parameters
    ----------
    import_string : :obj:`str`
        The import string representing the module containing the desired object.

    Returns
    -------
    ``object``
        The object identified by the import string.
    """
    module_name, obj_name = import_string.rsplit(".", 1)
    module = importlib.import_module(module_name)
    obj = getattr(module, obj_name)
    return obj


def common_elements(arr1: np.ndarray, arr2: np.ndarray) -> np.ndarray:
    r"""Provides only elements that are shared between two arrays.

    Parameters
    ----------
    arr1
        First array.
    arr2
        Second array.

    Returns
    -------
    :obj:`np.ndarray`
        Array containing only shared elements between ``arr1`` and ``arr2``.
    """
    u, c = np.unique(np.concatenate((arr1, arr2)), return_counts=True)
    return u[c > 1]
