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

import itertools
import numpy as np
from random import randrange, choice
from . import Saver
from .utils import center_structures as get_center_structures
from .utils import get_md5, gen_combs, exists_in_array, chunk_iterable
from .periodic import Cell
from . import __version__ as reptar_version


def entity_mask_gen(entity_ids, entities):
    """Generate an atom mask for a single entity.

    Parameters
    ----------
    entity_ids : :obj:`numpy.ndarray`
        ``entity_ids`` for the data you are going to mask.
    entities : ``iterable``
        A collection of ``entity_ids`` that we will mask individually.

    Yields
    ------
    :obj:`numpy.ndarray`
        Atom mask for some data such as ``Z``, ``R``, ``E``, ``G``, etc.

    Examples
    --------
    Suppose we want to slice each entity's structures from an array ``r``.

    >>> r = np.array(
    ...     [[ 1.62983581, -5.72097814,  3.33683543],
    ...      [ 1.18517275, -4.99524281,  2.8550542 ],
    ...      [ 2.54530238, -5.47454298,  3.2638527 ],
    ...      [ 3.29853071,  5.74803325, -5.55958208],
    ...      [ 2.39044514,  5.71619898, -5.6366701 ],
    ...      [ 3.44827586,  6.41505464, -4.85337475],
    ...      [-4.89465455, -0.18369004,  2.24627128],
    ...      [-5.57558564,  0.09906561,  1.66343969],
    ...      [-4.52877354,  0.7259792 ,  2.33425484]]
    ... )
    >>> entity_ids = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
    >>> mask_gen = entity_mask_gen(entity_ids, [0, 2, 1])
    >>> next(mask_gen)
    array([ True,  True,  True, False, False, False, False, False, False])
    >>> next(mask_gen)
    array([False, False, False, False, False, False,  True,  True,  True])
    >>> next(mask_gen)
    array([False, False, False,  True,  True,  True, False, False, False])

    """
    # Try to catch the error of passing a single number instead of iterable.
    try:
        len(entities)
    except TypeError as e:
        if "has no len()" in str(e):
            entities = [entities]
        else:
            raise
    for entity_id in entities:
        yield (entity_ids == entity_id)


def r_from_entities(R, entity_ids, entities):
    """Slice a geometries containing each entity in the same order as
    ``entities``.

    Parameters
    ----------
    R : :obj:`numpy.ndarray`, ndim: ``2``
        Coordinates of a single structure.
    entity_ids : :obj:`numpy.ndarray`
        Entity IDs of ``R``.
    entities : ``iterable``
        A collection of ``entity_ids`` that we will mask individually.
    """
    r_sel = []
    for entity_mask in entity_mask_gen(entity_ids, entities):
        r_sel.extend(R[entity_mask])
    return np.array(r_sel)


def _add_structures_to_R(R, n_to_add, n_atoms):
    R_to_add = np.empty((n_to_add, n_atoms, 3))
    R_to_add[:] = np.nan
    if R is None:
        idx_sel = 0  # Starting index of the newly sampled structures.
        R = R_to_add
    else:
        idx_sel = len(R)
        R = np.concatenate((R, R_to_add), axis=0)
    return idx_sel, R


def _add_structures_to_E(E, n_to_add):
    E_to_add = np.empty(n_to_add)
    E_to_add[:] = np.nan
    if E is None:
        E = E_to_add
    else:
        E = np.concatenate((E, E_to_add))
    return E


def _add_structures_to_G(G, n_to_add, n_atoms):
    G_to_add = np.empty((n_to_add, n_atoms, 3))
    G_to_add[:] = np.nan
    if G is None:
        G = G_to_add
    else:
        G = np.concatenate((G, G_to_add), axis=0)
    return G


def _add_structures_to_r_prov_specs(r_prov_specs, n_to_add, comp_labels):
    r_prov_specs_to_add = np.empty((n_to_add, 2 + len(comp_labels)), dtype=np.uint32)
    if r_prov_specs is None:
        r_prov_specs = r_prov_specs_to_add
    else:
        r_prov_specs = np.concatenate((r_prov_specs, r_prov_specs_to_add), axis=0)
    return r_prov_specs


def _generate_structure_samples(quantity, structure_idxs, entity_ids_samples):
    """Randomly generates structure and entity selections.

    The prov ID for the structure is not included here and must be
    inserted later.

    Parameters
    ----------
    quantity : :obj:`int` or :obj:`float`
        Number of structures to sample from the structure set. For example,
        ``'100'``, ``452``, or even ``'all'``.
    structure_idxs : :obj:`numpy.ndarray`, ndim: ``1``
        Indices of source structures we can sample from.
    entity_ids_samples : :obj:`tuple` (:obj:`tuple`)
        A nested tuple containing all possible entity selections for all desired
        components of the sampled structure. Has length ``n``, where ``n`` is
        the number of entities in the sampled structure. Each element is another
        :obj:`tuple` containing all possible source ``entity_ids`` for that
        destination ``entity_id``.

    Yields
    ------
    :obj:`list`
        Structure index and ``entity_id`` selections.

    Notes
    -----
    This could generate selections of multiple of the same entities.
    For example, ``[602, 0, 70, 70]`` selects entity ``70`` twice which is
    nonphysical. For now we will offload this check to ``sample_structures``.
    """
    # Sampling a specified number of structures.
    if isinstance(quantity, int) or str(quantity).isdigit():
        while True:
            R_selection_idx = randrange(len(structure_idxs))
            entity_selection = [None for entity_ids in entity_ids_samples]
            for i in range(len(entity_selection)):
                entity_selection[i] = choice(
                    tuple(
                        filter(
                            lambda x: x not in entity_selection[: i + 1],
                            entity_ids_samples[i],
                        )
                    )
                )
            selection = [R_selection_idx] + entity_selection
            yield selection
    # Sampling all possible structures.
    elif quantity == "all":
        for R_selection_idx in structure_idxs:
            combs = gen_combs(entity_ids_samples, replacement=False)
            for comb in combs:
                if len(comb) != len(set(comb)):
                    continue
                selection = [R_selection_idx] + list(comb)
                yield selection


def sampler_worker(
    selections,
    Z,
    entity_ids_dest,
    R_source,
    E_source,
    G_source,
    entity_ids_source,
    r_prov_specs_source,
    r_prov_id_source,
    criteria,
    periodic_cell,
):
    """Given generated selections will slice and sample structures from source."""
    if not isinstance(selections, np.ndarray):
        selections = np.array(selections)
    assert selections.ndim == 2

    n_Z = len(Z)
    n_selections = selections.shape[0]
    n_r_prov_spec_columns = int(len(selections[0])) + 1

    # Setup worker arrays
    r_prov_specs = np.empty((n_selections, n_r_prov_spec_columns), dtype=np.uint32)
    r_prov_specs_check = np.empty(
        (n_selections, n_r_prov_spec_columns), dtype=np.uint32
    )
    r_prov_spec_selection = np.empty(n_r_prov_spec_columns, dtype=np.uint32)
    R = np.empty((n_selections, n_Z, 3), dtype=np.float64)
    # If these are not None, then they were requested.
    if E_source is not None:
        E = np.empty(n_selections, dtype=np.float64)
    else:
        E = None
    if G_source is not None:
        G = np.empty((n_selections, n_Z, 3), dtype=np.float64)
    else:
        G = None

    # True or False to keep these selections based on uniqueness or criteria.
    keep_idxs = np.empty(n_selections, dtype=np.bool8)

    i_sel = 0
    for i_sel in range(len(selections)):
        selection = selections[i_sel]
        sel_source_idx = selection[0]
        entity_ids_selection = selection[1:]  # Used to slice source

        ###   Handle r_prov_specs   ###
        # Determines the correct prov specification based on the sample and source.
        # We have already determined r_prov_id_source for an original source.
        # Just need to retrieve the r_prov_id_source of a source that is
        # from sampled structures.
        if r_prov_specs_source is not None:
            orig_r_prov_specs = r_prov_specs_source[sel_source_idx]
            source_r_prov_id = orig_r_prov_specs[0]

            # Original r_prov_id and structure index
            r_prov_spec_selection[:2] = orig_r_prov_specs[:2]
            # Original entity_ids
            r_prov_spec_selection[2:] = orig_r_prov_specs[2:][selection[1:]]
        else:
            r_prov_spec_selection = np.array([r_prov_id_source, *selection])

        # Check if selection is already in the destination.
        # This is done by checking if the selection r_prov_spec is already included.
        # If it is, we do not include this selection.
        # We do this by adding the appropriately sorted r_prov_spec_selection to
        # a test array (we only sort the entity_ids).
        r_prov_specs_check[i_sel][:2] = r_prov_spec_selection[:2]
        r_prov_specs_check[i_sel][2:] = np.sort(r_prov_spec_selection[2:], axis=0)

        # Check if sample exists in this batch.
        duplicate_spec = exists_in_array(
            r_prov_specs_check[i_sel], r_prov_specs_check[:i_sel]
        )
        if duplicate_spec:
            keep_idxs[i_sel] = False
            continue

        # Slice r_sel
        R[i_sel] = r_from_entities(
            R_source[sel_source_idx], entity_ids_source, entity_ids_selection
        )

        # Enforce minimum image convention if periodic.
        if periodic_cell is not None:
            r_sel_periodic = periodic_cell.r_mic(R[i_sel])
            if r_sel_periodic is None:
                continue
            else:
                R[i_sel] = r_sel_periodic

        # Checks any structural criteria.
        if criteria is not None:
            accept_r, _ = criteria.accept(Z, R[i_sel], {"entity_ids": entity_ids_dest})
            # If descriptor is not met, will not include sample.
            if not accept_r:
                keep_idxs[i_sel] = False
                continue

        r_prov_specs[i_sel] = r_prov_spec_selection

        if E_source is not None:
            E[i_sel] = E_source[sel_source_idx]
        if G_source is not None:
            G[i_sel] = r_from_entities(
                G_source[sel_source_idx], entity_ids_source, entity_ids_selection
            )

        # Successful sample.
        keep_idxs[i_sel] = True

    R = R[keep_idxs]
    if E is not None:
        E = E[keep_idxs]
    if G is not None:
        G = G[keep_idxs]
    r_prov_specs = r_prov_specs[keep_idxs]

    return R, E, G, r_prov_specs


class Sampler(object):
    """Randomly sample structures."""

    def __init__(
        self,
        source_file,
        source_key,
        dest_file,
        dest_key,
        criteria=None,
        center_structures=False,
        E_key=None,
        G_key=None,
        dry_run=False,
        all_init_size=50000,
        use_ray=False,
        n_workers=2,
        ray_address="auto",
    ):
        """
        Parameters
        ----------
        source_file : :obj:`reptar.File`
            A loaded reptar File to sample structures from.
        source_key : :obj:`str`
            Key to the desired source group. The group must contain
            ``atomic_numbers``, ``geometry``, ``entity_ids``, and ``comp_ids``
            keys.
        dest_file : :obj:`reptar.File`
            A reptar File to add sampled structures to.
        dest_key : :obj:`str`
            Key to the desired destination group. If it does not exist then it
            will be created.
        criteria : ``reptar.descriptor.criteria``, default: ``None``
            Criteria object used to accept or reject a structure based on some
            descriptor.
        center_structures : :obj:`bool`, default: ``False``
            Center each structure by translating the center of mass to the
            origin.
        E_key : :obj:`str`, default: ``None``
            If not ``None``, this is the key to create and attempt to copy
            energies to ``dest_key``. For example, ``energy_ele_mp2.def2tzvp``.
        G_key : :obj:`str`, default: ``None``
            If not ``None``, this is the key to create and attempt to copy
            gradients to ``dest_key``. For example, ``grads_mp2.def2tzvp``.
        dry_run : :obj:`bool`, default: ``False``
            Sample but do not save structures.
        all_init_size : :obj:`int`, default: ``50000``
            Number of structures to initialize sampling arrays with when
            ``quantity`` is ``'all'``.
        use_ray : :obj:`bool`, default: ``False``
            **Not implemented yet.**
            Use ray to parallelize sampling. If ``False``, then sampling is
            done serially in batches. ``True`` is useful if tons of sampling
            on the order of thousands is desired.
        n_workers : :obj:`int`, default: ``2``
            Number of workers to use for sampling if ``use_ray`` is ``True``.
            Each worker will have one core.
        ray_address : :obj:`str`, default: ``'auto'``
            Ray cluster address to connect to.
        """
        self.source_file = source_file
        self.source_key = source_key
        self.dest_file = dest_file
        self.dest_key = dest_key
        self.criteria = criteria
        self.center_structures = center_structures
        self.E_key = E_key
        self.G_key = G_key
        self.dry_run = dry_run
        self.all_init_size = all_init_size

        self.use_ray = use_ray
        self.n_workers = n_workers
        if self.use_ray:
            global ray
            import ray

            if not ray.is_initialized():
                ray.init(address=ray_address)

            # We put criteria in ray object store here.
            self.criteria = ray.put(self.criteria)

        self.worker_chunk_size_for_all = 1000
        """Chunk size used when ``quantity`` is ``'all'``.

        :type: : :obj:`int`
        """

    def _prepare_destination(self):
        """Check for existing sampled structures in destinations and prepare."""
        dest_file = self.dest_file
        dest_key = self.dest_key

        # Will be None if it does not exist.
        self.Z = dest_file.get(f"{dest_key}/atomic_numbers", missing_is_none=True)
        self.R = dest_file.get(f"{dest_key}/geometry", missing_is_none=True)
        # Stores the index of the first newly sampled structure for checks.
        if self.R is None:
            self.n_R_initial = 0
        else:
            self.n_R_initial = self.R.shape[0]

        # Get energies and gradients if available and requested.
        # TODO: perform check if we should copy energies and gradients based on if the
        # sampled structure is the entire source structure.
        if self.E_key is not None:
            self.E = dest_file.get(f"{dest_key}/{self.E_key}", missing_is_none=True)
        else:
            self.E = None
        if self.G_key is not None:
            self.G = dest_file.get(f"{dest_key}/{self.G_key}", missing_is_none=True)
        else:
            self.G = None

        self.r_prov_ids = dest_file.get(f"{dest_key}/r_prov_ids", missing_is_none=True)
        self.r_prov_specs = dest_file.get(
            f"{dest_key}/r_prov_specs", missing_is_none=True
        )

        self.entity_ids = dest_file.get(f"{dest_key}/entity_ids", missing_is_none=True)
        self.comp_ids = dest_file.get(f"{dest_key}/comp_ids", missing_is_none=True)

    def _check_comp_ids(self, comp_labels):
        """Check if desired component IDs to sample is the same as the
        destination if they already exist.
        """
        if self.comp_ids is not None:
            try:
                assert np.array_equal(self.comp_ids, np.array(comp_labels))
            except AssertionError:
                e = f"Component IDs of destination ({comp_ids.tolist()}) do not match comp_labels ({comp_labels})."
                raise AssertionError(e)
        else:
            self.comp_ids = np.array(comp_labels)

    def _prepare_source(self, R_source_idxs):
        """Check for data from source."""
        source_file = self.source_file
        source_key = self.source_key

        self.Z_source = source_file.get(f"{source_key}/atomic_numbers")
        self.R_source = source_file.get(f"{source_key}/geometry")
        self.entity_ids_source = source_file.get(f"{source_key}/entity_ids")
        self.comp_ids_source = source_file.get(f"{source_key}/comp_ids")

        if self.E_key is not None:
            self.E_source = source_file.get(f"{source_key}/{self.E_key}")
        else:
            self.E_source = None
        if self.G_key is not None:
            self.G_source = source_file.get(f"{source_key}/{self.G_key}")
        else:
            self.G_source = None

        if R_source_idxs is None:
            self.R_source_idxs = tuple(range(0, len(self.R_source)))

        # Handle if source already has r_prov_ids
        # If these are None: This source is an original (was not created from sampling).
        self.r_prov_ids_source = source_file.get(
            f"{source_key}/r_prov_ids", missing_is_none=True
        )
        self.r_prov_specs_source = source_file.get(
            f"{source_key}/r_prov_specs", missing_is_none=True
        )

        # Handle periodic structures
        try:
            periodic_cell_vectors = source_file.get(f"{source_key}/periodic_cell")
            periodic_mic_cutoff = source_file.get(f"{source_key}/periodic_mic_cutoff")
            self.periodic_cell = Cell(periodic_cell_vectors, periodic_mic_cutoff)
        except RuntimeError:
            self.periodic_cell = None

    def _check_r_prov_ids(self):
        """ """
        # Prepare source r_prov_ids and r_prov_specs.
        # If source is an original, we just need to determine the next prov_id
        # and prov_specs stays None.
        # If the source was created with sampling, we need to adjust the
        # prov_ids and prov_specs to not overlap with the destination.
        # For example, change 0 from source to 1 if 0 is already taken in destination.

        r_prov_ids_source = self.r_prov_ids_source
        source_file = self.source_file
        source_key = self.source_key

        r_prov_ids = self.r_prov_ids

        # Original source (not from sampled)
        if r_prov_ids_source is None:
            try:
                md5_source = source_file.get(f"{source_key}/md5_structures")
            except RuntimeError:
                source_file.update_md5()
                md5_source = source_file.get(f"{source_key}/md5_structures")

            # Create pseudo r_prov_ids_source.
            r_prov_ids_source = {0: md5_source}

        # If is this not a new destination we cannot always reuse the source information.
        # Need to check for overlap of ids and specifications.
        # Just need to shift the new source IDs up by the maximum destination ID.
        if r_prov_ids is not None:

            # Check and remove any IDs already in the destination.
            present_ids = tuple(key for key in r_prov_ids.keys())
            present_md5s = tuple(value for value in r_prov_ids.values())
            remove_source_ids = []
            for source_id, source_md5 in r_prov_ids_source.items():
                if source_md5 in present_md5s:
                    update_id = present_ids[present_md5s.index(source_md5)]
                    if update_id != source_id:
                        if source_r_prov_specs is not None:
                            update_idx = np.argwhere(
                                source_r_prov_specs[:, 0] == source_id
                            )[:, 0]
                            source_r_prov_specs[update_idx, 0] = update_id
                    remove_source_ids.append(source_id)
            if len(remove_source_ids) > 0:
                for source_id in remove_source_ids:
                    del r_prov_ids_source[source_id]

            # Create any new source IDs that will be added to destination.
            if len(r_prov_ids_source) > 0:
                max_dest_prov_id = max(r_prov_ids.keys())
                next_id = max_dest_prov_id + 1

                for source_id in tuple(r_prov_ids_source.keys()):
                    r_prov_ids_source[next_id] = r_prov_ids_source.pop(source_id)
                    if source_r_prov_specs is not None:
                        update_idx = np.argwhere(
                            source_r_prov_specs[:, 0] == source_id
                        )[:, 0]
                        source_r_prov_specs[update_idx, 0] = next_id
                    next_id += 1

                # Merge the IDs
                self.r_prov_ids = {**r_prov_ids, **r_prov_ids_source}
            else:
                self.r_prov_ids = r_prov_ids
        else:
            self.r_prov_ids = r_prov_ids_source

    def _prepare_saver(self):
        dest_key = self.dest_key

        saver_keys = [f"{dest_key}/geometry", f"{dest_key}/r_prov_specs"]
        if not self.dry_run:
            if self.E_key is not None:
                saver_keys.append(f"{dest_key}/{self.E_key}")
            if self.G_key is not None:
                saver_keys.append(f"{dest_key}/{self.G_key}")
            self.saver = Saver(self.dest_file.fpath, saver_keys)

    def get_avail_entities(self, comp_ids_source, comp_labels, specific_entities=None):
        """Determines available ``entity_ids`` for each ``comp_id`` in
        requested sampling components..

        Parameters
        ----------
        comp_ids_source : :obj:`numpy.ndarray`, ndim: ``1``
            Component IDs of the source.
        comp_labels : :obj:`tuple`, ndim: ``1``
            The ``comp_id`` labels to include in sampled structures.
            ``entity_ids`` with the same ``comp_id`` labels are sampled for
            each :obj:`tuple` element. Thus, the order of each label
            **does** matter.
        specific_entities : :obj:`tuple`, default: ``None``
            Supersede entities determined from
            :meth:`reptar.Sampler.get_avail_entities` and use these instead. This
            needs to specify the ``comp_label`` index to supersede and the
            entities to use instead. For example, if you wanted to always sample
            entity ``43`` or ``78`` for the second ``comp_label`` then this
            would be ``((1, (43, 78)),)``

        Returns
        -------
        :obj:`list`
            All possible entity_ids for each ``comp_id`` provided in
            ``comp_labels``.
        """
        avail_entity_ids = []
        for i in range(len(comp_labels)):
            # Account for manually specified entity_ids
            if specific_entities is not None:
                for j in range(len(specific_entities)):
                    if specific_entities[j][0] == i:
                        avail_entity_ids.append(specific_entities[j][1])
                        continue
            # Find all possible entity_ids to sample from.
            comp_id = comp_labels[i]
            matching_entity_ids = np.where(comp_ids_source == comp_id)[0]
            avail_entity_ids.append(matching_entity_ids)
        return avail_entity_ids

    def _prepare_dest_const_data(self):
        """Take care of destination constant data like ``atomic_numbers`` and
        ``entity_ids``, ``comp_ids``. Also puts the data to the destination.
        """
        # Check that all atomic numbers of each entity is consistent.
        # Check that the z for these new samples match the destination.
        # Also build the entity_ids of these structures.
        Z_sample = []
        entity_ids = []
        entity_id = 0

        Z = self.Z
        entity_ids_source = self.entity_ids_source
        Z_source = self.Z_source

        for dest_entity in self.avail_entity_ids:
            entity_ref = dest_entity[0]  # an example entity_id

            Z_ref = Z_source[np.argwhere(entity_ids_source == entity_ref)[:, 0]]
            entity_ids.extend([entity_id for _ in range(len(Z_ref))])

            Z_sample.extend(Z_ref)

            # Check that all other entities have the same atomic numbers.
            if len(dest_entity) != 1:
                for other_entity in dest_entity[1:]:
                    Z_other = Z_source[
                        np.argwhere(entity_ids_source == other_entity)[:, 0]
                    ]
                    try:
                        assert np.array_equal(Z_ref, Z_other)
                    except AssertionError:
                        raise AssertionError(
                            f"Atomic numbers do not match for entity_id of {other_entity}"
                        )

            entity_id += 1

        Z_sample = np.array(Z_sample)

        # Check that, if the destination already contains atomic_numbers, that
        # they are the same as these structures to be sampled.
        if Z is not None:
            try:
                assert np.array_equal(Z, Z_sample)
                Z = np.array(Z_sample)
            except AssertionError:
                print(f"Destination atomic numbers: {Z}")
                print(f"Sample atomic numbers: {Z_sample}")
                e = "Atomic numbers of samples do not match destination."
                raise AssertionError(e)
        else:
            self.Z = Z_sample
        self.n_Z = len(self.Z)

        self.entity_ids = np.array(entity_ids)  # Destination entity_ids

        if not self.dry_run:
            self.dest_file.put(f"{self.dest_key}/r_prov_ids", self.r_prov_ids)
            self.dest_file.put(f"{self.dest_key}/atomic_numbers", self.Z)
            self.dest_file.put(f"{self.dest_key}/entity_ids", self.entity_ids)
            self.dest_file.put(f"{self.dest_key}/comp_ids", self.comp_ids)

    def _initialize_structure_sampling_arrays(self, comp_labels, quantity):
        """Creates or extends arrays for sampling structures.

        Parameters
        ----------
        comp_labels : :obj:`tuple`, ndim: ``1``
            The ``comp_id`` labels to include in sampled structures. ``entity_ids``
            with the same ``comp_id`` labels are sampled for each :obj:`tuple`
            element. Thus, the order of each label does matter.
        quantity : :obj:`int` or :obj:`str`
            Number of structures to sample from the source. This can be a
            specific number (e.g., ``5739`` or ``'8'``) or ``'all'`` which
            will sample every possible combination.
        """
        # Either creates new ones or concatenates to original destination data.
        # Trial selection to get shapes for R and G
        if quantity == "all":
            # Determine the maximum amount of structures we can sample from a single
            # source structure.
            combs_per_structure = sum(
                1
                for _ in _generate_structure_samples(
                    quantity, [0], self.avail_entity_ids
                )
            )
            size_quantity = int(combs_per_structure * len(self.R_source))

            # Sometimes we want to sample all with a criteria. Criteria should
            # drastically reduce the number of sampled structures. So, we just
            # initialize a large array with ``all_size``.
            if size_quantity > self.all_init_size:
                size_quantity = self.all_init_size
        else:
            size_quantity = int(quantity)
        idx_sel, R = _add_structures_to_R(self.R, size_quantity, self.n_Z)

        if self.E_key is not None:
            E = _add_structures_to_E(self.E, size_quantity)
        else:
            E = None

        if self.G_key is not None:
            G = _add_structures_to_G(self.G, size_quantity, self.n_Z)
        else:
            G = None

        r_prov_specs = _add_structures_to_r_prov_specs(
            self.r_prov_specs, size_quantity, comp_labels
        )

        return idx_sel, R, E, G, r_prov_specs

    def _post_process(self):
        # Put the final data.
        self.dest_file.put(f"{self.dest_key}/r_centered", self.center_structures)
        self.dest_file.put(f"{self.dest_key}/reptar_version", reptar_version)
        self.dest_file.update_md5(self.dest_key)

    def _cleanup(self):
        """Remove attributes after sampling to free memory and avoid using
        old data if new samples are taken.
        """
        # TODO: cleanup
        # Clean up class attributes.
        attrs = [
            "Z",
            "R",
            "n_R_initial",
            "E",
            "G",
            "r_prov_ids",
            "r_prov_specs",
            "entity_ids",
            "comp_ids",
        ]
        for attr in attrs:
            delattr(self, attr)

    def _update_sampling_arrays(
        self,
        quantity,
        i_start,
        R,
        E,
        G,
        r_prov_specs,
        R_selection,
        E_selection,
        G_selection,
        r_prov_spec_selection,
    ):
        do_break = False
        n_sampled = R_selection.shape[0]
        i_stop = i_start + n_sampled

        if i_stop > R.shape[0]:
            if quantity == "all":
                # If we are sampling all, we need to resize the arrays
                # an continue sampling.
                extra_to_add = 5000
                _, R = _add_structures_to_R(R, extra_to_add, self.n_Z)
                r_prov_specs = _add_structures_to_r_prov_specs(
                    r_prov_specs, extra_to_add, comp_labels
                )
                if self.E_key is not None:
                    E = _add_structures_to_E(E, extra_to_add)
                if self.G_key is not None:
                    G = _add_structures_to_G(G, extra_to_add, self.n_Z)
            else:
                # If we are sampling a specific quantity, then we only
                # keep a certain amount.
                i_stop = R.shape[0]
                n_keep = i_stop - i_start
                R_selection = R_selection[:n_keep]
                r_prov_spec_selection = r_prov_spec_selection[:n_keep]
                if self.E_key is not None:
                    E_selection = E_selection[:n_keep]
                if self.G_key is not None:
                    G_selection = G_selection[:n_keep]
                do_break = True

        # Add batch to sample arrays.
        R[i_start:i_stop] = R_selection
        r_prov_specs[i_start:i_stop] = r_prov_spec_selection
        if self.E_key is not None:
            E[i_start:i_stop] = E_selection
        if self.G_key is not None:
            G[i_start:i_stop] = G_selection

        return i_stop, do_break, R, E, G, r_prov_specs

    def sample(self, comp_labels, quantity, R_source_idxs=None, specific_entities=None):
        """

        Parameters
        ----------
        comp_labels : :obj:`tuple`, ndim: ``1``
            The ``comp_id`` labels to include in sampled structures.
            ``entity_ids`` with the same ``comp_id`` labels are sampled for
            each :obj:`tuple` element. Thus, the order of each label
            **does** matter.
        quantity : :obj:`int` or :obj:`str`
            Number of structures to sample from the source. This can be a
            specific number (e.g., ``5739`` or ``'8'``) or ``'all'`` which
            will sample every possible combination.
        R_source_idxs : :obj:`int` or :obj:`list`, default: ``None``
            Which structures (by index) we can sample from. If ``None``, then
            all structures are available.
        specific_entities : :obj:`tuple`, default: ``None``
            Supersede entities determined from ``Sampler.avail_entities`` and
            use these instead. This needs to specify the ``comp_label`` index
            to supersede and the entities to use instead. For example, if you
            wanted to always sample entity ``43`` or ``78`` for the second
            ``comp_label`` then this would be ``((1, (43, 78)),)``
        """
        self._prepare_destination()
        self._prepare_source(R_source_idxs)
        self._check_comp_ids(comp_labels)
        self._check_r_prov_ids()
        self._prepare_saver()

        # Explicitly loading data that we need to check and handle using ray.
        Z = self.Z

        # Determines what the r_prov_id will be for this sampling run.
        # If source_r_prov_specs is None, then this is an original source.
        # Which means there should only be one R_prov_id.
        if self.r_prov_specs_source is None:
            assert len(self.r_prov_ids) == 1
            r_prov_id_source = tuple(self.r_prov_ids.keys())[0]
        else:
            r_prov_id_source = None

        self.avail_entity_ids = self.get_avail_entities(
            self.comp_ids_source, comp_labels, specific_entities
        )
        self._prepare_dest_const_data()

        # We initialize all arrays that will be used to store sampled structures.
        i_start, R, E, G, r_prov_specs = self._initialize_structure_sampling_arrays(
            comp_labels, quantity
        )

        selection_gen = _generate_structure_samples(
            quantity, self.R_source_idxs, self.avail_entity_ids
        )

        # TODO: Print status?

        if quantity != "all":
            chunk_size = quantity
        else:
            chunk_size = self.worker_chunk_size_for_all
        
        stop_sampling = False

        if not self.use_ray:

            # Serial operation
            for selections in chunk_iterable(selection_gen, chunk_size):
                # Sample selections with the worker.
                # All selections are already checked and can be added.
                (
                    R_selection,
                    E_selection,
                    G_selection,
                    r_prov_spec_selection,
                ) = sampler_worker(
                    selections,
                    self.Z,
                    self.entity_ids,
                    self.R_source,
                    self.E_source,
                    self.G_source,
                    self.entity_ids_source,
                    self.r_prov_specs_source,
                    r_prov_id_source,
                    self.criteria,
                    self.periodic_cell,
                )

                (
                    i_start,
                    stop_sampling,
                    R,
                    E,
                    G,
                    r_prov_specs,
                ) = self._update_sampling_arrays(
                    quantity,
                    i_start,
                    R,
                    E,
                    G,
                    r_prov_specs,
                    R_selection,
                    E_selection,
                    G_selection,
                    r_prov_spec_selection,
                )

                if self.center_structures:
                    R[:i_start] = get_center_structures(self.Z, R[:i_start])

                if not self.dry_run:
                    data_to_save = [R[:i_start], r_prov_specs[:i_start]]
                    if self.E_key is not None:
                        data_to_save.append(E[:i_start])
                    if self.G_key is not None:
                        data_to_save.append(G[:i_start])
                    self.saver.save(*data_to_save)

                # Break sampling when quantity is reached.
                # 'all' sampling will break with the loop.
                if stop_sampling:
                    break
        else:
            # Parallel operation with ray
            chunk_size = int(chunk_size / self.n_workers)

            sampler_worker_ray = ray.remote(sampler_worker)

            # Put all data (other than criteria).
            self.entity_ids = ray.put(self.entity_ids)
            self.R_source = ray.put(self.R_source)
            self.E_source = ray.put(self.E_source)
            self.G_source = ray.put(self.G_source)
            self.entity_ids_source = ray.put(self.entity_ids_source)
            self.r_prov_specs_source = ray.put(self.r_prov_specs_source)
            self.periodic_cell = ray.put(self.periodic_cell)

            # Initialize ray workers
            workers = []
            selection_chunker = chunk_iterable(selection_gen, chunk_size)

            def add_worker(workers, chunker):
                try:
                    selection = next(selection_chunker)
                    workers.append(
                        sampler_worker_ray.options(num_cpus=1).remote(
                            selection,
                            self.Z,
                            self.entity_ids,
                            self.R_source,
                            self.E_source,
                            self.G_source,
                            self.entity_ids_source,
                            self.r_prov_specs_source,
                            r_prov_id_source,
                            self.criteria,
                            self.periodic_cell,
                        )
                    )
                except StopIteration:
                    pass

            for _ in range(self.n_workers):
                add_worker(workers, selection_chunker)

            # Start calculations
            while not stop_sampling:
                done_id, workers = ray.wait(workers)

                R_selection, E_selection, G_selection, r_prov_spec_selection = ray.get(
                    done_id
                )[0]

                (
                    i_start,
                    stop_sampling,
                    R,
                    E,
                    G,
                    r_prov_specs,
                ) = self._update_sampling_arrays(
                    quantity,
                    i_start,
                    R,
                    E,
                    G,
                    r_prov_specs,
                    R_selection,
                    E_selection,
                    G_selection,
                    r_prov_spec_selection,
                )

                if self.center_structures:
                    R[:i_start] = get_center_structures(self.Z, R[:i_start])

                if not self.dry_run:
                    data_to_save = [R[:i_start], r_prov_specs[:i_start]]
                    if self.E_key is not None:
                        data_to_save.append(E[:i_start])
                    if self.G_key is not None:
                        data_to_save.append(G[:i_start])
                    self.saver.save(*data_to_save)

                add_worker(workers, selection_chunker)

                # Stop sampling if we have exhausted our chunker and we have
                # no workers left.
                if len(workers) == 0:
                    stop_sampling = True

        if not self.dry_run:
            self._post_process()
        self._cleanup()
