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

from random import randrange, sample, choice
import itertools
import numpy as np
from .utils import center_structures as get_center_structures
from .utils import get_md5

def _get_R_prov_id(
    dest_R_prov_ids, source_md5
):
    """Determines the structure_prov_id for the sampling source.

    It is possible that 

    Parameters
    ----------
    dest_R_prov_ids : :obj:`dict`
        The r_prov_ids
    source_md5 : obj:`str`
        
    
    Returns
    -------
    :obj:`int`
        Numerical ID for the Rset.
    """
    R_prov_id = None
    if dest_R_prov_ids == {}:
        # No previous sampling, so this automatically becomes zero.
        R_prov_id = 0
    else:
        # Check if md5 is already stored.
        for key,value in dest_R_prov_ids.items():
            if value == source_md5:
                R_prov_id = key
                break
        # Generate a new R_prov_id because it was not already stored.
        if R_prov_id is None:
            R_prov_id = max(dest_R_prov_ids.keys()) + 1

    return R_prov_id

def _initialize_structure_sampling_arrays(
    Z, R, E, F, R_prov_specs, quantity, comp_labels, R_source, entity_ids_source,
    entity_ids_samples, copy_EF
):
    """Creates or extends arrays for sampling structures.

    Parameters
    ----------
        
    
    Returns
    -------

    """
    # Either creates new ones or concatenates to
    # original destination data.
    # Trial selection to get shapes for R and F
    if quantity == 'all':
        # Determine the maximum amount of structures we can sample from a single
        # source structure.
        combs_per_structure = sum(
            1 for _ in itertools.product(*entity_ids_samples)
        )
        size_quantity = int(combs_per_structure*len(R_source))
        sampling_all = True
    else:
        size_quantity = int(quantity)
        sampling_all = False
    
    # Determine number of atoms.
    entity_id_for_size = (i[0] for i in entity_ids_samples)
    number_of_atoms = 0
    for entity_id in entity_id_for_size:
        number_of_atoms += len(np.where(entity_ids_source==entity_id)[0])
    
    if Z is None:
        Z = np.zeros((number_of_atoms,))
    else:
        assert len(Z) == number_of_atoms
    
    sampled_R_empty = np.empty((size_quantity, number_of_atoms, 3))
    sampled_R_empty[:] = np.nan
    if R is None:
        idx_selection = 0  # Starting index of the newley sampled structures.
        R = sampled_R_empty
    else:
        idx_selection = len(R)
        R = np.concatenate((R, sampled_R_empty), axis=0)
    
    if copy_EF:
        sampled_E_empty = np.empty(sampled_R_empty.shape[0])
        sampled_E_empty[:] = np.nan
        if E is None:
            E = sampled_E_empty
        else:
            E = np.concatenate((E, sampled_E_empty))
        
        sampled_F_empty = np.empty((size_quantity, number_of_atoms, 3))
        sampled_F_empty[:] = np.nan
        if F is None:
            F = sampled_F_empty
        else:
            F = np.concatenate((F, sampled_F_empty), axis=0)
    else:
        E = None
        F = None
    
    sampled_R_prov_specs_empty = np.empty((size_quantity, 2+len(comp_labels)))
    sampled_R_prov_specs_empty[:] = np.nan
    if R_prov_specs is None:
        R_prov_specs = sampled_R_prov_specs_empty
    else:
        R_prov_specs = np.concatenate(
            (R_prov_specs, sampled_R_prov_specs_empty), axis=0
        )
    
    return (Z, R, E, F, R_prov_specs, sampling_all, idx_selection)



def _generate_structure_samples(
    quantity, structure_idxs, entity_ids_samples
):
    """Randomly generates structure and entity selections.

    The prov ID for the structure is not included here and must be
    inserted later.

    Parameters
    ----------
    quantity : :obj:`int` or :obj:`float`
        Number of structures to sample from the structure set. For example,
        ``'100'``, ``452``, or even ``'all'``.
    structure_idxs : :obj:`numpy.ndarray`
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
    nonphysical. For now we will offload this check to 
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
                            lambda x: x not in entity_selection[:i+1],
                            entity_ids_samples[i]
                        )
                    )
                )
            selection = [R_selection_idx] + entity_selection
            yield selection
    # Sampling all possible structures.
    elif quantity == 'all':
        for R_selection_idx in structure_idxs:
            #all_combinations = unique_prod(*entity_ids_samples)
            for comb in itertools.product(*entity_ids_samples):
                if len(comb) != len(set(comb)):
                    continue
                selection = [R_selection_idx] + list(comb)
                yield selection

def sample_structures(
    source_data, source_key, quantity, comp_labels, R_prov_ids, source_R_prov_specs,
    Z=None, R=None, R_prov_specs=None, structure_idxs=None,
    criteria=None, z_slice=None, cutoff=None, center_structures=False,
    sampling_updates=False, copy_EF=False, E=None, F=None,
    energy_label='energy_ele', force_label='forces'
):
    """Randomly samples structures from a source.

    If previous data are already provided (i.e., ``Z``, ``R``, and
    ``R_prov_specs`` are not ``None``) then the sampled structures are appended
    to these arrays

    Parameters
    ----------
    source_data : ``exdir.Group``
        A loaded ``exdir.Group`` to sample structures from. Must contain
        ``z``, ``r``, ``entity_ids``, and
        ``comp_ids``.
    source_key : :obj:`str`
        Key to the desired source group.
    quantity : :obj:`str`
        Number of structures to sample from the structure set. For example,
        ``'100'``, ``'452'``, or even ``'all'``.
    comp_labels : :obj:`tuple`
        The ``component_id`` labels to include in sample structures. The order
        of each label does matter.
    R_prov_ids : :obj:`dict`` {:obj:`int`: :obj:`str`}
        Structure prov IDs including any originally defined IDs and new
        ones being added from the source.
    source_R_prov_specs : :obj:`numpy.ndarray` or ``None``
        Structure prov specifications of the source if it was created by
        sampling from other sources. If the ``source`` is original, this should
        be ``None``.
    Z : :obj:`numpy.ndarray`, optional
        Atomic numbers of destination if they already exist. Defaults 
        to ``None``.
    R : :obj:`numpy.ndarray`, optional
        Cartesian atomic coordinates of destination if they already exist.
        Defaults to ``None``.
    R_prov_specs : :obj:`numpy.ndarray`, optional
        Structure prov specifications of destination. Defaults to
        ``None``.
    structure_idxs : :obj:`tuple`, optional
        Source indices of structures to sample from. Defaults to sampling from
        all structures.
    criteria : ``func``, optional
        Structure criteria during the sampling procedure. Defaults to
        ``None`` if no criteria should be used.
    z_slice : :obj:`numpy.ndarray`, optional
        Indices of the atoms to be used for the cutoff calculation. Defaults
        to ``[]`` is no criteria is selected or if it is not required for
        the selected criteria.
    cutoff : :obj:`list`, optional
        Distance cutoff between the atoms selected by ``z_slice``. Must be
        in the same units (e.g., Angstrom) as ``R``. Defaults to ``None`` if
        no criteria is selected or a cutoff is not desired.
    center_structures : :obj:`bool`, optional
        Move the center of mass of each structure to the origin thereby 
        centering each structure (and losing the actual coordinates of
        the structure). While not required for correct use of mbGDML this
        can be useful for other analysis or data set visualization. Defaults
        to ``False``.
    sampling_updates : :obj:`bool`, optional
        Will print something for every 100 successfully sampled structures.
        Defaults to ``False``.
    copy_EF : :obj:`bool`, optional
        Creates datasets for energies and forces (using ``energy_label`` and 
        ``force_label``) and attempts to copy data from the source if possible.
        If no compatible data is available it will just store :obj:`numpy.NaN`.
        Defaults to ``True``.
    E : :obj:`numpy.ndarray`, optional
        Original energies, if any, to append to. Defaults to ``None``.
    F : :obj:`numpy.ndarray`, optional
        Original atomic forces, if any, to append to. Defaults to ``None``.
    energy_label_source : :obj:`str`, optional
        Specifies source dataset name containing energies to copy.
        Defaults to ``'energy_ele'``. If ``copy_EF`` is ``True``, but
        there are no valid energy datasets then we create an empty array.
    force_label_source : :obj:`str`, optional
        Specifies source dataset name containing atomic forces to copy.
        Defaults to ``'forces'``. If ``copy_EF`` is ``True``, but
        there are no valid energy datasets then we create an empty array.
    
    Returns
    -------
    :obj:`numpy.ndarray`
        Atomic coordinates of sampled structures.
    :obj:`numpy.ndarray`
        Atomic coordinates of sampled structures appended to any original
        structures.
    :obj:`numpy.ndarray` or ``None``
        Energies of sampled structures appended to any original structures.
        ``None`` if ``copy_EF`` is ``False``.
    :obj:`numpy.ndarray` or ``None``
        Atomic forces of sampled structures appended to any original structures.
        ``None`` if ``copy_EF`` is ``False``.
    :obj:`numpy.ndarray`
        ``entity_ids`` of the structures.
    :obj:`numpy.ndarray`
        ``R_prov_specs`` of the structures.
    """
    # Get data from structure source.
    Z_source = source_data.get(f'{source_key}/atomic_numbers')
    R_source = source_data.get(f'{source_key}/geometry')
    entity_ids_source = source_data.get(f'{source_key}/entity_ids')
    comp_ids_source = source_data.get(f'{source_key}/comp_ids')

    # If source_R_prov_specs is None, then this is an original source.
    # Which means there should only be one R_prov_id (we store this for later).
    if source_R_prov_specs is None:
        assert len(R_prov_ids) == 1
        source_R_prov_id = tuple(R_prov_ids.keys())[0]

    if copy_EF:
        try:
            E_source = source_data.get(f'{source_key}/{energy_label_source}')
            F_source = source_data.get(f'{source_key}/{force_label_source}')
        except Exception:
            E_source = None
            F_source = None
    else:
        E_source = None
        F_source = None


    ### Setup sampling   ###

    # TODO: Setup option to have one entity be the same all the time.

    ###   Sampling setup   ###
    # Determine what entity IDs we can sample from for each component.
    comp_ids_source_array = np.array(comp_ids_source)
    entity_ids_samples = tuple(
        (
            tuple(np.argwhere(comp_ids_source_array == comp_id)[:,0]) \
            for comp_id in comp_labels
        )
    )
    
    # Check that all atomic numbers of each entity is consistent.
    # Check that the z for these new samples match the destination.
    # Also build the entity_ids of these structures.
    Z_sample = []
    entity_ids = []
    entity_id = 0
    for dest_entity in entity_ids_samples:
        entity_ref = dest_entity[0]
        Z_ref = Z_source[np.argwhere(entity_ids_source == entity_ref)[:,0]]
        entity_ids.extend([entity_id for _ in range(len(Z_ref))])
        entity_id += 1
        Z_sample.extend(Z_ref)
        if len(dest_entity) != 1:
            for other_entity in dest_entity[1:]:
                Z_other = Z_source[np.argwhere(entity_ids_source == other_entity)[:,0]]
                try:
                    assert np.array_equal(Z_ref, Z_other)
                except AssertionError:
                    e = f'Atomic numbers do not match for entity_id of {other_entity}'
                    raise AssertionError(e)
    entity_ids = np.array(entity_ids)
    
    Z_sample = np.array(Z_sample)
    if Z is not None:
        try:
            assert np.array_equal(Z, Z_sample)
        except AssertionError:
            print(f'Destination atomic numbers: {Z}')
            print(f'Sample atomic numbers: {Z_sample}')
            e = 'Atomic numbers of samples do not match destination.'
            raise AssertionError(e)
    else:
        Z = Z_sample


    ###   Initialize all arrays   ###
    Z, R, E, F, R_prov_specs, sampling_all, idx_selection = \
        _initialize_structure_sampling_arrays(
        Z, R, E, F, R_prov_specs, quantity, comp_labels, R_source,
        entity_ids_source, entity_ids_samples, copy_EF
    )

    ###   Sampling counters   ###
    # TODO: clean up progress printing (look into MDAnalysis rdf)
    num_generated = 0  # Number of structures we generated (includes rejected).
    num_accepted = 0  # Number of structures successfully sampled.
    prev_num_accepted_print = 0  # Workaround to control printing.
    
    ###   Begin sampling structures   ###
    if structure_idxs is None:
        structure_idxs = tuple(range(0, len(R_source)))
    for selection in _generate_structure_samples(
        quantity, structure_idxs, entity_ids_samples
    ):
        num_generated += 1

        ###   Sampling maintenance   ###
        # Ends sampling for number quantities.
        if not sampling_all:
            if num_accepted == int(quantity):
                break

        ###   Sampling updates   ###
        # Prints progress information every 500 successful samples.
        if sampling_updates:
            if not sampling_all:
                if num_accepted%500 == 0 \
                    and num_accepted != prev_num_accepted_print:
                    prev_num_accepted_print = num_accepted
                    print(f'Successfully sampled {num_accepted} structures')
            else:
                if (num_accepted+1)%500 == 0:
                    print(
                        f'Successfully sampled {num_accepted+1} structures'
                    )
        
        ###   Handle r_prov_specs   ###
        # Determines the correct prov specification based on the sample and source.
        # We have already determined source_R_prov_id for an original source.
        # Just need to retrieve the source_R_prov_id of a source that is
        # from sampled structures.
        if source_R_prov_specs is not None:
            source_R_prov_id = source_R_prov_specs[selection[0]][0]
        selection.insert(0, source_R_prov_id)


        # Check if selection is already in the destination.
        # If it is, we do not include this selection and try again.
        if (R_prov_specs[:idx_selection]==selection).all(1).any():
            continue
        
        ###   Checks structure criteria   ###
        # Creates mask for atoms in the selection.
        r_index_source = selection[1]
        atom_idx_mask_gen = (
            entity_ids_source == entity_id for entity_id in selection[2:]
        )
        atom_idx_mask = np.bitwise_or.reduce(
            np.array(list(atom_idx_mask_gen)), axis=0
        )

        entity_ids_selection = entity_ids_source[atom_idx_mask]
        r_selection = R_source[r_index_source][atom_idx_mask]
        if criteria is not None:
            # TODO: Update criteria to test for multiple structure (i.e., 3 dimensions)
            accept_r, _ = criteria(
                Z, r_selection, z_slice, entity_ids_selection, cutoff
            )
            # If criteria is not met, will not include sample.
            if not accept_r:
                continue
        
        
        ###   SUCCESSFUL SAMPLE   ###
        R[idx_selection] = r_selection
        R_prov_specs[idx_selection] = selection
        if copy_EF:
            if E_source is not None:
                E[idx_selection] = E_source[r_index_source]
            if F_source is not None:
                F[idx_selection] = F_source[r_index_source][atom_idx_mask]
        
        num_accepted += 1
        idx_selection += 1
    
    if copy_EF:
        E = E[:idx_selection]
        F = F[:idx_selection]

    # Center structures by moving the center of mass to the origin.
    if center_structures:
        R = get_center_structures(Z, R)
    
    return (
        Z, R[:idx_selection], E, F, entity_ids, R_prov_specs[:idx_selection]
    )

def add_structures_to_group(
    source_data, source_key, dest_data, dest_key, quantity,
    comp_labels, structure_idxs=None, consistent_entities=None, criteria=None,
    z_slice=[], cutoff=[], center_structures=False, sampling_updates=False,
    copy_EF=False, energy_label='energy_ele', force_label='forces', write=True
):
    """Adds randomly sampled structures to a ``exdir.Group``.

    Parameters
    ----------
    source_data : ``reptar.data``
        A reptar data object to sample structures from. Must contain
        ``atomic_numbers``, ``geometry``, ``entity_ids``, and ``comp_ids``.
        If it does not contain a ``md5_data`` one will be created and saved.
    source_key : :obj:`str`
        Key to the desired source group.
    dest_data : ``reptar.data``
        A reptar data object to add sampled structures to.
    dest_key : :obj:`str`
        Key to the desired destination group. If it does not
        exist then it will be created.
    quantity : :obj:`str` or :obj:`int`
        Number of structures to sample from the data. For example,
        ``'100'``, ``452``, or ``'all'``.
    comp_labels : :obj:`tuple`
        The component ID labels to include in sample structures. The order of
        each label does matter.
    structure_idxs : :obj:`tuple`, optional
        Specific structure indices from source to sample from. Defaults to all
        structures.
    consistent_entities : :obj:`list` [:obj:`int`], optional
        Molecule indices that will be in every selection. Not implemented
        yet. Defaults to ``None``.
    criteria : :obj:`mbgdml.criteria`, optional
        Structure criteria during the sampling procedure. Defaults to
        ``None`` if no criteria should be used.
    z_slice : :obj:`numpy.ndarray`, optional
        Indices of the atoms to be used for the cutoff calculation. Defaults
        to ``[]`` is no criteria is selected or if it is not required for
        the selected criteria.
    cutoff : :obj:`list`, optional
        Distance cutoff between the atoms selected by ``z_slice``. Must be
        in the same units (e.g., Angstrom) as ``R``. Defaults to ``[]`` if
        no criteria is selected or a cutoff is not desired.
    center_structures : :obj:`bool`, optional
        Move the center of mass of each structure to the origin thereby 
        centering each structure (and losing the actual coordinates of
        the structure). While not required for correct use of mbGDML this
        can be useful for other analysis or data set visualization. Defaults
        to ``False``.
    sampling_updates : :obj:`bool`, optional
        Will print something for every 100 successfully sampled structures.
        Defaults to ``False``.
    copy_EF : :obj:`bool`, optional
        Creates datasets for energies and forces (using ``energy_label`` and 
        ``force_label``) and attempts to copy data from the source if possible.
        If no compatible data is available it will just store :obj:`numpy.NaN`.
        Defaults to ``True``.
    energy_label : :obj:`str`, optional
        Specifies the name of the destination dataset with the desired energies.
        Defaults to ``'energy_ele'``. This should be the same as the source.
    force_label : :obj:`str`, optional
        Specifies the name of the destination dataset with the desired atomic forces.
        Defaults to ``'forces'``. This should be the same as the source.
    write : :obj:`bool`, optional
        Write newly sampled to ``Group``. Defaults to True
    
    Notes
    -----
    Nothing is saved until all checks are passed to avoid irreversible changes
    to already existing ``Group``.
    """
    # Grabs data from destination if exists.
    try:
        Z = dest_data.get(f'{dest_key}/atomic_numbers')
    except Exception:
        Z = None
    
    try:
        R = dest_data.get(f'{dest_key}/geometry')
    except Exception:
        R = None
    
    # Get energies and forces if available.
    # TODO: perform check if we should copy energies and forces based on if the
    # sampled structure is the entire source structure.
    if copy_EF:
        E = dest_data.get(f'{dest_key}/{energy_label}')
        F = dest_data.get(f'{dest_key}/{force_label}')
    else:
        E = None
        F = None

    try:
        R_prov_ids = dest_data.get(f'{dest_key}/r_prov_ids')
        R_prov_specs = dest_data.get(f'{dest_key}/r_prov_specs')
    except Exception as e:
        R_prov_ids = None
        R_prov_specs = None
    
    try:
        entity_ids = dest_data.get(f'{dest_key}/entity_ids')
        comp_ids = dest_data.get(f'{dest_key}/comp_ids')
    except Exception:
        entity_ids = None
        comp_ids = None
    
    # Handle if source already has r_prov_ids
    try:
        source_R_prov_ids = source_data.get(f'{source_key}/r_prov_ids')
        source_R_prov_specs = source_data.get(f'{source_key}/r_prov_specs')
    except Exception as e:
        # This source is an original (was not created from sampling).
        source_R_prov_ids = None
        source_R_prov_specs = None
    
    ###   Check component IDs   ###
    if comp_ids is not None:
        try:
            assert np.array_equal(comp_ids, np.array(comp_labels))
        except AssertionError:
            e = f'Component IDs of destination ({comp_ids.tolist()}) do not match comp_labels ({comp_labels}).'
            raise AssertionError(e)
    else:
        comp_ids = np.array(comp_labels)


    # Prepare source r_prov_ids and r_prov_specs.
    # If source is an original, we just need to determine the next prov_id
    # and prov_specs stays None.
    # If the source was created with sampling, we need to adjust the
    # prov_ids and prov_specs to not overlap with the destination.
    # For example, change 0 from source to 1 if 0 is already taken in destination.

    # Original source (not from sampled)
    if source_R_prov_ids is None:
        try:
            md5_source = source_data.get(f'{source_key}/md5_data')
        except Exception:
            md5_source = get_md5(source_data.get(source_key), only_dataset=True)
            source_data.add(f'{source_key}/md5_data', md5_source)
        
        # Create pseudo source_R_prov_ids.
        source_R_prov_ids = {0: md5_source}

    # If is this not a new destination we cannot always reuse the source information.
    # Need to check for overlap of ids and specifications.
    # Just need to shift the new source IDs up by the maximum destination ID.
    if R_prov_ids is not None:

        # Check and remove any IDs already in the destination.
        present_ids = tuple(key for key in R_prov_ids.keys())
        present_md5s = tuple(value for value in R_prov_ids.values())
        remove_source_ids = []
        for source_id,source_md5 in source_R_prov_ids.items():
            if source_md5 in present_md5s:
                update_id = present_ids[present_md5s.index(source_md5)]
                if update_id != source_id:
                    if source_R_prov_specs is not None:
                        update_idx = np.argwhere(
                            source_R_prov_specs[:,0] == source_id
                        )[:,0]
                        source_R_prov_specs[update_idx,0] = update_id
                remove_source_ids.append(source_id)
        if len(remove_source_ids) > 0:
            for source_id in remove_source_ids:
                del source_R_prov_ids[source_id]
        
        # Create any new source IDs that will be added to destination.
        if len(source_R_prov_ids) > 0:
            max_dest_prov_id = max(R_prov_ids.keys())
            next_id = max_dest_prov_id + 1

            for source_id in tuple(source_R_prov_ids.keys()):
                source_R_prov_ids[next_id] = source_R_prov_ids.pop(
                    source_id
                )
                if source_R_prov_specs is not None:
                    update_idx = np.argwhere(
                        source_R_prov_specs[:,0] == source_id
                    )[:,0]
                    source_R_prov_specs[update_idx,0] = next_id
                next_id += 1
    
            # Merge the IDs
            new_R_prov_ids = {**R_prov_ids, **source_R_prov_ids}
        else:
            new_R_prov_ids = R_prov_ids
    else:
        new_R_prov_ids = source_R_prov_ids
    
    # Stores the index of the first newly sampled structure for checks.
    if R is None:
        n_R_initial = 0
    else:
        n_R_initial = R.shape[0]
    
    # Begin sampling.
    Z, R, E, F, entity_ids_sampled, R_prov_specs = sample_structures(
        source_data, source_key, quantity, comp_labels, new_R_prov_ids,
        source_R_prov_specs, Z=Z, R=R, R_prov_specs=R_prov_specs,
        structure_idxs=structure_idxs, criteria=criteria, z_slice=z_slice,
        cutoff=cutoff, center_structures=center_structures,
        sampling_updates=sampling_updates, copy_EF=copy_EF, E=E, F=F,
        energy_label=energy_label, force_label=force_label
    )

    # Check that entity_ids are the same.
    if entity_ids is not None:
        try:
            assert np.array_equal(entity_ids, entity_ids_sampled)
        except AssertionError:
            e = f'Destination entity_ids ({entity_ids.tolist()}) do not '\
                f'match sampled entity_ids ({entity_ids_sampled.tolist()}).'
            raise AssertionError(e)
    else:
        entity_ids = entity_ids_sampled
    
    if write:
        dest_data.add(f'{dest_key}/atomic_numbers', Z)
        dest_data.add(f'{dest_key}/geometry', R)
        if copy_EF:
            dest_data.add(f'{dest_key}/{energy_label}', E)
            dest_data.add(f'{dest_key}/{force_label}', F)
        dest_data.add(f'{dest_key}/entity_ids', entity_ids)
        dest_data.add(f'{dest_key}/r_prov_ids', new_R_prov_ids)
        dest_data.add(f'{dest_key}/r_prov_specs', R_prov_specs)
        dest_data.add(f'{dest_key}/r_centered', center_structures)
        dest_data.add(f'{dest_key}/comp_ids', comp_ids)

        dest_data.add(
            f'{dest_key}/md5',
            get_md5(dest_data.get(dest_key))
        )
        dest_data.add(
            f'{dest_key}/md5_data', get_md5(
                dest_data.get(dest_key), only_dataset=True
            )
        )

    return dest_data
