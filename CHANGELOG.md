# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- Drivers and workers should only return guaranteed values.
`DriverOpt` and respective workers no longer return gradients.
- Renamed `xtb_engrad` to `xtb_python_engrad` to reflect the Python interface being used.

### Fixed

- Psi4 optimizations without ray would crash.
- Zarr arrays were classified as attributes when putting data.
Now reptar correctly identifies them as iterables and array-like data.
- Zarr data keys cannot start with ``.``, so ``clean_key`` has been modified to enforce no ``.`` in the beginning.

## [0.1.0] - 2023-04-21

### Added

- Logging capabilities with sparse implementation.
- Basic [zarr](https://zarr.readthedocs.io/en/stable/index.html) File handling.
- Geometry validation during sampling using ``validate_and_fill_geometry`` from
QCElemental. Also added `utils.validate_geometry` function that checks this for
multiple structures.
- Tests for xtb energy+gradient calculations.
- Tests for psi4 energy and energy+gradient calculations.
- Pylint section in documentation.

### Changed

- ``reptar-write-xyz`` now supports all file types.
- Handle 2D arrays in XYZ writers.
- References to ``string`` files are renamed to ``xyz``.
- System information (i.e., atomic numbers and coordinates) are moved to
``driver.run()`` instead of during initialization.
- Unified specification of ray workers and number of cores.
``n_workers`` always defaults to ``1`` which ``n_cpus_per_worker`` specifies the number
of CPU cores per worker (also defaults to ``1``).

### Fixed

- `utils.center_structures` sometimes repeated the array incorrectly.
- Explicitly use 3DMol v1.8.0 for documentation.
- ``psi4_energy`` returned wrong energy shape.
- ``File.put()`` did not create missing groups consistently in exdir files.

## [0.0.3] - 2022-11-30

### Added

- Ability to remove data with the ``File`` object.
- Option (``missing_is_none``) to return ``None`` from ``File.get()`` if the key does not exist.
- ``Saver`` with sampling.
- ``with_md5_update`` option for putting data in reptar files.
- ``File.copy()`` method to help copying data between groups.
- ``driverOPT`` and ``psi4_opt`` worker.
- ``reptar-write-xyz`` script for writing xyz files from command line.
- ``driverENERGY`` for just energy calculations.
- ``psi4_energy`` worker.
- Saver class for long calculations.
- Parse final electronic energies from ORCA.
- Support for CREST conformer and rotamer searches.
- Support for Atomic Simulation Environment (ASE) trajectories.
- Determine structure indices with respect to structure provenance specifications.
- Simple structural descriptors.
- Quick documentation for sampler.
- Example script for GDML npz files.
- Writer for standard XYZ files with optional comments.
- Writer for extended XYZ files used in the Gaussian approximation potentials (GAP) code.
- Documentation for writers.
- Method to store arrays in ASE and schnetpack databases.

### Changed

- Pylint and Black formatting; thus, lots of naming changes were made.
- Renamed ``creator`` class to ``Creator`` to avoid namespace conflicts.
- Moved ``reptar.calculators.save.Saver`` to ``reptar.Saver``.
- Condensed and cleaned sampling routines into ``Sampler`` class.
- Use ``qcelemental.periodictable`` for atomic numbers, symbols, and masses instead of manual dictionaries.
- All drivers have a ``use_ray`` parameter that defaults to ``False``.
- Calculators now use keyword arguments.
- Cleaned up ray calculations.
- Convert scf energies from cclib to Hartree by default.
- Moved example calculations to [reptar-data](https://github.com/aalexmmaldonado/reptar-data).
- Manual ``api`` documents with ``sphinx-multiversion``.
- Raise ``RuntimeError`` when key does not exist in exdir.
- PDB writer numbers each atom name.
- Renamed ``File.add()`` to ``File.put()``.
- Return data from exdir as ``np.ndarray`` instead of ``np.memmap``.
- ``creator.group`` is now ``creator.from_calc``.
- Specify array dimensions in docstrings.
- Do not import every submodule in reptar.
That way, optional calculation dependencies are not required when importing reptar.
- Upgrades to README.
- PDB writer is now a function instead of class.
- PDB writer requires arrays instead of reptar files and group keys.
- Rename ``data`` class to ``File``.
This clears up previous ambiguous usage of "data" to refer to both a file and value of a key.
- Rename reptarWriter to textWriter (more specific).
- Require setting the memory for Psi4 worker.

### Removed

- ``element_to_z``, ``z_to_element``, and ``z_to_mass`` in ``reptar.utils``.
- Many-body expansion routines were incorporated into [mbGDML](https://keithgroup.github.io/mbGDML/index.html).
- ``textWriter`` class (really had no purpose).

### Fixed

- Previously, coordinates were sliced based on a boolean mask during sampling.
The specific order of ``entity_ids`` in ``r_prov_specs`` is no longer correct.
For example, if the selection as ``[0, 382, 82, 45, 100]`` then the actual coordinates would be in order of ``45``, ``82``, and then ``100``, which is incorrect.
Any sampled structures before this commit do not have correct ``r_prov_specs`` with respect to the source.
- ``psi4_engrad`` no longer translates and rotates molecules to origin (which would result in incorrect gradients).
- ``_generate_structure_samples`` when ``quantity`` is all.
Would incorrectly generate combinations with replacement.

## [0.0.2] - 2022-05-03

### Added

- Parallel implementation of Psi4 and xTB energy and gradient calculations with ray.
- Zenodo DOI.

### Changed

- Updating README and documentation home page.
- Document pip install of reptar.
- Sampling structures copies gradients instead of forces.
- Write test files to temporary, untracked directory.

## [0.0.1] - 2022-04-30

- Initial release!
