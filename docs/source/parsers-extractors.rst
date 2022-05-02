======================
Parsers and Extractors
======================

Parsers and extractors work together to retrieve data from supported files and store it in :attr:`~reptar.parsers.parser.parser.parsed_info`.

Parsers
=======

A parser a package-specific driver for iterating through each line of output files and storing data.
Each package has its own parser that inherits the base :class:`~reptar.parsers.parser.parser` class which provides common attributes and methods.
For example, we have an xTB parser class:

.. autoclass:: reptar.parsers.parserXTB
    :noindex:

``parse``
---------
After initializing a parser object with the relevant paths, you can then call the :meth:`~reptar.parsers.parser.parser.parse` method.
Each package parser class must define its own :meth:`~reptar.parsers.parser.parser.parse` method.

.. automethod:: reptar.parsers.parser.parse
    :noindex:

``parsed_info``
---------------

Used to store any and all information that can be added to a :ref:`data <data>` object.

.. autoattribute:: reptar.parsers.parser.parsed_info
    :noindex:


Extractors
==========

As previously mentioned, the parser drives the extraction of information from files.
Each extractor is a :class:`~reptar.extractors.extractor.extractor` child class that must have two things.

Triggers
--------

Reptar has to know when to start extracting information while reading the output files.
This is accomplished with a tuple of triggers.

.. autoattribute:: reptar.extractors.extractor.triggers
    :noindex:

Each element is a tuple of length two with a lambda function and corresponding extractor method name as a string.
The lambda function returns ``True`` if a line contains specific words or formatting that signals reptar to parse certain information.
Reptar then calls the respective extractor method.
For example, our ``triggers`` attribute could be something like this:

.. code-block:: python

    triggers = (
        (lambda line: True if ('* xtb version' in line) else False, 'xtb_version'),
        (lambda line: True if (':                      SETUP                      :' == line.strip()) else False, 'gfn_setup'),
        (lambda line: True if ('::                     SUMMARY                     ::' in line.strip()) else False, 'summary_energies'),
        (lambda line: True if ('normal termination of xtb' == line.strip()) else False, 'success'),
    )

When reptar is reading the output file, it evaluates the lambda function and if ``True`` it calls the extractor method.

.. code-block:: python

    for i in range(len(extractor.triggers)):
        if extractor.triggers[i][0](line):
            getattr(extractor, extractor.triggers[i][1])(f, line)

Extractor methods
-----------------

Each trigger has an extractor method that is activated when its lambda function is ``True``.
It will parse information from the line and possibly looping through additional lines.
The extractor's ``parsed_info`` attribute is updated, but nothing is returned.

For example, xTB will print energy summaries like the one below.

.. code-block:: text

    :::::::::::::::::::::::::::::::::::::::::::::::::::::
    ::                     SUMMARY                     ::
    :::::::::::::::::::::::::::::::::::::::::::::::::::::
    :: total energy            -993.709652431952 Eh    ::
    :: gradient norm              0.028600907669 Eh/a0 ::
    :: HOMO-LUMO gap              9.439423838525 eV    ::
    ::.................................................::
    :: HOMO orbital eigv.       -10.046202443852 eV    ::
    :: LUMO orbital eigv.        -0.606778605327 eV    ::
    ::.................................................::
    :: SCC energy             -1002.538277162313 Eh    ::
    :: -> isotropic ES            2.849940828984 Eh    ::
    :: -> anisotropic ES         -0.333608012940 Eh    ::
    :: -> anisotropic XC          0.044384491983 Eh    ::
    :: -> dispersion             -0.696949048082 Eh    ::
    :: repulsion energy           8.590842870552 Eh    ::
    :: add. restraining           0.148264217251 Eh    ::
    :: total charge              -0.000000000001 e     ::
    ::.................................................::
    :: atomisation energy       153.792759436621 Eh    ::
    :::::::::::::::::::::::::::::::::::::::::::::::::::::

Once reptar encounters the ``SUMMARY`` line, it activates the :meth:`~reptar.extractors.xtb_extractor.extractorXTB.summary_energies` method.
Reptar then extracts information such as the total, SCC, and nuclear repulsion energy from this block of text.

Custom extractors
-----------------

Parsers will unfortunately not always extract all desired information.
The best option would be to contribute to reptar and add trigger(s), extractors, and tests.
However, this is not always conducive to a fast-paced research environment.
We would love for you to add functionality to reptar, but sometimes you just need to get your information and move on.

Reptar was designed to be easily accommodate custom extractors for a supported parser.
You can write your own :class:`~reptar.extractors.extractor.extractor` child class that implements its own triggers and extractor methods.

