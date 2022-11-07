============
Contributing
============

Reptar was designed to be easily expanded to support other research involving:

- parsing of other packages;
- export/writing to other formats;
- or helpful routines.

We welcome all contributions, and they are greatly appreciated!
Every little bit helps, and credit will always be given.



Types of contributions
======================

Report bugs
-----------

Report bugs at https://github.com/aalexmmaldonado/reptar/issues.

If you are reporting a bug, please include:

- Your operating system name and version.
- Any details about your local setup that might be helpful in troubleshooting.
- Detailed steps to reproduce the bug.

Fix bugs
--------

Look through the `GitHub issues <https://github.com/aalexmmaldonado/reptar/issues>`__ for bugs.
Anything tagged with ``bug`` and ``help wanted`` is open to whoever wants to implement it.

Implement features
------------------

Look through the `GitHub issues <https://github.com/aalexmmaldonado/reptar/issues>`__ for features.
Anything tagged with ``enhancement`` and ``help wanted`` is open to whoever wants to implement it.

Write documentation
-------------------

Reptar could always use more documentation, whether as part of the official reptar docs, in docstrings, or even on the web in blog posts, articles, and such.

Propose a new feature
---------------------

The best way to propose a new feature is by starting a discussion at https://github.com/aalexmmaldonado/reptar/discussions.

- Create a discussion in the |:bulb:| Ideas category.
- Explain in detail how it would work.
- Keep the scope as narrow as possible, to make it easier to implement.
- Remember that this is a volunteer-driven project, and that contributions are welcome :)



Discussions
===========

If you have any questions, comments, concerns, or criticisms please start a `discussion <https://github.com/aalexmmaldonado/reptar/discussions>`__ so we can improve reptar!



Get Started!
============

Ready to contribute?
Here's how to set up ``reptar`` for local development.

1. Fork the `reptar repo on GitHub <https://github.com/aalexmmaldonado/reptar>`__.
2. Clone your fork locally.

.. code-block:: bash

    $ git clone git@github.com:your_name_here/reptar.git

3. Install your local copy.

.. code-block:: bash

    $ cd reptar/
    $ pip install .

4. Create a branch for local development.

.. code-block:: bash

    $ git checkout -b name-of-your-branch

Now you can make your changes locally.

5. When you're done making changes, check that your changes pass the tests.

Because parsing is a large part of reptar, we have stored the large output files in a separate `GitHub repo <https://github.com/aalexmmaldonado/reptar-data>`__.
The ``tests/`` is setup to search for files in the ``examples/reptar-data/`` directory.
So, to properly set this up you need to clone this repo into your local repo.

.. code-block:: bash

    $ cd examples
    $ clone https://github.com/aalexmmaldonado/reptar-data

Files in ``examples/reptar-data/``, and changes made to it, are not tracked by reptar.
Then, while in the repo root directory, you can run all of the tests with the ``pytest`` command (after running ``pip install .``).

.. code-block:: bash

    $ pytest
    =========================== test session starts ============================
    platform linux -- Python 3.8.10, pytest-7.1.2, pluggy-1.0.0
    rootdir: /home/alex/repos/reptar, configfile: pytest.ini
    plugins: anyio-3.6.1, cov-3.0.0, order-1.0.1
    collected 19 items                                                         

    tests/test_creator_ase.py .                                          [  5%]
    tests/test_creator_crest.py ..                                       [ 15%]
    tests/test_creator_orca.py ...                                       [ 31%]
    tests/test_creator_xtb.py ....                                       [ 52%]
    tests/test_descriptors.py ..                                         [ 63%]
    tests/test_sampling.py ..                                            [ 73%]
    tests/test_writer_ase_db.py .                                        [ 78%]
    tests/test_writer_pdb.py .                                           [ 84%]
    tests/test_writer_schnetpack_db.py .                                 [ 89%]
    tests/test_writer_xyz.py .                                           [ 94%]
    tests/test_writer_xyz_gap.py .                                       [100%]

    =========================== 19 passed in 20.02s ============================

.. hint::

    We run all of the ``creator`` tests first to build the required files for other tests.
    If one of these fail, it will likely cause others to fail as well.
    Always debug these first.

.. attention::

    If you are implementing new parsers or calculation types you need to include output files for your tests.
    Locally, you can store the files in your cloned ``reptar-data`` directory and run tests that way.
    Once you are ready for merge your changes, you need to add new data to ``reptar-data`` by forking and creating a `pull request <https://github.com/aalexmmaldonado/reptar-data>`__.
    If you need any help doing this, please search the `discussions <https://github.com/aalexmmaldonado/reptar/discussions>`__ or start a new one. 

6. Write any additional documentation in ``docs/source/``.
You can easily build and view the documentation locally by running the ``docs/branch-build-docs.sh`` script then opening ``docs/html/index.html`` in your favorite browser.

.. code-block:: bash

    $ ./docs/branch-build-docs.sh 
    Running Sphinx v5.1.1
    making output directory... done
    loading intersphinx inventory from https://urllib3.readthedocs.io/en/latest/objects.inv...
    loading intersphinx inventory from https://docs.python.org/3/objects.inv...
    loading intersphinx inventory from https://numpy.org/doc/stable/objects.inv...
    building [mo]: targets for 0 po files that are out of date
    building [html]: targets for 67 source files that are out of date
    updating environment: [new config] 67 added, 0 changed, 0 removed
    reading sources... [100%] writers                                                                                         
    /home/alex/repos/reptar/docs/docstring of reptar.parsers.parser.parsed_info:4: WARNING: duplicate term description of system_info, other instance in api/parsers-api/parser
    /home/alex/repos/reptar/docs/docstring of reptar.parsers.parser.parsed_info:9: WARNING: duplicate term description of runtime_info, other instance in api/parsers-api/parser
    /home/alex/repos/reptar/docs/docstring of reptar.parsers.parser.parsed_info:13: WARNING: duplicate term description of outputs, other instance in api/parsers-api/parser
    looking for now-outdated files... none found
    pickling environment... done
    checking consistency... done
    preparing documents... done
    writing output... [100%] writers                                                                                          
    generating indices... genindex done
    highlighting module code... [100%] reptar.writers.xyz_gap                                                                 
    writing additional pages... search done
    copying images... [100%] files/30h2o-md/30h2o.2h2o-com.sum-distribution-13457.png                                         
    copying downloadable files... [100%] files/30h2o-md/30h2o-gfn2-md.exdir.zip                                               
    copying static files... done
    copying extra files... done
    dumping search index in English (code: en)... done
    dumping object inventory... done
    build succeeded, 3 warnings.

    The HTML pages are in html.

7. Add a description of the changes in the ``CHANGELOG.md``.
Please follow the general format specified `here <https://keepachangelog.com/en/1.0.0/>`__.

8. Commit your changes and push your branch to GitHub.

.. code-block:: bash

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-branch

9. Submit a pull request through the GitHub website.



Pull Request Guidelines
=======================

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated.
   Put your new functionality into a function with a docstring, and add the feature to the list in ``CHANGELOG.md``.

.. tip::

    You can open a draft pull request first to check that GitHub actions pass for all supported Python versions.
