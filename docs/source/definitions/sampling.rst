========
sampling
========

identification
--------------

.. glossary::

    r_centered
        If structures (likely sampled) were translated such that the center of mass is at the origin.

        **Type:** :obj:`bool`

    r_prov_ids
        Species an ID (int) to uniquely identifying labels for each structure if it originated from another reptar file. Labels should always be ``md5_structures``. For example, ``{0: '6038e101da7fc0085978741832ebc7ad', 1: 'eeaf93dec698de3ecb55e9292bd9dfcb'}``.

        **Type:** :obj:`dict`

    r_prov_specs
        Specifies ``structure_prov_id``, the ``structure index`` and original ``entity_id`` of the reptar prov. Take the following array for example.

        .. code-block:: python

            [[0, 483,   0, 12],
             [0,  89, 194,  2],
             [1,   0,   0,  1]]

        This indicates that we should have three structures in our ``Group``.
        Structure one and two come from the ``Group`` that has a ``structure_prov_id`` of ``0`` whereas the last one is from a different ``Group`` (``structure_prov_id`` of ``1``).
        The first structure in this ``Group`` comes from the 484th (index of ``483``) structure in the source ``Group``.
        But this structure does not contain every atom from the original ``Group`` structure; it only contains the atoms from ``0`` and ``12`` ``entity_id``.

        **Type:** :obj:`numpy.ndarray`


