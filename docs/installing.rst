Installing
==========

Standalone
----------

.. code-block:: bash

   pip install -e .

Inside the shared ``wsinsight`` conda environment
-------------------------------------------------

Install **without dependencies** so ``pip`` cannot upgrade the locked
``numpy<2`` / ``zarr<3`` / ``fsspec`` generation that WSInsight depends on
(every runtime dependency is already present there):

.. code-block:: bash

   pip install --no-deps -e .

Optional extras
---------------

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Extra
     - Adds
     - Note
   * - ``zarr``
     - ``zarr<3``
     - Read ``.zarr`` samples in the shared env.
   * - ``spatialdata``
     - ``spatialdata``
     - Needs ``numpy>=2`` / ``zarr>=3`` — dedicated env only.
   * - ``scanpy``
     - ``scanpy``
     - Same ``numpy>=2`` constraint.
   * - ``harmony``
     - ``harmonypy``
     - Required for ``niche --batch-correct harmony``.
   * - ``mcp``
     - ``fastmcp>=2.0``
     - Model Context Protocol server.
   * - ``kurtorank``
     - ``kurtorank>=3.1.0``
     - Enables KurtoRank-backed annotation commands (for example ``annotate``).

Install with KurtoRank-powered annotation enabled:

.. code-block:: bash

   pip install -e '.[kurtorank]'

Annotation input behavior
-------------------------

``sptxinsight annotate`` accepts either:

* an exact Xenium directory containing files such as
  ``cell_by_feature_matrix.h5`` and ``cells.csv.gz``
* or AnnData inputs (``.h5ad``, ``.zarr``, or an ``sptx-list://`` manifest)

It does **not** auto-append ``/outs``. If your Xenium export lives in
``sample/outs/``, pass ``-i sample/outs`` explicitly.

When ``--markers-csv`` is omitted, the annotation path uses Kurtorank's bundled
default panel, ``markers-v6.csv``.

Building these docs
-------------------

.. code-block:: bash

   pip install sphinx
   cd docs
   make html      # output in docs/_build/html
