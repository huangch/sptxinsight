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
     - Required for ``cme --cme-batch-correct harmony``.
   * - ``mcp``
     - ``fastmcp>=2.0``
     - Model Context Protocol server.

Building these docs
-------------------

.. code-block:: bash

   pip install sphinx
   cd docs
   make html      # output in docs/_build/html
