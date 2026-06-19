sptxinsight
===========

Cell-typing, spatial-heterogeneity (H-Plot), and cellular-microenvironment
(CME / niche) analysis for **spatial transcriptomics**, built as a lightweight
sibling of `WSInsight <https://github.com/huangch/wsinsight>`_.

Where WSInsight ingests whole-slide images, ``sptxinsight`` ingests AnnData
spatial samples (``.h5ad`` / ``.zarr``) whose coordinates are already in
microns. It reuses WSInsight's H-Plot engine (vendored under
``sptxinsight.insightlib``) but needs none of the heavy perception stack.

- **Repository:** https://github.com/huangch/sptxinsight
- **License:** Apache-2.0
- **Python:** 3.11+
- **Status:** Alpha

.. toctree::
   :maxdepth: 2
   :caption: Contents

   installing
   cme

Quick start
-----------

.. code-block:: bash

   pip install -e .                    # standalone
   # or, inside the shared wsinsight conda env:
   pip install --no-deps -e .

   sptxinsight run \
     -i ./samples \
     -o ./results \
     --base-type tumor --target-type lymphocyte

Each sample must provide ``adata.obsm["spatial"]`` (N×2 micron coordinates) and
a categorical ``adata.obs["cell_type"]``. Gene-expression and CME workflows
additionally use an expression matrix (``X``, ``raw``, or a named layer).

Indices and tables
------------------

* :ref:`genindex`
* :ref:`search`
