Cellular microenvironments (CME / niches)
=========================================

A **cellular microenvironment** (CME, or *niche*) is a recurring local cell
mixture — e.g. a tumor core, an immune-infiltrated rim, or a stromal band.
``sptxinsight cme`` discovers them unsupervised: it builds per-sample Delaunay
cell graphs, gathers k-hop composition features, trains a global DGI encoder,
clusters the embeddings, and writes a one-hot ``cme_<n>`` label per cell. Run it
**after** ``ingest`` / ``run``.

.. code-block:: bash

   # Cell-type niches across all ingested samples:
   sptxinsight cme -o ./results --cme-clusters 8 --cme-k-hops 3

   # Gene-expression niches (k-hop mean expression):
   sptxinsight cme -o ./results --cme-mode expression --cme-batch-correct center

   # Feed every gene to the encoder instead of PCA-reduced expression:
   sptxinsight cme -o ./results --cme-mode expression --disable-pca

Modes
-----

``--cme-mode`` selects what drives the niches and namespaces the outputs so the
families coexist on the same cells:

.. list-table::
   :header-rows: 1
   :widths: 18 32 25 25

   * - ``--cme-mode``
     - features
     - output folder
     - one-hot columns
   * - ``celltype`` (default)
     - k-hop cell-type composition
     - ``cme-outputs-csv/``
     - ``cme_<n>``
   * - ``expression``
     - k-hop mean gene expression (``expr_``)
     - ``cme-gex-outputs-csv/``
     - ``gexcme_<n>``
   * - ``both``
     - composition + expression (fused)
     - ``cme-hybrid-outputs-csv/``
     - ``hcme_<n>``

Run the command once per mode to get **parallel** cell-type and gene niches on
the same cells; ``celltype`` stays byte-identical to earlier releases.

Shared PCA of expression features (default on)
----------------------------------------------

For ``expression`` / ``both`` modes the per-cell gene panel is **reduced to a
shared set of principal components before the k-hop aggregation** (default
``--cme-pca-components 50``). The basis is fit **once on the pooled cohort** with
an ``IncrementalPCA`` and applied identically to every sample, which:

- denoises the sparse gene panel,
- shrinks the encoder input (``50 · hops`` instead of ``n_genes · hops``), and
- keeps niches comparable across samples.

PCA only affects the **encoder input** — the interpretable ``expr_`` columns are
kept for ``cme-profile`` markers. The k-hop feature / embedding checkpoints are
tagged with the PCA setting (e.g. ``slide-graphs-gex-pca50.joblib``) so toggling
PCA never reuses a stale cache.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Option
     - Effect
   * - ``--cme-pca-components N``
     - Number of shared components (min 2; default 50). The effective dimension
       is ``min(N, n_genes)``.
   * - ``--disable-pca``
     - Feed all genes to the encoder instead of the PCA-reduced features.

Other key options
------------------

``--cme-clusters`` (KMeans k; omit for an automatic Leiden sweep),
``--cme-k-hops``, ``--cme-max-edge-len-um``, ``--cme-soft`` (probability instead
of argmax composition), ``--cme-batch-correct`` (``none`` / ``center`` /
``harmony`` cross-sample correction of the embeddings — use a **technical** unit
such as sample/run as the batch, never a biological condition), and
``--cme-regions`` (merge cells into annotation-level regions).
``--cme-expression`` is a deprecated alias for ``--cme-mode both``.

Naming niches
-------------

.. code-block:: bash

   sptxinsight cme-profile -o ./results --top-types 5 --top-genes 10
   sptxinsight cme-profile -o ./results --cme-mode expression

This writes ``cme-profile-composition.csv`` (mean cell-type fractions per CME)
and, when ``expr_`` columns are present, ``cme-profile-markers.csv``. When both
cell-type and gene-expression niches exist, ``cme-profile`` also reports their
**agreement** (normalized mutual information plus a cross-tab) and writes
``cme-agreement.csv``.

Niches as an H-Plot axis
------------------------

.. code-block:: bash

   SPTXINSIGHT_EXPERIMENTAL=1 sptxinsight hplot -o ./results \
     --base-type tumor --base-by celltype \
     --target-type 7 --target-by cme        # fraction of cells in cme_7 per layer

``--base-by`` / ``--target-by`` accept ``celltype`` (default), ``gene``, or a
CME niche family — ``cme``, ``cmegex``, or ``cmehybrid``. When both axes are CME
families they must be the same family.

H-Plot output follows the WSInsight distance contract: ``layer`` is the signed
graph-hop index from the base-region border, and ``distance_um`` is the
cumulative physical distance in microns derived from Delaunay edge lengths.
