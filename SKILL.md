---
name: sptxinsight
description: Install and operate sptxinsight for spatial-transcriptomics cell-typing, H-Plot spatial heterogeneity, and CME / niche discovery on AnnData samples
---

# sptxinsight — Agentic AI Skill File

> **Purpose**: Enable an agentic AI (Claude, OpenClaw, Hermes, or any
> tool-using LLM agent) to autonomously install and operate `sptxinsight` for
> spatial-transcriptomics cell-typing, H-Plot spatial-heterogeneity analysis,
> and cellular-microenvironment (CME / niche) discovery on AnnData samples.

---

## 1. What Is sptxinsight?

`sptxinsight` is a lightweight sibling of
[WSInsight](https://github.com/huangch/wsinsight) for **spatial
transcriptomics**. Where WSInsight ingests whole-slide images, `sptxinsight`
ingests AnnData spatial samples (`.h5ad` / `.zarr`) whose coordinates are
already in microns. It reuses WSInsight's H-Plot engine (vendored under
`sptxinsight.insightlib`) but needs **none** of the heavy perception stack
(no torch / tensorflow / openslide for the core; CME niche discovery adds an
optional PyTorch-Geometric encoder).

- **Repository**: <https://github.com/huangch/sptxinsight>
- **License**: Apache-2.0
- **Python**: ≥ 3.11
- **Status**: Alpha
- **Entry point**: `sptxinsight` (installed via `pip install -e .`)

Each sample must provide `adata.obsm["spatial"]` (N×2 micron coordinates) and a
categorical `adata.obs["cell_type"]`. Gene-expression workflows additionally use
an expression matrix (`X`, `raw`, or a named layer).

---

## 2. Install

Standalone:

```bash
pip install -e .
```

Inside the shared `wsinsight` conda environment, install without dependencies so
`pip` cannot upgrade the locked `numpy<2` / `zarr<3` / `fsspec` generation that
WSInsight depends on:

```bash
pip install --no-deps -e .
```

Optional extras: `zarr` (read `.zarr` in the shared env), `spatialdata` /
`scanpy` (dedicated `numpy>=2` env), `mcp` (`fastmcp>=2.0` MCP server),
`harmony` (`harmonypy`, for `--cme-batch-correct harmony`).

---

## 3. CLI Overview

```bash
sptxinsight --help
```

Global options precede the subcommand: `--backend {anndata,zarr,spatialdata}`
selects the sample loader; `--log-level` sets verbosity.

| Command | Purpose |
|---|---|
| `run` | Ingest → adapt → H-Plot, end to end. |
| `ingest` | Read samples and write the per-sample H-Plot CSV contract. |
| `annotate` | Verify samples are cell-typed and report per-type counts. |
| `export` | Print the path to the aggregated H-Plot table. |
| `describe` | Emit a JSON schema of every subcommand (for tooling / MCP). |
| `cme` | Discover cellular microenvironments (niches) across ingested samples. |
| `cme-profile` | Summarise each CME's cell composition and marker genes to name niches. |
| `hplot`, `hplot-finalize` | Experimental; hidden unless `SPTXINSIGHT_EXPERIMENTAL=1`. |
| `agg` | Experimental; detect cell-type aggregates (e.g. TLS) namespaced by `--agg-name`; hidden unless `SPTXINSIGHT_EXPERIMENTAL=1`. |

Typical end-to-end run:

```bash
sptxinsight run \
  -i ./samples \                      # dir of .h5ad, or sptx-list:///list.txt
  -o ./results \                      # local dir or s3://bucket/prefix
  --base-type tumor --target-type lymphocyte
```

---

## 4. Cellular Microenvironments (CME / niches)

A **CME** (or *niche*) is a recurring local cell mixture — a tumor core, an
immune-infiltrated rim, a stromal band. `sptxinsight cme` discovers them
unsupervised: it builds per-sample Delaunay cell graphs, gathers k-hop
composition features, trains a global DGI encoder, clusters the embeddings, and
writes a one-hot `cme_<n>` label per cell. Run it **after** `ingest`/`run`.

```bash
# Cell-type niches across all ingested samples:
sptxinsight cme -o ./results --cme-clusters 8 --cme-k-hops 3

# Gene-expression niches (k-hop mean expression):
sptxinsight cme -o ./results --cme-mode expression --cme-batch-correct center

# Feed every gene to the encoder instead of PCA-reduced expression:
sptxinsight cme -o ./results --cme-mode expression --disable-pca
```

`--cme-mode` selects what drives the niches and namespaces the outputs:

| `--cme-mode` | features | output folder | one-hot columns |
| --- | --- | --- | --- |
| `celltype` (default) | k-hop cell-type composition | `cme-outputs-csv/` | `cme_<n>` |
| `expression` | k-hop mean gene expression (`expr_`) | `cme-gex-outputs-csv/` | `gexcme_<n>` |
| `both` | composition + expression (fused) | `cme-hybrid-outputs-csv/` | `hcme_<n>` |

Run the command once per mode to get **parallel** niche families on the same
cells; `celltype` stays byte-identical to earlier releases.

### 4.1 Shared PCA of expression features (default ON)

For `expression` / `both` modes the per-cell gene panel is **reduced to a shared
set of principal components before the k-hop aggregation** (default
`--cme-pca-components 50`). The basis is fit **once on the pooled cohort** with
an `IncrementalPCA` and applied identically to every sample, which:

- denoises the sparse gene panel,
- shrinks the encoder input (50·hops instead of n_genes·hops), and
- keeps niches comparable across samples.

PCA only affects the **encoder input** — the interpretable `expr_` columns are
kept for `cme-profile` markers. The checkpoint files are tagged with the PCA
setting (e.g. `slide-graphs-gex-pca50.joblib`) so toggling PCA never reuses a
stale cache.

- `--cme-pca-components N` — number of shared components (min 2; default 50).
  Effective dimension is `min(N, n_genes)`.
- `--disable-pca` — feed all genes to the encoder instead (no reduction).

### 4.2 Other key options

`--cme-clusters` (KMeans k; omit for an automatic Leiden sweep), `--cme-k-hops`,
`--cme-max-edge-len-um`, `--cme-soft` (probability instead of argmax
composition), `--cme-batch-correct` (`none`/`center`/`harmony` cross-sample
correction of the embeddings — use a **technical** unit such as sample/run as
the batch, never a biological condition), `--cme-regions` (merge cells into
annotation-level regions). `--cme-expression` is a deprecated alias for
`--cme-mode both`.

### 4.3 Naming niches

```bash
sptxinsight cme-profile -o ./results --top-types 5 --top-genes 10
sptxinsight cme-profile -o ./results --cme-mode expression
```

Writes `cme-profile-composition.csv` (mean cell-type fractions per CME) and,
when `expr_` columns are present, `cme-profile-markers.csv` (the `expression`
and `both` modes append a `-gex`/`-hybrid` suffix). When both cell-type and
gene-expression niches exist, `cme-profile` also reports their **agreement**
(normalized mutual information plus a cross-tab) and writes `cme-agreement.csv`.

### 4.4 Niches as an H-Plot axis

```bash
SPTXINSIGHT_EXPERIMENTAL=1 sptxinsight hplot -o ./results \
  --base-type tumor --base-by celltype \
  --target-type 7 --target-by cme        # fraction of cells in cme_7 per layer
```

`--base-by`/`--target-by` accept `celltype` (default), `gene`, or a CME niche
family — `cme`, `cmegex`, or `cmehybrid`. When both axes are CME families they
must be the same family.

---

## 5. Outputs

```
results/
  model-outputs-csv/<id>.csv     # center_x, center_y, prob_<type>, expr_<gene> ...
  graphs/<id>.h5                 # cached Delaunay graph (shared across cme modes)
  slide-graphs-gex-pca50.joblib  # k-hop feature cache, tagged by mode + PCA setting
  cme-outputs-csv/cells/<id>.csv     # per-cell cme_<n> labels (--cme-mode celltype)
  cme-gex-outputs-csv/cells/<id>.csv # per-cell gexcme_<n> labels (--cme-mode expression)
  cme-profile-composition.csv    # per-niche cell-type fractions (after cme-profile)
  cme-profile-markers.csv        # per-niche marker genes (gene-mode only)
  cme-agreement.csv              # celltype-vs-gene niche cross-tab (when both exist)
  hplot-outputs-csv/hplots/...   # per-sample layer curves
  hplot-outputs.csv              # aggregated, gap-filled layer table
```

H-Plot uses WSInsight's naming pattern: `layer` is the signed graph-hop index
from the base-region border, and `distance_um` is the cumulative physical
distance in microns. Per-cell outputs expose hop semantics explicitly through
`hplot_layer` / `signed_distance_to_border_hops`; legacy
`signed_distance_to_border` is only a backward-compatible hop alias.

---

## 6. License

Apache License, Version 2.0.
