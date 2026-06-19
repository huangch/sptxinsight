# sptxinsight

Cell-typing and spatial-heterogeneity (H-Plot) analysis for spatial
transcriptomics, built as a lightweight sibling of
[WSInsight](https://github.com/huangch/wsinsight).

Where WSInsight ingests whole-slide images, `sptxinsight` ingests AnnData
spatial samples (`.h5ad` / `.zarr`) whose coordinates are already in microns. It
reuses WSInsight's H-Plot engine (vendored under `sptxinsight.insightlib`) but
needs **none** of the heavy perception stack (no torch / tensorflow /
openslide).

- **Python:** 3.11+
- **License:** Apache-2.0
- **Status:** Alpha

## Install

Standalone:

```bash
pip install -e .
```

Inside the shared `wsinsight` conda environment, install without dependencies so
`pip` cannot upgrade the locked `numpy<2` / `zarr<3` / `fsspec` generation that
WSInsight depends on (every runtime dependency is already present there):

```bash
pip install --no-deps -e .
```

Optional extras (see `pyproject.toml` for compatibility caveats):

| Extra | Adds | Note |
|---|---|---|
| `zarr` | `zarr<3` | Read `.zarr` samples in the shared env. |
| `spatialdata` | `spatialdata` | Needs `numpy>=2`/`zarr>=3` — dedicated env only. |
| `scanpy` | `scanpy` | Same `numpy>=2` constraint. |
| `mcp` | `fastmcp>=2.0` | Model Context Protocol server. |

## CLI

```bash
sptxinsight --help
```

Global options apply before the subcommand: `--backend {anndata,zarr,spatialdata}`
selects the sample loader and `--log-level` sets logging verbosity.

| Command | Purpose |
|---|---|
| `run` | Ingest → adapt → H-Plot, end to end. |
| `ingest` | Read samples and write the per-sample H-Plot CSV contract. |
| `annotate` | Verify samples are cell-typed and report per-type counts. |
| `export` | Print the path to the aggregated H-Plot table. |
| `describe` | Emit a JSON schema of every subcommand (for tooling / MCP). |
| `cme` | Discover cellular microenvironments (niches) across ingested samples. |
| `cme-profile` | Summarise each CME's cell composition and marker genes to help name niches. |
| `hplot`, `hplot-finalize` | Experimental: run/aggregate H-Plot over ingested CSVs. Hidden unless `SPTXINSIGHT_EXPERIMENTAL=1`. |

### Example

```bash
sptxinsight run \
  -i ./samples \                      # dir of .h5ad, or sptx-list:///list.txt
  -o ./results \                      # local dir or s3://bucket/prefix
  --base-type tumor --target-type lymphocyte
```

Input requirements per sample: `adata.obsm["spatial"]` (N×2 micron coordinates)
and a categorical `adata.obs["cell_type"]`. Cloud `-i`/`-o` (`s3://`, `gs://`)
use the same `S3_STORAGE_OPTIONS` / `GS_STORAGE_OPTIONS` env JSON and
`SPTXINSIGHT_REMOTE_CACHE_DIR` as the URIPath layer.

### `sptx-list://` manifest

One sample path per line; blank lines and `#` comments are ignored. An optional
second column (TAB- or comma-separated) sets an explicit **sample id** — use it
when files share a stem (e.g. Xenium exports every project as `cells.h5ad`),
otherwise the id defaults to the file stem.

```text
# path                              <TAB>  sample_id
/data/projX/cells.h5ad	XENIUM_X01
/data/projY/cells.h5ad	XENIUM_Y02
/data/other/sample_C.h5ad                  # no 2nd column -> id = "sample_C"
```

```bash
sptxinsight run -i sptx-list:///data/manifest.txt -o ./results \
  --base-type tumor --target-type lymphocyte
```

### Gene-expression H-Plots

By default `--base-type` / `--target-type` name **cell types** and the H-Plot
y-value (`target_prop`) is the target-type **proportion** per layer. Switch
either axis to **genes** with `--base-type-by gene` / `--target-type-by gene`;
the listed names are then gene symbols and the y-value becomes the target
gene's **mean expression** per layer (`target_count` becomes the count of
expressing cells). `--expression-matrix` selects the source — `X` (default),
`raw` (`adata.raw.X`), or a layer name in `adata.layers`.

```bash
# Tumor regions by cell type; plot CD8A expression across the layers:
sptxinsight run -i ./samples -o ./results \
  --base-type tumor \
  --target-type CD8A --target-type-by gene \
  --expression-matrix X

# Define the region by gene too (EPCAM-high), threshold the membership:
sptxinsight run -i ./samples -o ./results \
  --base-type EPCAM --base-type-by gene --base-gene-threshold 1.0 \
  --target-type CD8A --target-type-by gene
```

When a gene mode is active the per-sample contract CSV also carries
`expr_<gene>` columns. Cell-type mode (the default) is unchanged and produces
byte-identical results to before.

## Cellular microenvironments (CME / niches)

A **cellular microenvironment** (CME, or *niche*) is a recurring local cell
mixture — e.g. a tumor core, an immune-infiltrated rim, or a stromal band.
`sptxinsight cme` discovers them unsupervised: it builds per-sample Delaunay
cell graphs, gathers k-hop composition features, trains a global DGI encoder,
clusters the embeddings, and writes a one-hot `cme_<n>` label per cell.

```bash
# Discover niches across all ingested samples (run after `ingest`/`run`):
sptxinsight cme -o ./results

# Fix the number of niches, widen the neighborhood, and merge annotation regions:
sptxinsight cme -o ./results --cme-clusters 8 --cme-k-hops 3 --cme-regions

# Fold k-hop mean gene expression into the features (gene-mode samples only):
sptxinsight cme -o ./results --cme-expression
```

Key options: `--cme-clusters` (KMeans k; omit for an automatic Leiden sweep),
`--cme-k-hops`, `--cme-max-edge-len-um`, `--cme-soft` (probability instead of
argmax composition), `--cme-expression` (add `expr_` features), and
`--cme-regions` (merge cells into annotation-level regions).

### Naming niches

`cme-profile` turns the bare `cme_<n>` ids into interpretable profiles — the
dominant cell types per niche plus, for gene-mode runs, the top enriched marker
genes:

```bash
sptxinsight cme-profile -o ./results --top-types 5 --top-genes 10
```

It writes `cme-profile-composition.csv` (mean cell-type fractions per CME) and,
when `expr_` columns are present, `cme-profile-markers.csv`.

### Niches as an H-Plot axis

Once niches exist, the experimental `hplot` subcommand can use a CME as the
**base** region or the **target** quantity via `--base-by cme` / `--target-by
cme` (alongside the existing `celltype` and `gene` modes). The y-value is then
the per-layer **fraction of cells belonging to that niche**:

```bash
SPTXINSIGHT_EXPERIMENTAL=1 sptxinsight hplot -o ./results \
  --base-type tumor --base-by celltype \
  --target-type 7 --target-by cme        # fraction of cells in cme_7 per layer
```

`--base-by`/`--target-by` accept `celltype` (default), `gene`, or `cme`;
`--base-gene-threshold` applies only to `--base-by gene`. CME ids may be given
as `7` or `cme_7`.

## Outputs

```
results/
  model-outputs-csv/<id>.csv     # center_x, center_y, prob_<type>, expr_<gene> ...
  graphs/<id>.h5                 # cached Delaunay graph
  cme-outputs-csv/cells/<id>.csv # per-cell cme_<n> one-hot labels (after `cme`)
  cme-profile-composition.csv    # per-niche cell-type fractions (after `cme-profile`)
  cme-profile-markers.csv        # per-niche marker genes (gene-mode only)
  hplot-outputs-csv/hplots/...   # per-sample layer curves
  hplot-outputs.csv              # aggregated, gap-filled layer table
```

## License

Apache License, Version 2.0.
