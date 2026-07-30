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
| `kurtorank` | `kurtorank>=3.1.0` | Enables KurtoRank-backed annotation command(s), such as `annotate`. |

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
| `annotate` | Run KurtoRank cell-type annotation (Xenium dirs and AnnData inputs). Requires installing `sptxinsight[kurtorank]`; hidden otherwise. |
| `marker-init` | Build a skeleton marker panel from DISCO atlases (KurtoRank wrapper). Requires installing `sptxinsight[kurtorank]`; hidden otherwise. |
| `marker-rerank` | Rerank marker genes against CELLxGENE Census (KurtoRank wrapper). Requires installing `sptxinsight[kurtorank]`; hidden otherwise. |
| `verify` | Verify samples are cell-typed and report per-type counts. |
| `export` | Print the path to the aggregated H-Plot table. |
| `describe` | Emit a JSON schema of every subcommand (for tooling / MCP). |
| `niche` | Discover niches across ingested samples. |
| `niche-profile` | Summarise each niche's cell composition and marker genes to help name niches. |
| `hplot`, `hplot-finalize` | Experimental: run/aggregate H-Plot over ingested CSVs. Hidden unless `SPTXINSIGHT_EXPERIMENTAL=1`. |
| `cci` | Experimental: score per-cell ligand-receptor cell-cell interactions over the distance-pruned Delaunay graph (`--d-max` / `--lambda`). Hidden unless `SPTXINSIGHT_EXPERIMENTAL=1`. |
| `agg` | Experimental: detect density-gated cell-type aggregates (e.g. TLS from T+B cells), namespaced by `--name`; usable as `hplot --target-by aggregate`. Hidden unless `SPTXINSIGHT_EXPERIMENTAL=1`. |

Install with KurtoRank features enabled:

```bash
pip install -e '.[kurtorank]'
```

For wsinsight compatibility, install via the command above (single resolver
pass) instead of installing `kurtorank` standalone first. The extra pins
`numpy<2` and the spatial stack (`zarr<3`, `spatialdata/squidpy` caps) so the
environment remains co-installable with wsinsight.

### Example

```bash
sptxinsight run \
  -s ./samples \                      # dir of .h5ad, or sptx-list:///list.txt
  -o ./results \                      # local dir or s3://bucket/prefix
  --base-type tumor --target-type lymphocyte
```

Input requirements per sample: `adata.obsm["spatial"]` (N×2 micron coordinates)
and a categorical `adata.obs["cell_type"]`. Cloud `-s`/`-o` (`s3://`, `gs://`)
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
sptxinsight run -s sptx-list:///data/manifest.txt -o ./results \
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

## Niches (niche / niches)

A **niche** (niche, or *niche*) is a recurring local cell
mixture — e.g. a tumor core, an immune-infiltrated rim, or a stromal band.
`sptxinsight niche` discovers them unsupervised: it builds per-sample Delaunay
cell graphs, gathers k-hop composition features, trains a global DGI encoder,
clusters the embeddings, and writes a one-hot `niche_<n>` label per cell.

```bash
# Discover niches across all ingested samples (run after `ingest`/`run`):
sptxinsight niche -o ./results

# Fix the number of niches, widen the neighborhood, and merge annotation regions:
sptxinsight niche -o ./results --clusters 8 --k-hops 3 --regions

# Gene-expression niches (k-hop mean expression) instead of cell-type niches:
sptxinsight niche -o ./results --mode expression --batch-correct center

# Same, but feed every gene to the encoder instead of PCA-reduced expression:
sptxinsight niche -o ./results --mode expression --disable-pca
```

`--mode` selects what drives the niches and namespaces the outputs so the
families coexist on the same cells:

| `--mode` | features | output folder | one-hot columns |
| --- | --- | --- | --- |
| `celltype` (default) | k-hop cell-type composition | `niche-outputs-csv/` | `niche_<n>` |
| `expression` | k-hop mean gene expression (`expr_`) | `niche-gex-outputs-csv/` | `gexniche_<n>` |
| `both` | composition + expression (fused) | `niche-hybrid-outputs-csv/` | `hniche_<n>` |

Run the command twice (once per mode) to get **parallel** cell-type and gene
niches on the same cells; `celltype` stays byte-identical to earlier releases.

For `expression`/`both` modes the per-cell gene panel is **reduced to a shared
set of principal components before the k-hop aggregation** (default
`--pca-components 50`). The basis is fit once on the pooled cohort and
applied identically to every sample, which denoises the sparse panel, shrinks
the encoder input, and keeps niches comparable across samples. Pass
`--disable-pca` to feed all genes in instead. PCA only affects the encoder
input — the interpretable `expr_` columns are kept for `niche-profile` markers.

Key options: `--clusters` (KMeans k; omit for an automatic Leiden sweep),
`--k-hops`, `--max-edge-len-um`, `--soft` (probability instead of
argmax composition), `--mode` (`celltype`/`expression`/`both`),
`--epochs` (upper bound on DGI training epochs — early stopping is always
active, so training may finish sooner), `--patience` (consecutive epochs
without a mean-loss improvement > `--min-delta` before stopping; default 20),
`--min-delta` (minimum relative mean-loss improvement to reset patience;
default 1e-4), `--min-epochs` (never stop before this many epochs; default 50),
`--amp` (CUDA automatic mixed precision for DGI training — faster, lower GPU
memory; off by default; no effect on CPU/MPS),
`--pca-components` / `--disable-pca` (shared PCA of expression features),
`--batch-correct` (`none`/`center`/`harmony` cross-sample correction of the
embeddings — use a technical unit such as sample/run as the batch, never a
biological condition), and `--regions` (merge cells into annotation-level
regions). `--expression` is a deprecated alias for `--mode both`.
The `harmony` method needs the optional `harmonypy` extra
(`pip install 'sptxinsight[harmony]'`); `center` needs no extra dependency.

### Naming niches

`niche-profile` turns the bare niche ids into interpretable profiles — the
dominant cell types per niche plus, for gene-mode runs, the top enriched marker
genes:

```bash
sptxinsight niche-profile -o ./results --top-types 5 --top-genes 10

# Profile the gene-expression niches instead of the cell-type ones:
sptxinsight niche-profile -o ./results --mode expression
```

It writes `niche-profile-composition.csv` (mean cell-type fractions per niche) and,
when `expr_` columns are present, `niche-profile-markers.csv` (the `expression` and
`both` modes append a `-gex`/`-hybrid` suffix). When both cell-type and
gene-expression niches exist, `niche-profile` also reports their **agreement** —
the normalized mutual information plus a cross-tab — and writes
`niche-agreement.csv`, showing where the two definitions of "niche" diverge.

### Niches as an H-Plot axis

Once niches exist, the experimental `hplot` subcommand can use a niche as the
**base** region or the **target** quantity via `--base-by` / `--target-by`. The
y-value is then the per-layer **fraction of cells belonging to that niche**:

```bash
SPTXINSIGHT_EXPERIMENTAL=1 sptxinsight hplot -o ./results \
  --base-type tumor --base-by celltype \
  --target-type 7 --target-by niche        # fraction of cells in niche_7 per layer

# Gene-expression niches as the target axis:
SPTXINSIGHT_EXPERIMENTAL=1 sptxinsight hplot -o ./results \
  --base-type tumor --base-by celltype \
  --target-type 3 --target-by nichegex     # fraction of cells in gexniche_3 per layer
```

`--base-by`/`--target-by` accept `celltype` (default), `gene`, or a niche
family — `niche` (cell-type niches), `nichegex` (gene-expression niches), or
`nichehybrid` (fused). When both axes are niche families they must be the same
family. `--base-gene-threshold` applies only to `--base-by gene`. Niche ids may
be given as `7` or `niche_7`/`gexniche_7`.

## Outputs

```
results/
  model-outputs-csv/<id>.csv     # center_x, center_y, prob_<type>, expr_<gene> ...
  graphs/<id>.h5                 # cached Delaunay graph (shared across niche modes)
  niche-outputs-csv/cells/<id>.csv     # per-cell niche_<n> labels (niche --mode celltype)
  niche-gex-outputs-csv/cells/<id>.csv # per-cell gexniche_<n> labels (--mode expression)
  niche-profile-composition.csv    # per-niche cell-type fractions (after `niche-profile`)
  niche-profile-markers.csv        # per-niche marker genes (gene-mode only)
  niche-agreement.csv              # celltype-vs-gene niche cross-tab (when both exist)
  hplot-outputs-csv/hplots/...   # per-sample layer curves
  hplot-outputs.csv              # aggregated, gap-filled layer table
```

H-Plot distances follow the WSInsight contract: `layer` is the signed graph-hop
index from the base-region border, while `distance_um` is the cumulative spatial
distance in microns derived from Delaunay edge lengths. Per-cell H-Plot outputs
also include explicit hop columns such as `hplot_layer` and
`signed_distance_to_border_hops`; legacy `signed_distance_to_border` is kept as a
backward-compatible hop-distance alias.

## License

Apache License, Version 2.0.
