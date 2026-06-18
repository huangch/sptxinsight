# sptxinsight

Cell-typing and spatial-heterogeneity (H-Plot) analysis for spatial
transcriptomics, built as a lightweight sibling of
[WSInsight](https://github.com/huangch/wsinsight).

Where WSInsight ingests whole-slide images, `sptxinsight` ingests AnnData
spatial samples (`.h5ad` / `.zarr`) whose coordinates are already in microns. It
reuses WSInsight's H-Plot engine (vendored under `sptxinsight.insightlib`) but
needs **none** of the heavy perception stack (no torch / tensorflow /
openslide).

## Install

```bash
pip install -e .
```

## CLI

```bash
sptxinsight --help
```

| Command | Purpose |
|---|---|
| `run` | Ingest → adapt → H-Plot, end to end. |
| `ingest` | Read samples and write the per-sample H-Plot CSV contract. |
| `annotate` | Verify samples are cell-typed and report per-type counts. |
| `export` | Print the path to the aggregated H-Plot table. |
| `describe` | Emit a JSON schema of every subcommand. |
| `hplot`, `hplot-finalize` | Experimental: run/aggregate H-Plot over ingested CSVs (set `SPTXINSIGHT_EXPERIMENTAL=1`). |

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

## Outputs

```
results/
  model-outputs-csv/<id>.csv     # center_x, center_y, prob_<type>, expr_<gene> ...
  graphs/<id>.h5                 # cached Delaunay graph
  hplot-outputs-csv/hplots/...   # per-sample layer curves
  hplot-outputs.csv              # aggregated, gap-filled layer table
```
