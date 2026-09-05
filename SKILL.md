---
name: sptxinsight
description: Operate sptxinsight for spatial-transcriptomics cell typing, niche discovery, H-Plot spatial heterogeneity, and cell-cell interaction scoring on AnnData samples
---

# sptxinsight — Agentic AI Skill File

> **Purpose**: Enable an agentic AI (Claude, OpenClaw, Hermes, or any
> tool-using LLM agent) to autonomously operate `sptxinsight` for
> spatial-transcriptomics cell typing, niche discovery, H-Plot
> spatial-heterogeneity analysis, and ligand-receptor scoring on AnnData
> samples.

---

## 1. What Is sptxinsight?

`sptxinsight` is the **spatial-transcriptomics** sibling of
[WSInsight](https://github.com/huangch/wsinsight). Where WSInsight ingests
whole-slide images, `sptxinsight` ingests AnnData spatial samples
(`.h5ad` / `.zarr`) whose coordinates are already in microns. It reuses
WSInsight's H-Plot engine (vendored under `sptxinsight.insightlib`) but needs
**none** of the heavy perception stack — no torch / tensorflow / openslide for
the core. Niche discovery adds an optional PyTorch-Geometric encoder.

- **Repository**: <https://github.com/huangch/sptxinsight>
- **License**: Apache-2.0
- **Python**: ≥ 3.11
- **Status**: Alpha
- **Entry points**: `sptxinsight` (CLI), `sptxinsight-mcp` (MCP server)

### 1.1 Sample data contract

Every sample must provide:

| Field | Requirement |
| ----- | ----------- |
| `adata.obsm["spatial"]` | N×2 array of **micron** coordinates. Key overridable with `--spatial-key`. |
| `adata.obs["cell_type"]` | Categorical per-cell type label. Key overridable with `--cell-type-key`. |

Gene-expression workflows (`--base-type-by gene`, `--mode expression`, `cci`)
additionally need an expression matrix — `X`, `adata.raw.X`, or a named entry in
`adata.layers`, selected with `--expression-matrix`.

Run `sptxinsight verify -s <dir>` first; it is the cheapest way to confirm a
cohort satisfies this contract before committing to a long run.

---

## 2. Assumed Install

Assume `sptxinsight` is already installed, and verify with
`sptxinsight --version` before anything else. Installing is the fallback, not
the normal path: only when that verification fails do you build an environment.

| Situation | Command |
| --------- | ------- |
| Standalone environment | `sh ./conda-setup.sh sptxinsight [-m\|--mcp]` |
| Into the shared `wsinsight` env | `pip install --no-deps -e .` |

`--no-deps` in the shared environment is mandatory: it stops `pip` from
upgrading the locked `numpy<2` / `zarr<3` / `fsspec` generation that WSInsight
depends on. The version caps in `pyproject.toml` exist to keep the two packages
co-installable — do not loosen them.

Optional extras that change the **command surface**, so check them before
reporting a command as missing:

| Extra | Adds | Effect if absent |
| ----- | ---- | ---------------- |
| `kurtorank` | KurtoRank annotation backend | `annotate`, `marker-init`, `marker-rerank` are not registered at all — they do not appear in `--help` and calling them fails as an unknown command. |
| `mcp` | `fastmcp>=4.0,<5` | `sptxinsight-mcp` is unavailable. |
| `spatialdata` | `--backend spatialdata` loader | Needs `numpy>=2` + `zarr>=3`, which is **incompatible** with the shared `wsinsight` env. Install only in a dedicated environment. |

`zarr`, `scanpy`, `harmony` / `harmonypy` are core dependencies; there is no
extra to install for them.

---

## 3. Environment Variables

Export these **before** the `sptxinsight` command; none of them is a CLI flag.

| Variable | Effect |
| -------- | ------ |
| `SPTXINSIGHT_EXPERIMENTAL` | `1` / `true` / `yes` / `on` unhides and enables `hplot`, `hplot-finalize`, `cci`, `agg`. Without it those four are hidden from `--help` and refuse to run. |
| `SPTXINSIGHT_FORCE_CPU` | Force CPU execution even when CUDA is visible (niche encoder). |
| `SPTXINSIGHT_REMOTE_CACHE_DIR` | Local cache directory used when materialising remote (`s3://`, `gs://`) sample URIs. |

---

## 4. CLI Reference

### 4.1 Command index

Global options precede the subcommand:
`sptxinsight [--backend ...] [--log-level ...] <command> [options]`.

| Command | Availability | Purpose |
| ------- | ------------ | ------- |
| `run` | stable | Ingest samples → adapt → H-Plot, end to end. |
| `verify` | stable | Check samples are cell-typed; report per-type counts. |
| `export` | stable | Print the path to the aggregated H-Plot table. |
| `niche` | stable | Discover recurring local cell mixtures (niches). |
| `niche-profile` | stable | Summarise each niche's composition and markers so you can name it. |
| `schema` | stable | Emit a machine-readable JSON schema of every subcommand. |
| `annotate` | needs `[kurtorank]` | KurtoRank cell-type annotation of a Xenium directory or AnnData file. |
| `marker-init` | needs `[kurtorank]` | Build a skeleton marker panel from DISCO atlases. |
| `marker-rerank` | needs `[kurtorank]` | Rerank marker genes against CELLxGENE Census. |
| `hplot` | experimental | Compute H-Plot layer curves from already-ingested CSVs. |
| `hplot-finalize` | experimental | Aggregate per-sample H-Plot CSVs into `hplot-outputs.csv`. |
| `cci` | experimental | Score per-cell ligand-receptor cell-cell interactions. |
| `agg` | experimental | Detect cell-type aggregates (e.g. TLS), namespaced by `--name`. |

> **There is no `ingest` subcommand.** It was deliberately unregistered. To
> convert an annotated AnnData into the per-cell CSV contract, use
> `wsinsight import --platform xenium-h5ad`, which reads sptxinsight
> `annotated.h5ad` files natively. `sptxinsight run` still ingests internally,
> so a plain end-to-end run needs no separate step.

### 4.2 Global options

| Option | Type | Default | Description |
| ------ | ---- | ------- | ----------- |
| `--backend` | choice | auto | Sample loader: `anndata`, `zarr`, `spatialdata`. |
| `--log-level` | choice | `info` | `debug`, `info`, `warning`, `error`, `critical`. |
| `--version` | flag | — | Print the version and exit. |

### 4.3 `sptxinsight run` — ingest → H-Plot, end to end

```bash
sptxinsight run \
  -s ./samples \
  -o ./results \
  --base-type tumor \
  --target-type lymphocyte
```

`-s` accepts a directory of `.h5ad` / `.zarr`, or a
`sptx-list:///abs/path/list.txt` URI with one sample path per line. `-o`
accepts a local directory or a remote prefix (`s3://bucket/prefix`); it is
created if it does not exist.

| Option | Type | Default | Description |
| ------ | ---- | ------- | ----------- |
| `--sptx-dir / -s` | path/URI | *required* | Directory of samples, or a `sptx-list://` URI. |
| `--results-dir / -o` | path/URI | *required* | Output directory (`model-outputs-csv/`, `graphs/`, `hplot-outputs.csv`). |
| `--cell-type-key` | string | `cell_type` | Column in `adata.obs` holding the per-cell type label. |
| `--spatial-key` | string | `spatial` | Key in `adata.obsm` holding micron coordinates. |
| `--base-type` | csv string | none | Base cell type(s) forming the region, e.g. `tumor`. Comma-separated for several. |
| `--target-type` | csv string | none | Target cell type(s) whose layer-wise proportion is computed. Comma-separated. |
| `--base-type-by` | choice | `celltype` | Interpret `--base-type` as `celltype` or `gene`. |
| `--target-type-by` | choice | `celltype` | Interpret `--target-type` as `celltype` (proportion) or `gene` (mean expression). |
| `--expression-matrix` | string | `X` | Expression source for gene mode: `X`, `raw`, or a layer name. |
| `--base-gene-threshold` | float | `0.0` | With `--base-type-by gene`, mean base-gene expression strictly above this marks a cell base-positive. |
| `--hplot-max-neighbor-distance` | float | `25.0` | Maximum distance (µm) to a neighboring cell. |
| `--hplot-k` | int | `2` | k-hop radius defining a cell's neighborhood. |
| `--hplot-n` | int | `8` | Minimum neighborhood size for region determination. |
| `--hplot-r` | float | `0.5` | Minimum ratio of base cells in a neighborhood to include a cell in a region. |
| `--hplot-range-min` | int | none | Minimum layer index toward the **inside** of regions. |
| `--hplot-range-max` | int | none | Maximum layer index toward the **outside** of regions. |
| `--hplot-samples-with-valid-range-only` | flag | off | Keep only samples whose layer range is valid. |
| `--num-workers` | int | `8` | Samples processed concurrently. |
| `--overwrite` | flag | off | Recompute instead of skipping samples that already have outputs. |

> **Flag-name warning.** `run` prefixes the H-Plot stage options
> (`--hplot-k`, `--hplot-n`, `--hplot-r`,
> `--hplot-max-neighbor-distance`, `--hplot-range-min/max`) because it
> orchestrates several stages and must keep their namespaces apart. The
> standalone `sptxinsight hplot` command uses the **bare** forms (`--k`, `--n`,
> `--r`, `--max-neighbor-distance`, `--range-min/max`). Likewise `run` spells
> the axis selectors `--base-type-by` / `--target-type-by`, while `hplot`
> spells them `--base-by` / `--target-by`. Emitting the wrong form is the most
> common failure with this CLI — check which command you are calling first.

### 4.4 `sptxinsight verify` — check the data contract

```bash
sptxinsight verify -s ./samples -v
```

| Option | Type | Default | Description |
| ------ | ---- | ------- | ----------- |
| `--sptx-dir / -s` | path/URI | *required* | Directory of samples, or a `sptx-list://` URI. |
| `--cell-type-key` | string | `cell_type` | Column expected to hold the per-cell type label. |
| `--verbose / -v` | count | `0` | `-v` lists the available cell types; `-vv` additionally lists the gene panel. |

Run this before `run` or `niche`. It is the only command that reports *why* a
sample would be skipped (missing `cell_type`, missing `spatial`, wrong key
names) without spending the ingest time first.

### 4.5 `sptxinsight export` — locate the aggregated table

```bash
sptxinsight export -o ./results
```

| Option | Type | Default | Description |
| ------ | ---- | ------- | ----------- |
| `--results-dir / -o` | path/URI | *required* | Results directory produced by `sptxinsight run`. |

Prints the path to `hplot-outputs.csv` if present. It does **not** compute
anything; if the file is absent, run `run` (or `hplot` + `hplot-finalize`)
first.

### 4.6 `sptxinsight niche` — unsupervised niche discovery

A **niche** is a recurring local cell mixture — a tumor core, an
immune-infiltrated rim, a stromal band. `niche` discovers them without
supervision: per-sample Delaunay cell graphs → k-hop composition features →
a global DGI encoder → clustering → a one-hot label per cell.

**Run it after `run`** (or after `wsinsight import`); it reads
`model-outputs-csv/` and will fail without it.

```bash
# cell-type niches
sptxinsight niche -o ./results --clusters 8 --k-hops 3

# gene-expression niches
sptxinsight niche -o ./results --mode expression --batch-correct center

# feed every gene to the encoder instead of PCA-reduced expression
sptxinsight niche -o ./results --mode expression --disable-pca
```

| Option | Type | Default | Description |
| ------ | ---- | ------- | ----------- |
| `--results-dir / -o` | path/URI | *required* | Results directory containing `model-outputs-csv/`. |
| `--mode` | choice | `celltype` | Feature source and output namespace: `celltype`, `expression`, `both`. |
| `--clusters` | int | auto | Number of KMeans clusters. Omit for an automatic Leiden-resolution sweep. |
| `--k-hops` | int | `2` | Neighborhood hops used for the composition features. |
| `--max-edge-len-um` | float | `25.0` | Maximum Delaunay edge length (µm) when building the cell graph. |
| `--max-cell-radius-um` | float | `15.0` | Maximum cell radius (µm) used when merging annotation-level regions. |
| `--soft` | flag | off | Use soft (probability) composition features instead of hard argmax labels. |
| `--pca-components` | int | `50` | Shared PCA components the expression features are reduced to before k-hop aggregation. Minimum 2. |
| `--disable-pca` | flag | off | Feed all genes to the encoder instead of the PCA reduction. |
| `--batch-correct` | choice | `none` | Cross-sample correction of embeddings before clustering: `none`, `center`, `harmony`. |
| `--epochs` | int | `300` | **Upper bound** on DGI training epochs; early stopping may finish sooner. |
| `--patience` | int | `20` | Consecutive epochs without a mean-loss improvement > `--min-delta` before stopping. |
| `--min-delta` | float | `1e-4` | Minimum relative improvement that resets the patience counter. |
| `--min-epochs` | int | `50` | Never stop before this many epochs. |
| `--amp` | flag | off | CUDA automatic mixed precision for DGI training. Faster, lower GPU memory; no effect on CPU/MPS. |
| `--regions` | flag | off | Also merge per-cell labels into annotation-level regions. Needs the optional geopandas/shapely extra. |
| `--expression` | flag | off | **Deprecated** alias for `--mode both`. Prefer `--mode`. |
| `--overwrite` | flag | off | Delete cached checkpoints and recompute every niche output. |

`--batch-correct` must be given a **technical** unit (sample, run, instrument)
as the batch. Never pass a biological condition — that regresses out the effect
you are trying to measure.

#### 4.6.1 `--mode` decides everything downstream

| `--mode` | Features | Output folder | One-hot columns |
| -------- | -------- | ------------- | --------------- |
| `celltype` (default) | k-hop cell-type composition | `niche-outputs-csv/` | `niche_<n>` |
| `expression` | k-hop mean gene expression (`expr_`) | `niche-gex-outputs-csv/` | `gexniche_<n>` |
| `both` | composition + expression, fused | `niche-hybrid-outputs-csv/` | `hniche_<n>` |

Run the command once per mode to get **parallel** niche families over the same
cells. `celltype` output stays byte-identical to earlier releases.

#### 4.6.2 Shared PCA of expression features (on by default)

For `expression` / `both`, the per-cell gene panel is reduced to a shared set of
principal components **before** k-hop aggregation. The basis is fit once on the
pooled cohort with an `IncrementalPCA` and applied identically to every sample,
which denoises the sparse panel, shrinks the encoder input to `50·hops` instead
of `n_genes·hops`, and keeps niches comparable across samples.

PCA affects only the encoder input — the interpretable `expr_` columns are kept
for `niche-profile` markers. Checkpoints are tagged with the PCA setting
(`slide-graphs-gex-pca50.joblib`), so toggling PCA never reuses a stale cache.

### 4.7 `sptxinsight niche-profile` — name the niches

```bash
sptxinsight niche-profile -o ./results --top-types 5 --top-genes 10
sptxinsight niche-profile -o ./results --mode expression
```

| Option | Type | Default | Description |
| ------ | ---- | ------- | ----------- |
| `--results-dir / -o` | path/URI | *required* | Results directory containing `niche-outputs-csv/cells/`. |
| `--mode` | choice | `celltype` | Which niche family to profile. Must match the `niche --mode` run. |
| `--top-types` | int | `5` | Top cell types to summarise per niche. |
| `--top-genes` | int | `10` | Top enriched marker genes to report per niche. |
| `--agreement / --no-agreement` | flag | auto | Report celltype-vs-gene niche agreement (NMI + cross-tab). Auto-on when both families exist. |

Writes `niche-profile-composition.csv` and, when `expr_` columns are present,
`niche-profile-markers.csv` (the `expression` and `both` modes append a
`-gex` / `-hybrid` suffix). With both families present it also writes
`niche-agreement.csv`.

### 4.8 `sptxinsight annotate` — KurtoRank cell typing

Requires `sptxinsight[kurtorank]`. Pass the **exact** directory containing the
Xenium output files (`cell_by_feature_matrix.h5`, `cells.csv.gz`, …). The
command does **not** rewrite `-i sample` to `sample/outs`; if your data live
under `outs/`, pass that path.

```bash
sptxinsight annotate \
  -i /data/xenium/sample1/outs \
  -o ./annotated/sample1 \
  --tissue-type breast \
  --use-top-k-markers 25 \
  --n-jobs 8
```

| Option | Type | Default | Description |
| ------ | ---- | ------- | ----------- |
| `--input-path / -i` | path/URI | *required* | Xenium output directory, or an AnnData file. |
| `--output-dir / -o` | path | *required* | Output directory for annotation artifacts. Non-destructive. |
| `--tissue-type` | string | *required* | KurtoRank tissue type, e.g. `bladder`, `lung`, `breast`. |
| `--markers-csv` | path | bundled | Custom marker CSV. Omit to use KurtoRank's bundled `markers-v6.csv`. |
| `--cell-type-key` | string | `cell_type` | Output cell-type column in `adata.obs` (AnnData mode). |
| `--spatial-key` | string | `spatial` | Coordinates key in `adata.obsm` (AnnData mode). |
| `--common-only / --no-common-only` | flag | `--common-only` | Restrict to marker rows flagged common. |
| `--normal-only / --include-cancer` | flag | `--include-cancer` | Include or exclude malignant subtypes. |
| `--include-immune / --no-include-immune` | flag | `--include-immune` | Add shared `immune` / `circulating` marker rows alongside `--tissue-type`. |
| `--use-graphclust / --use-leiden` | flag | `--use-graphclust` | Cluster source: 10x graph clustering or Leiden. |
| `--chosen-leiden-res` | float | `0.5` | Leiden resolution (with `--use-leiden`). |
| `--use-top-k-markers` | int | all | Keep the K most discriminative genes per row (the CSV is stored in atlas-specificity order). |
| `--method` | string (repeatable) | all | Select individual KurtoRank test methods. |
| `--n-perm` | int | `1000` | Permutations for the empirical null. |
| `--n-top-genes` | int | `3500` | Highly-variable genes retained. |
| `--min-cells` | int | `5` | Minimum cells per cluster. |
| `--min-genes` | int | `10` | Minimum genes per cell. |
| `--upper-percentile` | float | `99.0` | Upper clipping percentile. |
| `--lower-percentile` | float | `1.0` | Lower clipping percentile. |
| `--n-jobs` | int | `32` | Parallel workers. |
| `--allow-cuda-parallel / --no-allow-cuda-parallel` | flag | off | Permit `n-jobs > 1` while CUDA is visible. Off by default because duplicate CUDA contexts OOM. |
| `--generate-plots / --no-generate-plots` | flag | on | Write diagnostic plots. |
| `--regenerate-plots / --no-regenerate-plots` | flag | off | Redraw plots from an existing result. |
| `--plot-format` | choice | `png` | `png` or `svg`. |
| `--plot-dpi` | int | `600` | Raster DPI. |
| `--overwrite / --no-overwrite` | flag | off | Recompute existing outputs. |
| `--verbose / --quiet` | flag | `--quiet` | Verbosity. |

### 4.9 `sptxinsight marker-init` — skeleton marker panel from DISCO

Requires `sptxinsight[kurtorank]`. Output is a **skeleton**: it carries
`tissue_type, subtype, markers, n_cells, source, added_at` only. The biology
columns `annotate` needs (`major_type`, `pannuke_label`, `hne_type`,
`hne_label`, `common`, `malignant`) must be filled in by hand. Do not treat it
as a drop-in replacement for the bundled panel.

```bash
sptxinsight marker-init --list-atlases
sptxinsight marker-init --atlases blood,lung -o panel-seed.csv
```

| Option | Type | Default | Description |
| ------ | ---- | ------- | ----------- |
| `--output / -o` | path | *required unless* `--list-atlases` | Output CSV path. |
| `--atlases` | csv string | none | DISCO atlas **slugs** (e.g. `blood`, `adipose_cell`) — not display labels. |
| `--all-atlases` | flag | off | Fetch every atlas matching the type filter. |
| `--list-atlases` | flag | off | Print the atlas catalog and exit without downloading. |
| `--include-disease` | flag | off | Also include atlases with `type == "disease"`. |
| `--include-celltype` | flag | off | Also include atlases with `type == "cell type"`. |
| `--logfc-min` | float | `1.0` | Minimum log fold-change to keep a DEG row. |
| `--pct1-min` | float | `0.25` | Minimum fraction of target-type cells expressing the gene. |
| `--max-markers` | int | none | Per-cell-type cap on marker count (top-N by logfc). |
| `--cache-dir` | path | user cache | Directory for the DISCO HTTP response cache. |
| `--skip-missing` | flag | off | Warn and skip unknown atlas names instead of failing. |
| `--timeout` | int | `30` | HTTP timeout in seconds per request. |

Defaults are Seurat-standard (`logfc >= 1.0` **and** `pct1 >= 0.25`) and the
atlas filter is `type == "tissue"`. One tissue can span several atlases
(`adipose_cell` vs `adipose_nucleus`); merging them is a judgement call, so ask
the user rather than merging silently.

### 4.10 `sptxinsight marker-rerank` — rerank markers against Census

Requires `sptxinsight[kurtorank]`. Reorders the `markers` column of a panel CSV
using CELLxGENE Census statistics and writes `rank_source` / `low_support`
columns back.

```bash
sptxinsight marker-rerank \
  --census-uri /path/to/census-soma \
  --tissues breast,colorectal \
  --parallel 4 \
  --checkpoint checkpoint.csv \
  --log-file rank.log
```

| Option | Type | Default | Description |
| ------ | ---- | ------- | ----------- |
| `--input` | path | bundled panel | Input marker CSV. |
| `--output` | path | in place | Output marker CSV. Defaults to overwriting `--input`. |
| `--qc-output` | path | none | Per-gene QC audit CSV. |
| `--tissues` | csv string | all | Tissue types to process. Takes precedence over `--tissue`. |
| `--tissue` | string | all | Single tissue type (debug). |
| `--census-uri` | path/URI | hosted | Local Census SOMA mirror. Omit to stream over HTTPS. |
| `--census-version` | string | latest | Census release tag. |
| `--parallel` | int | `4` | Parallel tissue workers. `0` or `1` = sequential. |
| `--checkpoint` | path | none | Write an intermediate CSV after each tissue. |
| `--log-file` | path | none | Append logs to this file. |
| `--seed` | int | `1234` | RNG seed. |
| `--dry-run` | flag | off | Compute without writing. |
| `--verbose` | flag | off | Show per-query Census timing. |

Operational notes that materially change runtime:

- One worker runs per tissue, so effective parallelism is
  `min(--parallel, n_tissues)`. Raising it beyond the tissue count does nothing.
- With `--census-uri` (local SOMA) expect roughly 15–25 min for a four-tissue
  subset. Streaming from the hosted Census is **20–30× slower**.
- When streaming, keep `--parallel` at 2–4. More workers thrash the shared HTTPS
  connection and trigger S3 throttling.
- Pin `--census-version` for reproducibility; the default drifts.
- Always pass `--checkpoint` when streaming — transient S3 errors then resume
  cleanly.
- Do not pipe stderr through `tee`; it destroys the live status line. Use
  `--log-file`.

### 4.11 `sptxinsight schema` — machine-readable command surface

```bash
sptxinsight schema                      # JSON to stdout
sptxinsight schema --output cli.json    # JSON to a file
```

| Option | Type | Default | Description |
| ------ | ---- | ------- | ----------- |
| `--output` | path | stdout | Write the schema JSON here instead of stdout. |

Emits `{"schema_version": 1, "commands": {...}}` covering every subcommand
except `schema` itself, **including the experimental ones** regardless of
`SPTXINSIGHT_EXPERIMENTAL`. This is the authoritative source for exact flag
names and defaults — prefer it over guessing when a command's options are
unclear.

### 4.12 Experimental commands

These four are hidden from `--help` and refuse to run unless
`SPTXINSIGHT_EXPERIMENTAL=1` is exported. Their names, defaults, and semantics
are unstable. Do not use them unless the user explicitly asks for that stage.

```
Error: 'hplot' is an experimental sptxinsight subcommand.
Set SPTXINSIGHT_EXPERIMENTAL=1 to enable it.
```

#### `hplot` — layer curves from already-ingested CSVs

```bash
SPTXINSIGHT_EXPERIMENTAL=1 sptxinsight hplot -o ./results \
  --base-type tumor --base-by celltype \
  --target-type 7 --target-by niche
```

| Option | Type | Default | Description |
| ------ | ---- | ------- | ----------- |
| `--results-dir / -o` | path/URI | *required* | Results directory containing `model-outputs-csv/`. |
| `--base-type` | csv string | none | Base cell type(s) / gene(s) / niche id(s) forming the cluster. |
| `--target-type` | csv string | none | Target whose layer-wise value is computed. |
| `--base-by` | choice | `celltype` | `celltype`, `gene`, `niche`, `nichegex`, `nichehybrid`, `cci`, `aggregate`. |
| `--target-by` | choice | `celltype` | Same choices as `--base-by`. |
| `--base-gene-threshold` | float | `0.0` | Mean expression above which a cell counts as base (`--base-by gene` only). |
| `--max-neighbor-distance` | float | `25.0` | Maximum distance (µm) to a neighboring cell. |
| `--k` | int | `2` | k-hop radius defining a cell's neighborhood. |
| `--n` | int | `8` | Minimum neighborhood size for region determination. |
| `--r` | float | `0.5` | Minimum ratio of base cells in a neighborhood. |
| `--range-min` | int | none | Minimum layer index toward the inside. |
| `--range-max` | int | none | Maximum layer index toward the outside. |
| `--samples-with-valid-range-only` | flag | off | Keep only samples with a valid layer range. |
| `--num-workers` | int | `8` | Samples processed concurrently. |
| `--overwrite` | flag | off | Recompute instead of skipping. |

When both axes name a niche family they must be the **same** family — you
cannot put `niche` on one axis and `nichegex` on the other.

#### `hplot-finalize` — aggregate the per-sample curves

| Option | Type | Default | Description |
| ------ | ---- | ------- | ----------- |
| `--results-dir / -o` | path/URI | *required* | Results directory whose per-sample H-Plot CSVs are aggregated. |
| `--overwrite` | flag | off | Overwrite an existing `hplot-outputs.csv`. |

#### `cci` — ligand-receptor interaction scoring

| Option | Type | Default | Description |
| ------ | ---- | ------- | ----------- |
| `--results-dir / -o` | path/URI | *required* | Results directory containing `model-outputs-csv/`. |
| `--lr-pairs` | path | built-in | CSV/TSV of ligand-receptor pairs (`ligand`/`receptor` or `ligand_gene_symbol`/`receptor_gene_symbol`). |
| `--genes` | csv string | all | Keep only pairs whose ligand **and** receptor are in this list. |
| `--d-max` | float | `25.0` | Maximum neighbour distance (µm); prunes the Delaunay graph. |
| `--kernel` | choice | `exponential` | Distance decay: `exponential`, `gaussian`, `binary` (no decay). |
| `--lambda` | float | `25.0` | Decay length (µm). Ignored for `--kernel binary`. |
| `--num-workers` | int | `4` | Samples processed concurrently (threads). |
| `--overwrite` | flag | off | Recompute existing per-sample CCI CSVs. |

#### `agg` — cell-type aggregate detection (e.g. TLS)

| Option | Type | Default | Description |
| ------ | ---- | ------- | ----------- |
| `--results-dir / -o` | path/URI | *required* | Results directory containing `model-outputs-csv/`. |
| `--name` | string | *required* | Product label, lower-case `[a-z0-9_]+`, e.g. `tls`. Namespaces every artifact. |
| `--types` | csv string | *required* | Ingredient cell types, e.g. `t_cell,b_cell`. |
| `--max-neighbor-distance` | float | `25.0` | Maximum Delaunay edge length (µm). |
| `--k` | int | `2` | k-hop radius for the density gate. |
| `--n` | int | `8` | Minimum neighborhood size for membership. |
| `--r` | float | `0.5` | Minimum ingredient-type fraction for membership. |
| `--min-size` | int | `10` | Drop aggregates with fewer cells than this. |
| `--num-workers` | int | `8` | Samples processed concurrently. |
| `--overwrite` | flag | off | Recompute existing per-sample outputs for this `--name`. |

---

## 5. Results Directory Layout

```
results/
  model-outputs-csv/<id>.csv           # center_x, center_y, prob_<type>, expr_<gene> ...
  graphs/<id>.h5                       # cached Delaunay graph (shared across niche modes)
  slide-graphs-gex-pca50.joblib        # k-hop feature cache, tagged by mode + PCA setting
  niche-outputs-csv/cells/<id>.csv     # per-cell niche_<n>      (--mode celltype)
  niche-gex-outputs-csv/cells/<id>.csv # per-cell gexniche_<n>   (--mode expression)
  niche-hybrid-outputs-csv/cells/<id>.csv  # per-cell hniche_<n> (--mode both)
  niche-profile-composition.csv        # per-niche cell-type fractions
  niche-profile-markers.csv            # per-niche marker genes (expression/both modes)
  niche-agreement.csv                  # celltype-vs-gene niche cross-tab (when both exist)
  hplot-outputs-csv/hplots/...         # per-sample layer curves
  hplot-outputs.csv                    # aggregated, gap-filled layer table
```

### 5.1 Layer semantics

H-Plot uses WSInsight's naming: `layer` is the **signed graph-hop index** from
the base-region border (0 = border, positive = inside, negative = outside), and
`distance_um` is the cumulative physical distance in microns.

Per-cell outputs expose hop semantics explicitly through `hplot_layer` and
`signed_distance_to_border_hops`. The legacy column
`signed_distance_to_border` is a backward-compatible **hop** alias — it is not
a micron distance, despite the name. Read `distance_um` when you want microns.

---

## 6. Common Workflows

### 6.1 Annotated cohort → niches → named niches

```bash
sptxinsight verify -s ./samples -v
sptxinsight run -s ./samples -o ./results --base-type tumor --target-type lymphocyte
sptxinsight niche -o ./results --clusters 8 --k-hops 3
sptxinsight niche-profile -o ./results --top-types 5 --top-genes 10
sptxinsight export -o ./results
```

### 6.2 Raw Xenium → annotation → analysis

```bash
sptxinsight annotate -i /data/xenium/s1/outs -o ./annotated/s1 --tissue-type breast
# repeat per sample, then point run at the annotated .h5ad files
sptxinsight run -s ./annotated -o ./results --base-type tumor --target-type t_cell
```

### 6.3 Two parallel niche families and their agreement

```bash
sptxinsight niche -o ./results --mode celltype   --clusters 8
sptxinsight niche -o ./results --mode expression --clusters 8 --batch-correct center
sptxinsight niche-profile -o ./results --mode expression   # writes niche-agreement.csv
```

### 6.4 Gene-defined regions instead of cell-type regions

```bash
sptxinsight run -s ./samples -o ./results \
  --base-type EPCAM --base-type-by gene --base-gene-threshold 0.5 \
  --target-type CD8A --target-type-by gene \
  --expression-matrix raw
```

---

## 7. Troubleshooting

| Symptom | Cause | Fix |
| ------- | ----- | --- |
| `Error: No such command 'annotate'` | The `kurtorank` extra is not installed, so the command is never registered | Install `sptxinsight[kurtorank]`, or use a different annotation source |
| `'hplot' is an experimental sptxinsight subcommand` | `SPTXINSIGHT_EXPERIMENTAL` not set | `export SPTXINSIGHT_EXPERIMENTAL=1` before invoking |
| `Error: No such option: --base-by` | Used the `hplot` spelling on `run` | On `run` the options are `--base-type-by` / `--target-type-by` |
| `Error: No such option: --k` | Used the `hplot` spelling on `run` | On `run` the H-Plot options are prefixed: `--hplot-k`, `--hplot-n`, `--hplot-r` |
| `Error: No such command 'ingest'` | `ingest` is unregistered | Use `wsinsight import --platform xenium-h5ad`, or just `sptxinsight run` |
| Samples silently skipped | Missing / mis-keyed `cell_type` or `spatial` | `sptxinsight verify -s <dir> -v`; then set `--cell-type-key` / `--spatial-key` |
| `niche` fails immediately | `model-outputs-csv/` absent | Run `run` (or `wsinsight import`) first |
| `export` prints nothing useful | `hplot-outputs.csv` not produced yet | Run `run`, or `hplot` + `hplot-finalize` |
| Niches unchanged after toggling PCA | Reused checkpoint | Checkpoints are PCA-tagged, but pass `--overwrite` if you also changed the cohort |
| Niches track sample identity, not biology | Batch effect | `--batch-correct center` (or `harmony`); use a technical batch, never a biological condition |
| `numpy` / `zarr` upgraded and WSInsight broke | Installed without `--no-deps` in the shared env | Reinstall with `pip install --no-deps -e .` |
| CUDA OOM during `annotate` | Duplicate CUDA contexts across workers | Leave `--allow-cuda-parallel` off, or drop `--n-jobs` to 1 |
| `marker-rerank` progress bar garbled | stderr piped through `tee` | Use `--log-file` instead |

---

## 8. Agent Decision Guide

```text
Is sptxinsight installed?  (sptxinsight --version)
├─ Yes → continue
└─ No  → fallback install (§2), then continue

Are the samples cell-typed?
├─ Unknown → sptxinsight verify -s <dir> -v
├─ No, raw Xenium → sptxinsight annotate  (needs [kurtorank])
└─ Yes → continue

What does the user want?
├─ Layer curves vs a tissue boundary → sptxinsight run -s ... -o ...
│        --base-type <region> --target-type <readout>
│        then sptxinsight export -o ...
├─ Recurring local cell mixtures → sptxinsight run first, then
│        sptxinsight niche -o ... → sptxinsight niche-profile -o ...
├─ Marker panel work → marker-init (build) → marker-rerank (order)
└─ Ligand-receptor / aggregates → experimental; confirm with the user, then
         export SPTXINSIGHT_EXPERIMENTAL=1
```

### Key constraints for agents

1. **Check the option spelling against the command.** `run` uses
   `--hplot-*` prefixes and `--base-type-by`; `hplot` uses bare `--k` / `--n` /
   `--r` and `--base-by`. When unsure, read `sptxinsight schema` rather than
   guessing.
2. **`niche` and `hplot` read `model-outputs-csv/`.** Verify it exists and is
   non-empty for the target samples before calling them; if it is missing or
   partial, run `run` (or `wsinsight import`) for the missing samples first
   rather than proceeding.
3. **`--overwrite` is needed to recompute.** Without it, samples with existing
   outputs are skipped, which makes runs idempotent and resumable.
4. **Environment variables are exported, never passed as flags.**
5. **Experimental commands need explicit user intent.** Do not export
   `SPTXINSIGHT_EXPERIMENTAL=1` on your own initiative.
6. **In the shared `wsinsight` environment, never install without `--no-deps`.**
7. **`marker-init` output is a skeleton.** Do not feed it to `annotate` or merge
   it into the bundled panel without hand-curated biology columns.
8. **`--base-type` / `--target-type` are comma-separated lists**, not repeatable
   flags. Pass `--base-type tumor,epithelial`, not `--base-type tumor
   --base-type epithelial`.

---

## 9. MCP Server (`sptxinsight-mcp`)

Requires the `mcp` extra.

```bash
sptxinsight-mcp                              # stdio (default)
sptxinsight-mcp --http 127.0.0.1:8766        # streamable HTTP, loopback only
sptxinsight-mcp --experimental               # also expose hplot/hplot-finalize/cci/agg
sptxinsight-mcp --max-concurrent 2
```

| Option | Default | Description |
| ------ | ------- | ----------- |
| `--http HOST:PORT` | stdio | Serve over streamable HTTP instead of stdio. Suggested port 8766 (wsinsight uses 8765, hplot 8767). |
| `--experimental` | off | Also register `hplot`, `hplot-finalize`, `cci`, `agg`. |
| `--max-concurrent N` | server default | Cap on simultaneously running jobs. |

Tools are registered from the **live** Click command tree — the same data
`sptxinsight schema` serialises — so the MCP surface can never drift from the
CLI. There is no committed JSON to keep in sync.

**Job model.** `run`, `annotate`, `marker-init`, `marker-rerank`, `verify`,
`niche`, `hplot`, `cci`, `agg` are long-running: the tool returns a `job_id`
immediately. Poll `job_status` until it reports a terminal state, use
`job_logs` to stream output, and `cancel_job` to stop. On failure, fetch
`job_logs` before retrying, and never re-issue the same job while a prior one
is still running. `export`, `hplot-finalize` and `niche-profile` run
synchronously and return their result directly.

The adapter translates snake_case arguments to kebab-case flags
(`results_dir` → `--results-dir`). **No positional arguments are supported** —
every parameter is flag-based.

---

## 10. License

Apache License, Version 2.0.
