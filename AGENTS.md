# sptxinsight — Agent Guide

Spatial-transcriptomics sibling of WSInsight: ingests AnnData samples (`.h5ad`/`.zarr`, micron coords) → cell-typing + H-Plot spatial heterogeneity + niche discovery. Python >=3.11, Apache-2.0, **Alpha**.

## Environment (read first)

- **Co-installable with the shared `wsinsight` conda env** — that's the intended setup. In a shared env, install with `pip install --no-deps -e .` so pip cannot upgrade the locked `numpy<2` / `zarr<3` / `fsspec<2026` generation that WSInsight depends on.
- Version caps in `pyproject.toml` are **deliberate** (they keep the two packages co-installable). Don't "fix" or loosen them. anndata is `>=0.11,<0.13` — 0.12 is needed to read 0.12-format `annotated.h5ad` (older anndata raises `IORegistryError` on `encoding_type='null'`); the hard floor is stardist→`numpy<2` in the shared env, not anndata.
- **Exception:** the `spatialdata` extra needs `numpy>=2` + `zarr>=3` — INCOMPATIBLE with the shared env. Install it only in a dedicated numpy-2 env.
- Optional extras: `mcp` (fastmcp), `kurtorank` (annotate/marker-* commands), `harmony` (batch correction), `zarr` (read `.zarr` in shared env).
- Standalone env: `sh ./conda-setup.sh -n sptxinsight [-m|--mcp]` — fastmcp is **not** installed by default (the `-m`/`--mcp` flag adds it, matching wsinsight's convention).
- Core needs **no** torch/tensorflow/openslide — niche discovery adds torch + torch_geometric + geopandas/shapely.

## CLI

Entry point `sptxinsight` (global `--backend {anndata,zarr,spatialdata}` precedes the subcommand).

Stable: `run` (ingest → H-Plot end-to-end), `ingest`, `verify`, `export`, `describe`, `niche`, `niche-profile`, plus KurtoRank-backed `annotate` / `marker-init` / `marker-rerank` (hidden without `sptxinsight[kurtorank]`).

Experimental (hidden unless `SPTXINSIGHT_EXPERIMENTAL=1`): `hplot`, `hplot-finalize`, `cci` (ligand-receptor), `agg` (cell-type aggregates, e.g. TLS).

- `niche` runs **after** `ingest`/`run`: Delaunay cell graphs → k-hop composition features → global DGI encoder → Leiden clusters → one-hot `niche_<n>` labels.
- Typical run: `sptxinsight run -s ./samples -o ./results --base-type tumor --target-type lymphocyte`.
- Samples need `adata.obsm["spatial"]` (N×2 microns) + categorical `adata.obs["cell_type"]`.
- `annotate` takes the **exact** Xenium dir (no auto `sample → sample/outs` rewrite).

## MCP server (`sptxinsight-mcp`)

- Entry point `sptxinsight.mcp.__main__:main`; extra `mcp = ["fastmcp>=2.0"]`. stdio by default; `--http` (default port **8766**).
- Unlike WSInsight, tools are registered from the **live** CLI schema (`sptxinsight describe` at startup) — no bundled JSON to keep in sync.
- Long-running commands (`run`, `ingest`, `verify`, `niche` …) return a `job_id`; poll `job_status`/`job_logs`/`cancel_job`. Short commands (`export`, `niche-profile`) run synchronously.
- Adapter (`sptxinsight/mcp/adapters.py`) translates snake_case args → kebab-case `--flags`; no positional args supported.

## Vendored H-Plot engine

- The H-Plot engine is **vendored under `sptxinsight.insightlib`** (reused from WSInsight, not imported from the hplot package). When the upstream engine changes in `hplot`/`wsinsight`, the vendored copy must be updated by hand.

## Tests & lint

- Tests: `python -m pytest tests/` (currently: aggregate, CLI option naming, H-Plot distance contract, niche modes).
- Lint: ruff (config in `pyproject.toml`); pre-commit config present (`.pre-commit-config.yaml`).

## Sibling repos (same ecosystem)

- `wsinsight` — WSI-image sibling (primary runtime; shared conda env; MCP port 8765).
- `hplot` — upstream of the H-Plot stats/plotting core that this repo vendors.
- `clawsight` / `clawpyter` — client-side agent plugins driving these MCP servers / Jupyter.
