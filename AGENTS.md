# sptxinsight — Agent Guide

Spatial-transcriptomics sibling of WSInsight: ingests AnnData samples (`.h5ad`/`.zarr`, micron coords) → cell-typing + H-Plot spatial heterogeneity + niche discovery. Python >=3.11, Apache-2.0, **Alpha**.

## Environment (read first)

- **Co-installable with the shared `wsinsight` conda env** — that's the intended setup. In a shared env, install with `pip install --no-deps -e .` so pip cannot upgrade the locked `numpy<2` / `zarr<3` / `fsspec<2026` generation that WSInsight depends on.
- Version caps in `pyproject.toml` are **deliberate** (they keep the two packages co-installable). Don't "fix" or loosen them. anndata is `>=0.11,<0.13` — 0.12 is needed to read 0.12-format `annotated.h5ad` (older anndata raises `IORegistryError` on `encoding_type='null'`); the hard floor is stardist→`numpy<2` in the shared env, not anndata.
- **Exception:** the `spatialdata` extra needs `numpy>=2` + `zarr>=3` — INCOMPATIBLE with the shared env. Install it only in a dedicated numpy-2 env.
- Optional extras: `mcp` (fastmcp), `kurtorank` (annotate/marker-* commands). `zarr`, `scanpy`, `harmony`/`harmonypy` are core dependencies — no extra to install.
- Standalone env: `sh ./conda-setup.sh sptxinsight [-m|--mcp] [-d|--dev] [-r|--reset]` — fastmcp is **not** installed by default (the `-m`/`--mcp` flag adds it, matching wsinsight's convention). Add `-d`/`--dev` to also install pytest/pytest-cov/ruff/pre_commit for running the test suite; add `-r`/`--reset` to nuke and recreate the env. Run `./conda-setup.sh --help` for the full CLI.
- Core needs **no** torch/tensorflow/openslide — niche discovery adds torch + torch_geometric + geopandas/shapely.

## CLI

Entry point `sptxinsight` (global `--backend {anndata,zarr,spatialdata}` precedes the subcommand).

Stable: `run` (ingest → H-Plot end-to-end), `ingest`, `verify`, `export`, `schema`, `niche`, `niche-profile`, plus KurtoRank-backed `annotate` / `marker-init` / `marker-rerank` (hidden without `sptxinsight[kurtorank]`).

Experimental (hidden unless `SPTXINSIGHT_EXPERIMENTAL=1`): `hplot`, `hplot-finalize`, `cci` (ligand-receptor), `agg` (cell-type aggregates, e.g. TLS).

- `niche` runs **after** `ingest`/`run`: Delaunay cell graphs → k-hop composition features → global DGI encoder → Leiden clusters → one-hot `niche_<n>` labels.
- Typical run: `sptxinsight run -s ./samples -o ./results --base-type tumor --target-type lymphocyte`.
- Samples need `adata.obsm["spatial"]` (N×2 microns) + categorical `adata.obs["cell_type"]`.
- `annotate` takes the **exact** Xenium dir (no auto `sample → sample/outs` rewrite).

## MCP server (`sptxinsight-mcp`)

- Entry point `sptxinsight.mcp.__main__:main`; extra `mcp = ["fastmcp>=2.0"]`. stdio by default; `--http` (default port **8766**).
- Unlike WSInsight, tools are registered from the **live** CLI schema (`sptxinsight schema` at startup) — no bundled JSON to keep in sync.
- Long-running commands (`run`, `ingest`, `verify`, `niche` …) return a `job_id`; poll `job_status`/`job_logs`/`cancel_job`. Short commands (`export`, `niche-profile`) run synchronously.
- Adapter (`sptxinsight/mcp/adapters.py`) translates snake_case args → kebab-case `--flags`; no positional args supported.

## insightlib/ is a copy of WSInsight's pipeline layer

- `sptxinsight/insightlib/` was **copied from `wsinsight/insightlib/`**, then patched for
  AnnData input (see `io/_wsi_stub.py`, which stubs out the WSI-only helpers, and
  `adapt.py`). When wsinsight fixes something in its `insightlib/`, port it here by hand.
- It is **not** a copy of the `hplot` package, and its `compute_hplot` is a different
  implementation from that package's `HPlot` class. The shared name is a collision.
- The `hplot` package is the **downstream, user-facing** analysis layer: it runs in the
  Jupyter environment (driven by `clawpyter`) on outputs this repo has written. Do not
  make this repo import it — that would invert the layering.

## Tests & lint

- Tests: `python -m pytest tests/` (currently: aggregate, CLI option naming, H-Plot distance contract, niche modes).
- Lint: ruff (config in `pyproject.toml`); pre-commit config present (`.pre-commit-config.yaml`).

## Sibling repos (same ecosystem)

- `wsinsight` — WSI-image sibling (primary runtime; shared conda env; MCP port 8765).
- `hplot` — downstream, user-facing analysis package used from Jupyter; not a dependency of this repo.
- `clawsight` / `clawpyter` — client-side agent plugins driving these MCP servers / Jupyter.
