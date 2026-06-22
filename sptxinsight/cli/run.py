"""``sptxinsight run``: ingest → adapt → H-Plot, end to end."""

from __future__ import annotations

from typing import List

import click

from ..pipeline import run_hplot
from ..uri_path import URIPath
from ..uri_path import URIPathType
from ._common import _STORAGE_KWARGS
from ._common import csv_to_list
from ._common import enumerate_sample_uris


@click.command()
@click.option(
    "-i",
    "--sptx-dir",
    type=URIPathType(exists=True, **_STORAGE_KWARGS),
    required=True,
    help="Directory of spatial samples (.h5ad/.zarr), or a "
    "sptx-list:///path/to/list.txt URI with one sample path per line "
    "(optional TAB/comma 2nd column sets an explicit sample id).",
)
@click.option(
    "-o",
    "--results-dir",
    type=URIPathType(exists=False, **_STORAGE_KWARGS),
    required=True,
    help="Directory to store results (model-outputs-csv/, graphs/, hplot-outputs.csv).",
)
@click.option(
    "--cell-type-key",
    default="cell_type",
    show_default=True,
    help="Column in adata.obs holding the per-cell type label.",
)
@click.option(
    "--spatial-key",
    default="spatial",
    show_default=True,
    help="Key in adata.obsm holding spatial coordinates (microns).",
)
@click.option(
    "--base-type",
    "base_types",
    callback=csv_to_list,
    default=None,
    help="Base cell type(s) (or gene symbol(s) when --base-type-by gene) "
    "forming the region, e.g. tumor cells.",
)
@click.option(
    "--target-type",
    "target_types",
    callback=csv_to_list,
    default=None,
    help="Target cell type(s) whose layer-wise proportion is computed (or "
    "gene symbol(s) when --target-type-by gene, giving mean expression per layer).",
)
@click.option(
    "--base-type-by",
    type=click.Choice(["celltype", "gene"]),
    default="celltype",
    show_default=True,
    help="Interpret --base-type as cell types or gene symbols. In gene mode a "
    "cell is base-positive when its mean base-gene expression exceeds "
    "--base-gene-threshold.",
)
@click.option(
    "--target-type-by",
    type=click.Choice(["celltype", "gene"]),
    default="celltype",
    show_default=True,
    help="Interpret --target-type as cell types (layer-wise proportion) or "
    "gene symbols (layer-wise mean expression).",
)
@click.option(
    "--expression-matrix",
    default="X",
    show_default=True,
    help="Expression source for gene mode: 'X', 'raw' (adata.raw.X), or a layer "
    "name in adata.layers.",
)
@click.option(
    "--base-gene-threshold",
    default=0.0,
    show_default=True,
    type=float,
    help="With --base-type-by gene, mean base-gene expression strictly above "
    "this value marks a cell as base-positive.",
)
@click.option(
    "--hplot-max-neighbor-distance",
    default=25.0,
    type=click.FloatRange(min=0),
    help="Maximal distance (um) to a neighboring cell.",
)
@click.option(
    "--hplot-k",
    default=2,
    type=click.IntRange(min=0),
    help="Maximal edge distance defining a cell's neighborhood.",
)
@click.option(
    "--hplot-n",
    default=8,
    type=click.IntRange(min=0),
    help="Minimal neighborhood size for tumor-region determination.",
)
@click.option(
    "--hplot-r",
    default=0.5,
    type=click.FloatRange(min=0, max=1),
    help="Minimal ratio of base cells in a neighborhood to include a cell in a region.",
)
@click.option(
    "--hplot-range-max",
    default=None,
    type=click.IntRange(min=1),
    help="Maximal layer index toward OUTSIDE of regions for the H-Plot window.",
)
@click.option(
    "--hplot-range-min",
    default=None,
    type=click.IntRange(max=0),
    help="Minimal layer index toward INSIDE of regions for the H-Plot window.",
)
@click.option(
    "--hplot-samples-with-valid-range-only",
    is_flag=True,
    default=False,
    show_default=True,
    help="Use only samples with a valid range of cellular-wise layers.",
)
@click.option(
    "--num-workers",
    default=8,
    show_default=True,
    type=click.IntRange(min=1),
    help="Number of samples to process concurrently.",
)
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    show_default=True,
    help="Overwrite existing results instead of skipping samples with outputs.",
)
def run(
    *,
    sptx_dir: URIPath,
    results_dir: URIPath,
    cell_type_key: str,
    spatial_key: str,
    base_types: List[str] | None,
    target_types: List[str] | None,
    base_type_by: str,
    target_type_by: str,
    expression_matrix: str,
    base_gene_threshold: float,
    hplot_max_neighbor_distance: float,
    hplot_k: int,
    hplot_n: int,
    hplot_r: float,
    hplot_range_max: int | None,
    hplot_range_min: int | None,
    hplot_samples_with_valid_range_only: bool,
    num_workers: int,
    overwrite: bool,
) -> None:
    """Ingest spatial samples and compute aggregated H-Plot outputs."""
    results_dir.mkdir(parents=True, exist_ok=True)
    samples = enumerate_sample_uris(sptx_dir)
    if not samples:
        raise click.ClickException(f"No .h5ad/.zarr samples found under {sptx_dir}")

    failed = run_hplot(
        samples,
        results_dir,
        base_types,
        target_types,
        cell_type_key=cell_type_key,
        spatial_key=spatial_key,
        base_by=base_type_by,
        target_by=target_type_by,
        expression_matrix=expression_matrix,
        base_gene_threshold=base_gene_threshold,
        max_neighbor_distance_um=hplot_max_neighbor_distance,
        hplot_k=hplot_k,
        hplot_N=hplot_n,
        hplot_R=hplot_r,
        hplot_range_max=hplot_range_max,
        hplot_range_min=hplot_range_min,
        samples_with_valid_range_only=hplot_samples_with_valid_range_only,
        num_workers=num_workers,
        overwrite=overwrite,
    )
    ok = len(samples) - len(failed)
    msg = f"H-Plot complete: {ok}/{len(samples)} sample(s); wrote {results_dir / 'hplot-outputs.csv'}"
    if failed:
        msg += f" ({len(failed)} failed: {failed})"
    click.echo(msg)
