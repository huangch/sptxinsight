"""``sptxinsight cme``: cellular microenvironment (CME) / niche discovery.

Operates on a results directory whose ``model-outputs-csv/`` was populated by
``sptxinsight ingest`` (or ``run``). For each sample it builds a Delaunay cell
graph, computes k-hop cell-type composition features, trains one shared DGI
encoder across the cohort, clusters the embeddings into recurring
microenvironments, and writes per-cell CME labels.

Outputs written to ``<results-dir>/``::

    cme-outputs-csv/cells/<id>.csv   per-cell CME labels + k-hop features
    cme-outputs-csv/cmes/<id>.csv    annotation-level merged regions (with --cme-regions)

The ``cme_0 .. cme_{K-1}`` one-hot columns in the cells CSV can be fed straight
into ``sptxinsight hplot`` as ``prob_`` columns to plot niche proportion over
distance.
"""

from __future__ import annotations

import click
import pandas as pd

from ..uri_path import URIPath, URIPathType
from ._common import _STORAGE_KWARGS


def _slide_paths_from_results(results_dir: URIPath):
    model_dir = results_dir / "model-outputs-csv"
    if not model_dir.exists():
        raise click.ClickException(
            f"{model_dir} not found; run `sptxinsight ingest` first."
        )
    stems = sorted(p.stem for p in model_dir.iterdir() if p.suffix.lower() == ".csv")
    if not stems:
        raise click.ClickException(f"No model-output CSVs under {model_dir}.")
    slide_paths = [URIPath(f"{s}.h5ad") for s in stems]
    mpp_lookup = {s: 1.0 for s in stems}
    return slide_paths, mpp_lookup


@click.command()
@click.option(
    "-o",
    "--results-dir",
    type=URIPathType(exists=True, **_STORAGE_KWARGS),
    required=True,
    help="Results directory containing model-outputs-csv/ from a prior ingest.",
)
@click.option(
    "--cme-clusters",
    default=None,
    type=click.IntRange(min=2),
    help=(
        "Number of microenvironment clusters (KMeans). When omitted, the "
        "optimal number is chosen automatically via a Leiden sweep."
    ),
)
@click.option("--cme-k-hops", default=2, show_default=True, type=click.IntRange(min=0),
              help="Number of neighborhood hops for the composition features.")
@click.option("--cme-max-edge-len-um", default=25.0, show_default=True, type=click.FloatRange(min=0),
              help="Maximal Delaunay edge length (um) when building the cell graph.")
@click.option("--cme-max-cell-radius-um", default=15.0, show_default=True, type=click.FloatRange(min=0),
              help="Maximal cell radius (um) used when merging annotation-level regions.")
@click.option("--cme-epochs", default=300, show_default=True, type=click.IntRange(min=1),
              help="DGI encoder training epochs.")
@click.option("--cme-soft", is_flag=True, default=False, show_default=True,
              help="Use soft (probability) composition features instead of hard argmax labels.")
@click.option("--cme-expression", is_flag=True, default=False, show_default=True,
              help="Augment composition features with k-hop mean gene expression "
                   "(uses the expr_ columns; only meaningful for gene-mode samples).")
@click.option("--cme-regions", is_flag=True, default=False, show_default=True,
              help="Also merge per-cell labels into annotation-level regions "
                   "(requires the optional geopandas/shapely extra).")
@click.option("--overwrite", is_flag=True, default=False, show_default=True,
              help="Delete cached checkpoints and recompute all CME outputs from scratch.")
def cme(
    *,
    results_dir: URIPath,
    cme_clusters: int | None,
    cme_k_hops: int,
    cme_max_edge_len_um: float,
    cme_max_cell_radius_um: float,
    cme_epochs: int,
    cme_soft: bool,
    cme_expression: bool,
    cme_regions: bool,
    overwrite: bool,
) -> None:
    """Discover cellular microenvironments (CMEs) across ingested samples."""
    # Imported inside the callback (not at module import) so the heavy torch /
    # torch_geometric stack is only loaded when CME analysis is actually run,
    # keeping `sptxinsight --help` and other subcommands fast.
    from ..insightlib.cme_generation import cme_generation

    slide_paths, mpp_lookup = _slide_paths_from_results(results_dir)

    click.secho("\nRunning cellular microenvironment (CME) analysis.\n", fg="green")

    cme_generation(
        wsi_dir=None,
        wsi_paths=slide_paths,
        results_dir=str(results_dir),
        max_edge_len_um=cme_max_edge_len_um,
        max_cell_radius_um=cme_max_cell_radius_um,
        k_hops=cme_k_hops,
        alpha=1.0,
        hidden=64,
        out_dim=32,
        epochs=cme_epochs,
        cme_cellular=True,
        cme_annotation=cme_regions,
        cme_clustering_k=cme_clusters,
        cme_clustering_resolutions=[0.5, 1.0, 2.0],
        cme_soft_mode=cme_soft,
        use_expression=cme_expression,
        overwrite=overwrite,
        slide_mpp_lookup=mpp_lookup,
    )

    ncells = results_dir / "cme-outputs-csv" / "cells"
    click.secho(f"\nCME analysis completed. Per-cell labels in {ncells}\n", fg="green")


@click.command(name="cme-profile")
@click.option(
    "-o",
    "--results-dir",
    type=URIPathType(exists=True, **_STORAGE_KWARGS),
    required=True,
    help="Results directory containing cme-outputs-csv/cells/ from `sptxinsight cme`.",
)
@click.option("--top-genes", default=10, show_default=True, type=click.IntRange(min=1),
              help="Number of top enriched marker genes to report per CME.")
@click.option("--top-types", default=5, show_default=True, type=click.IntRange(min=1),
              help="Number of top cell types to summarise per CME.")
def cme_profile_cmd(*, results_dir: URIPath, top_genes: int, top_types: int) -> None:
    """Summarise each CME's cell composition and marker genes to help name niches."""
    from ..insightlib.cme_profile import cme_profile

    comp, markers = cme_profile(
        str(results_dir), top_genes=top_genes, top_types=top_types, write=True,
    )

    click.secho("\nCME composition (mean cell-type fractions):\n", fg="green")
    cols = [c for c in ("n_cells", "frac", "top_types") if c in comp.columns]
    with pd.option_context("display.max_colwidth", 80, "display.width", 200):
        click.echo(comp[cols].to_string())

    if markers is not None:
        click.secho("\nTop enriched marker genes per CME:\n", fg="green")
        for cme_id, grp in markers.groupby("cme", sort=False):
            top = ", ".join(f"{r.gene}({r.log2_enrichment:+.1f})" for r in grp.itertuples())
            click.echo(f"  {cme_id}: {top}")
    else:
        click.secho(
            "\n(No expr_ columns found; run `sptxinsight cme` on gene-mode samples "
            "for marker-gene fingerprints.)\n", fg="yellow")

    click.secho(f"\nWrote cme-profile-composition.csv (and markers, if any) to {results_dir}\n",
                fg="green")
