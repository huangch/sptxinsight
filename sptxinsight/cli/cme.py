"""``sptxinsight cme``: cellular microenvironment (CME) / niche discovery.

Operates on a results directory whose ``model-outputs-csv/`` was populated by
``sptxinsight ingest`` (or ``run``). For each sample it builds a Delaunay cell
graph, computes k-hop cell-type composition features, trains one shared DGI
encoder across the cohort, clusters the embeddings into recurring
microenvironments, and writes per-cell CME labels.

Outputs written to ``<results-dir>/`` (namespaced by ``--cme-mode``)::

    cme-outputs-csv/cells/<id>.csv       celltype niches (cme_*; default)
    cme-gex-outputs-csv/cells/<id>.csv   gene-expression niches (gexcme_*)
    cme-hybrid-outputs-csv/cells/<id>.csv fused niches (hcme_*)
    <subdir>/cmes/<id>.csv               annotation-level merged regions (--cme-regions)

The ``cme_0 .. cme_{K-1}`` one-hot columns in the cells CSV can be fed straight
into ``sptxinsight hplot`` as ``prob_`` columns to plot niche proportion over
distance. Because each mode writes to its own folder/column prefix, cell-type and
gene-expression niches coexist on the same cells and can be compared with
``cme-profile``'s agreement report.
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
@click.option("--cme-mode", default="celltype", show_default=True,
              type=click.Choice(["celltype", "expression", "both"]),
              help="Feature source for the niches and output namespace: 'celltype' "
                   "(k-hop cell-type composition -> cme-outputs-csv/, cme_ columns), "
                   "'expression' (k-hop mean gene expression -> cme-gex-outputs-csv/, "
                   "gexcme_ columns), or 'both' (fused -> cme-hybrid-outputs-csv/, "
                   "hcme_ columns). Modes write to separate folders so they coexist.")
@click.option("--cme-batch-correct", default="none", show_default=True,
              type=click.Choice(["none", "center", "harmony"]),
              help="Cross-sample correction of the DGI embeddings before clustering: "
                   "'center' (per-sample mean-centering, no extra deps) or 'harmony' "
                   "(needs the optional harmonypy package). Use the technical unit "
                   "(sample/run), never a biological condition, as the batch.")
@click.option("--cme-expression", is_flag=True, default=False, show_default=True,
              help="[deprecated] Alias for --cme-mode both (augment composition with "
                   "k-hop mean gene expression). Prefer --cme-mode.")
@click.option("--cme-pca-components", default=50, show_default=True, type=click.IntRange(min=2),
              help="Number of shared PCA components the expression features are "
                   "reduced to before k-hop aggregation (expression/both modes). "
                   "Ignored for celltype mode and when --disable-pca is set.")
@click.option("--disable-pca", is_flag=True, default=False, show_default=True,
              help="Disable the shared PCA reduction of expression features and feed "
                   "all genes into the encoder. PCA is on by default because the raw "
                   "gene panel is high-dimensional and redundant.")
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
    cme_mode: str,
    cme_batch_correct: str,
    cme_expression: bool,
    cme_pca_components: int,
    disable_pca: bool,
    cme_regions: bool,
    overwrite: bool,
) -> None:
    """Discover cellular microenvironments (CMEs) across ingested samples."""
    # Imported inside the callback (not at module import) so the heavy torch /
    # torch_geometric stack is only loaded when CME analysis is actually run,
    # keeping `sptxinsight --help` and other subcommands fast.
    from ..insightlib.cme_generation import _CME_MODE_SPEC, cme_generation

    # Backward-compat: --cme-expression is an alias for --cme-mode both.
    if cme_expression and cme_mode == "celltype":
        cme_mode = "both"

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
        cme_mode=cme_mode,
        batch_correct=cme_batch_correct,
        expression_pca=0 if disable_pca else cme_pca_components,
        overwrite=overwrite,
        slide_mpp_lookup=mpp_lookup,
    )

    ncells = results_dir / _CME_MODE_SPEC[cme_mode]["subdir"] / "cells"
    click.secho(f"\nCME analysis completed. Per-cell labels in {ncells}\n", fg="green")


@click.command(name="cme-profile")
@click.option(
    "-o",
    "--results-dir",
    type=URIPathType(exists=True, **_STORAGE_KWARGS),
    required=True,
    help="Results directory containing cme-outputs-csv/cells/ from `sptxinsight cme`.",
)
@click.option("--cme-mode", default="celltype", show_default=True,
              type=click.Choice(["celltype", "expression", "both"]),
              help="Which niche family to profile (matches the `cme --cme-mode` run).")
@click.option("--top-genes", default=10, show_default=True, type=click.IntRange(min=1),
              help="Number of top enriched marker genes to report per CME.")
@click.option("--top-types", default=5, show_default=True, type=click.IntRange(min=1),
              help="Number of top cell types to summarise per CME.")
@click.option("--agreement/--no-agreement", "agreement", default=None,
              help="Also report celltype-vs-gene niche agreement (NMI + cross-tab). "
                   "Default: auto when both cme-outputs-csv and cme-gex-outputs-csv exist.")
def cme_profile_cmd(*, results_dir: URIPath, cme_mode: str, top_genes: int,
                    top_types: int, agreement: bool | None) -> None:
    """Summarise each CME's cell composition and marker genes to help name niches."""
    from ..insightlib.cme_profile import cme_agreement, cme_profile

    comp, markers = cme_profile(
        str(results_dir), top_genes=top_genes, top_types=top_types, write=True,
        mode=cme_mode,
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

    click.secho(f"\nWrote cme-profile-composition*.csv (and markers, if any) to {results_dir}\n",
                fg="green")

    # ---- celltype-vs-gene niche agreement (auto when both families exist) ----
    if agreement is not False:
        result = cme_agreement(str(results_dir), write=True)
        if result is None:
            if agreement is True:
                click.secho(
                    "\n(Agreement needs both cme-outputs-csv/ and cme-gex-outputs-csv/; "
                    "run `cme --cme-mode celltype` and `cme --cme-mode expression` first.)\n",
                    fg="yellow")
        else:
            nmi, crosstab = result
            click.secho(
                f"\nCelltype-vs-gene niche agreement: NMI = {nmi:.3f}\n"
                f"(0 = independent labelings, 1 = identical). Cross-tab "
                f"(rows=celltype niche, cols=gene niche):\n", fg="green")
            with pd.option_context("display.width", 200):
                click.echo(crosstab.to_string())
            click.secho(f"\nWrote cme-agreement.csv to {results_dir}\n", fg="green")
