"""``sptxinsight niche``: niche / niche discovery.

Operates on a results directory whose ``model-outputs-csv/`` was populated by
``sptxinsight ingest`` (or ``run``). For each sample it builds a Delaunay cell
graph, computes k-hop cell-type composition features, trains one shared DGI
encoder across the cohort, clusters the embeddings into recurring
microenvironments, and writes per-cell niche labels.

Outputs written to ``<results-dir>/`` (namespaced by ``--mode``)::

    niche-outputs-csv/cells/<id>.csv       celltype niches (niche_*; default)
    niche-gex-outputs-csv/cells/<id>.csv   gene-expression niches (gexniche_*)
    niche-hybrid-outputs-csv/cells/<id>.csv fused niches (hniche_*)
    <subdir>/niches/<id>.csv               annotation-level merged regions (--regions)

The ``niche_0 .. niche_{K-1}`` one-hot columns in the cells CSV can be fed straight
into ``sptxinsight hplot`` as ``prob_`` columns to plot niche proportion over
distance. Because each mode writes to its own folder/column prefix, cell-type and
gene-expression niches coexist on the same cells and can be compared with
``niche-profile``'s agreement report.
"""

from __future__ import annotations

import click
import pandas as pd

from ..uri_path import URIPath
from ..uri_path import URIPathType
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
    "--clusters",
    "niche_clusters",
    default=None,
    type=click.IntRange(min=2),
    help=(
        "Number of microenvironment clusters (KMeans). When omitted, the "
        "optimal number is chosen automatically via a Leiden sweep."
    ),
)
@click.option(
    "--k-hops",
    "niche_k_hops",
    default=2,
    show_default=True,
    type=click.IntRange(min=0),
    help="Number of neighborhood hops for the composition features.",
)
@click.option(
    "--max-edge-len-um",
    "niche_max_edge_len_um",
    default=25.0,
    show_default=True,
    type=click.FloatRange(min=0),
    help="Maximal Delaunay edge length (um) when building the cell graph.",
)
@click.option(
    "--max-cell-radius-um",
    "niche_max_cell_radius_um",
    default=15.0,
    show_default=True,
    type=click.FloatRange(min=0),
    help="Maximal cell radius (um) used when merging annotation-level regions.",
)
@click.option(
    "--epochs",
    "niche_epochs",
    default=300,
    show_default=True,
    type=click.IntRange(min=1),
    help="Upper bound on DGI encoder training epochs.  Early stopping is always "
    "active, so training may finish sooner (see --patience, "
    "--min-delta and --min-epochs).",
)
@click.option(
    "--patience",
    "niche_patience",
    default=20,
    show_default=True,
    type=click.IntRange(min=1),
    help="Early-stopping patience: stop after this many consecutive epochs "
    "without a mean-loss improvement greater than --min-delta.",
)
@click.option(
    "--min-delta",
    "niche_min_delta",
    default=1e-4,
    show_default=True,
    type=click.FloatRange(min=0),
    help="Minimum relative mean-loss improvement required to reset the "
    "early-stopping patience counter.",
)
@click.option(
    "--min-epochs",
    "niche_min_epochs",
    default=50,
    show_default=True,
    type=click.IntRange(min=1),
    help="Never trigger early stopping before this many epochs have elapsed.",
)
@click.option(
    "--amp",
    "niche_amp",
    is_flag=True,
    default=False,
    show_default=True,
    help="Enable CUDA automatic mixed precision for DGI training (faster, lower "
    "GPU memory).  Off by default.  No effect on CPU/MPS.  Note: FP16 math "
    "changes results very slightly versus full FP32.",
)
@click.option(
    "--soft",
    "niche_soft",
    is_flag=True,
    default=False,
    show_default=True,
    help="Use soft (probability) composition features instead of hard argmax labels.",
)
@click.option(
    "--mode",
    "niche_mode",
    default="celltype",
    show_default=True,
    type=click.Choice(["celltype", "expression", "both"]),
    help="Feature source for the niches and output namespace: 'celltype' "
    "(k-hop cell-type composition -> niche-outputs-csv/, niche_ columns), "
    "'expression' (k-hop mean gene expression -> niche-gex-outputs-csv/, "
    "gexniche_ columns), or 'both' (fused -> niche-hybrid-outputs-csv/, "
    "hniche_ columns). Modes write to separate folders so they coexist.",
)
@click.option(
    "--batch-correct",
    "niche_batch_correct",
    default="none",
    show_default=True,
    type=click.Choice(["none", "center", "harmony"]),
    help="Cross-sample correction of the DGI embeddings before clustering: "
    "'center' (per-sample mean-centering, no extra deps) or 'harmony' "
    "(needs the optional harmonypy package). Use the technical unit "
    "(sample/run), never a biological condition, as the batch.",
)
@click.option(
    "--expression",
    "niche_expression",
    is_flag=True,
    default=False,
    show_default=True,
    help="[deprecated] Alias for --mode both (augment composition with "
    "k-hop mean gene expression). Prefer --mode.",
)
@click.option(
    "--pca-components",
    "niche_pca_components",
    default=50,
    show_default=True,
    type=click.IntRange(min=2),
    help="Number of shared PCA components the expression features are "
    "reduced to before k-hop aggregation (expression/both modes). "
    "Ignored for celltype mode and when --disable-pca is set.",
)
@click.option(
    "--disable-pca",
    is_flag=True,
    default=False,
    show_default=True,
    help="Disable the shared PCA reduction of expression features and feed "
    "all genes into the encoder. PCA is on by default because the raw "
    "gene panel is high-dimensional and redundant.",
)
@click.option(
    "--regions",
    "niche_regions",
    is_flag=True,
    default=False,
    show_default=True,
    help="Also merge per-cell labels into annotation-level regions "
    "(requires the optional geopandas/shapely extra).",
)
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    show_default=True,
    help="Delete cached checkpoints and recompute all niche outputs from scratch.",
)
def niche(
    *,
    results_dir: URIPath,
    niche_clusters: int | None,
    niche_k_hops: int,
    niche_max_edge_len_um: float,
    niche_max_cell_radius_um: float,
    niche_epochs: int,
    niche_patience: int,
    niche_min_delta: float,
    niche_min_epochs: int,
    niche_amp: bool,
    niche_soft: bool,
    niche_mode: str,
    niche_batch_correct: str,
    niche_expression: bool,
    niche_pca_components: int,
    disable_pca: bool,
    niche_regions: bool,
    overwrite: bool,
) -> None:
    """Discover niches across ingested samples."""
    # Imported inside the callback (not at module import) so the heavy torch /
    # torch_geometric stack is only loaded when niche analysis is actually run,
    # keeping `sptxinsight --help` and other subcommands fast.
    from ..insightlib.niche_generation import _NICHE_MODE_SPEC
    from ..insightlib.niche_generation import niche_generation

    # Backward-compat: --expression is an alias for --mode both.
    if niche_expression and niche_mode == "celltype":
        niche_mode = "both"

    slide_paths, mpp_lookup = _slide_paths_from_results(results_dir)

    click.secho("\nRunning niche analysis.\n", fg="green")

    niche_generation(
        wsi_dir=None,
        wsi_paths=slide_paths,
        results_dir=str(results_dir),
        max_edge_len_um=niche_max_edge_len_um,
        max_cell_radius_um=niche_max_cell_radius_um,
        k_hops=niche_k_hops,
        alpha=1.0,
        hidden=64,
        out_dim=32,
        epochs=niche_epochs,
        early_stop_patience=niche_patience,
        early_stop_min_delta=niche_min_delta,
        early_stop_min_epochs=niche_min_epochs,
        niche_cellular=True,
        niche_annotation=niche_regions,
        niche_clustering_k=niche_clusters,
        niche_clustering_resolutions=[0.5, 1.0, 2.0],
        niche_soft_mode=niche_soft,
        niche_mode=niche_mode,
        batch_correct=niche_batch_correct,
        expression_pca=0 if disable_pca else niche_pca_components,
        overwrite=overwrite,
        slide_mpp_lookup=mpp_lookup,
        amp=niche_amp,
    )

    ncells = results_dir / _NICHE_MODE_SPEC[niche_mode]["subdir"] / "cells"
    click.secho(f"\nniche analysis completed. Per-cell labels in {ncells}\n", fg="green")


@click.command(name="niche-profile")
@click.option(
    "-o",
    "--results-dir",
    type=URIPathType(exists=True, **_STORAGE_KWARGS),
    required=True,
    help="Results directory containing niche-outputs-csv/cells/ from `sptxinsight niche`.",
)
@click.option(
    "--mode",
    "niche_mode",
    default="celltype",
    show_default=True,
    type=click.Choice(["celltype", "expression", "both"]),
    help="Which niche family to profile (matches the `niche --mode` run).",
)
@click.option(
    "--top-genes",
    default=10,
    show_default=True,
    type=click.IntRange(min=1),
    help="Number of top enriched marker genes to report per niche.",
)
@click.option(
    "--top-types",
    default=5,
    show_default=True,
    type=click.IntRange(min=1),
    help="Number of top cell types to summarise per niche.",
)
@click.option(
    "--agreement/--no-agreement",
    "agreement",
    default=None,
    help="Also report celltype-vs-gene niche agreement (NMI + cross-tab). "
    "Default: auto when both niche-outputs-csv and niche-gex-outputs-csv exist.",
)
def niche_profile_cmd(
    *,
    results_dir: URIPath,
    niche_mode: str,
    top_genes: int,
    top_types: int,
    agreement: bool | None,
) -> None:
    """Summarise each niche's cell composition and marker genes to help name niches."""
    from ..insightlib.niche_profile import niche_agreement
    from ..insightlib.niche_profile import niche_profile

    comp, markers = niche_profile(
        str(results_dir),
        top_genes=top_genes,
        top_types=top_types,
        write=True,
        mode=niche_mode,
    )

    click.secho("\nniche composition (mean cell-type fractions):\n", fg="green")
    cols = [c for c in ("n_cells", "frac", "top_types") if c in comp.columns]
    with pd.option_context("display.max_colwidth", 80, "display.width", 200):
        click.echo(comp[cols].to_string())

    if markers is not None:
        click.secho("\nTop enriched marker genes per niche:\n", fg="green")
        for niche_id, grp in markers.groupby("niche", sort=False):
            top = ", ".join(
                f"{r.gene}({r.log2_enrichment:+.1f})" for r in grp.itertuples()
            )
            click.echo(f"  {niche_id}: {top}")
    else:
        click.secho(
            "\n(No expr_ columns found; run `sptxinsight niche` on gene-mode samples "
            "for marker-gene fingerprints.)\n",
            fg="yellow",
        )

    click.secho(
        f"\nWrote niche-profile-composition*.csv (and markers, if any) to {results_dir}\n",
        fg="green",
    )

    # ---- celltype-vs-gene niche agreement (auto when both families exist) ----
    if agreement is not False:
        result = niche_agreement(str(results_dir), write=True)
        if result is None:
            if agreement is True:
                click.secho(
                    "\n(Agreement needs both niche-outputs-csv/ and niche-gex-outputs-csv/; "
                    "run `niche --mode celltype` and `niche --mode expression` first.)\n",
                    fg="yellow",
                )
        else:
            nmi, crosstab = result
            click.secho(
                f"\nCelltype-vs-gene niche agreement: NMI = {nmi:.3f}\n"
                f"(0 = independent labelings, 1 = identical). Cross-tab "
                f"(rows=celltype niche, cols=gene niche):\n",
                fg="green",
            )
            with pd.option_context("display.width", 200):
                click.echo(crosstab.to_string())
            click.secho(f"\nWrote niche-agreement.csv to {results_dir}\n", fg="green")
