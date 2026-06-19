"""``sptxinsight hplot`` / ``hplot-finalize``: H-Plot over already-ingested CSVs.

These operate on a results directory whose ``model-outputs-csv/`` was populated
by ``sptxinsight ingest`` (or ``run``). No spatial samples are re-read; slide ids
and the unit micron spacing are derived from the existing CSV stems.
"""

from __future__ import annotations

from typing import List

import click

from ..insightlib.hplot_generation import hplot_finalize, hplot_generation
from ..uri_path import URIPath, URIPathType
from ._common import _STORAGE_KWARGS, csv_to_list


def _slide_paths_from_results(results_dir: URIPath) -> tuple[list[URIPath], dict[str, float]]:
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
@click.option("--base-type", "base_types", callback=csv_to_list, default=None,
              help="Base cell type(s)/gene(s)/CME id(s) forming the cluster(s).")
@click.option("--target-type", "target_types", callback=csv_to_list, default=None,
              help="Target cell type(s)/gene(s)/CME id(s) whose layer-wise proportion is computed.")
@click.option("--base-by", default="celltype", show_default=True,
              type=click.Choice(["celltype", "gene", "cme", "cmegex", "cmehybrid", "cci"]),
              help="Interpret --base-type as cell types, genes, a CME niche family "
                   "(cme=celltype niches, cmegex=gene niches, cmehybrid=fused), or "
                   "cci (ligand-receptor score columns from `sptxinsight cci`).")
@click.option("--target-by", default="celltype", show_default=True,
              type=click.Choice(["celltype", "gene", "cme", "cmegex", "cmehybrid", "cci"]),
              help="Interpret --target-type as cell types, genes, a CME niche family "
                   "(cme=celltype niches, cmegex=gene niches, cmehybrid=fused), or "
                   "cci (ligand-receptor score columns from `sptxinsight cci`).")
@click.option("--base-gene-threshold", default=0.0, show_default=True, type=float,
              help="Mean expression above which a cell counts as base (only for --base-by gene).")
@click.option("--hplot-max-neighbor-distance", default=25.0, type=click.FloatRange(min=0),
              help="Maximal distance (um) to a neighboring cell.")
@click.option("--hplot-k", default=2, type=click.IntRange(min=0),
              help="Maximal edge distance defining a cell's neighborhood.")
@click.option("--hplot-n", default=8, type=click.IntRange(min=0),
              help="Minimal neighborhood size for region determination.")
@click.option("--hplot-r", default=0.5, type=click.FloatRange(min=0, max=1),
              help="Minimal ratio of base cells in a neighborhood to include a cell.")
@click.option("--hplot-range-max", default=None, type=click.IntRange(min=1),
              help="Maximal layer index toward OUTSIDE for the H-Plot window.")
@click.option("--hplot-range-min", default=None, type=click.IntRange(max=0),
              help="Minimal layer index toward INSIDE for the H-Plot window.")
@click.option("--hplot-samples-with-valid-range-only", is_flag=True, default=False,
              show_default=True, help="Use only samples with a valid layer range.")
@click.option("--num-workers", default=8, show_default=True, type=click.IntRange(min=1),
              help="Number of samples to process concurrently.")
@click.option("--overwrite", is_flag=True, default=False, show_default=True,
              help="Overwrite existing results instead of skipping samples.")
def hplot(
    *,
    results_dir: URIPath,
    base_types: List[str] | None,
    target_types: List[str] | None,
    base_by: str,
    target_by: str,
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
    """Compute H-Plot layer curves from already-ingested model-output CSVs."""
    slide_paths, mpp_lookup = _slide_paths_from_results(results_dir)
    # CME (niche) one-hot columns live in the cme*-outputs-csv/cells/ folder of
    # the matching --cme-mode run (a superset of model-outputs-csv that also
    # carries prob_/expr_ columns). Read from the right family folder whenever an
    # axis is a CME family; both CME axes must use the same family (one file).
    from ..insightlib.hplot_generation import _CME_FAMILY_SUBDIR, _CCI_SUBDIR

    # Axes that read from a non-default cells folder (CME families or CCI). All
    # such axes in one run must agree on a single folder (one file is read).
    _AXIS_SUBDIR = {**_CME_FAMILY_SUBDIR, "cci": _CCI_SUBDIR}
    special_axes = [_AXIS_SUBDIR[b] for b in (base_by, target_by) if b in _AXIS_SUBDIR]
    if special_axes:
        if len(set(special_axes)) > 1:
            raise click.ClickException(
                "--base-by and --target-by select different cells folders "
                f"({base_by!r} vs {target_by!r}); they must read the same file. "
                "Combine a cci/cme axis only with celltype, or with the same family."
            )
        model_output_subdir = special_axes[0]
    else:
        model_output_subdir = "model-outputs-csv"
    failed = hplot_generation(
        wsi_dir=None,
        slide_paths=slide_paths,
        results_dir=results_dir,
        base_type_list=base_types,
        target_type_list=target_types,
        base_by=base_by,
        target_by=target_by,
        base_gene_threshold=base_gene_threshold,
        max_neighbor_distance_um=hplot_max_neighbor_distance,
        hplot_k=hplot_k,
        hplot_N=hplot_n,
        hplot_R=hplot_r,
        hplot_range_max=hplot_range_max,
        hplot_range_min=hplot_range_min,
        hplot_samples_with_valid_range_only=hplot_samples_with_valid_range_only,
        num_workers=num_workers,
        slide_mpp_lookup=mpp_lookup,
        overwrite=overwrite,
        model_output_subdir=model_output_subdir,
    )
    ok = len(slide_paths) - len(failed)
    msg = f"H-Plot complete: {ok}/{len(slide_paths)} sample(s)."
    if failed:
        msg += f" ({len(failed)} failed: {failed})"
    click.echo(msg)


@click.command(name="hplot-finalize")
@click.option(
    "-o",
    "--results-dir",
    type=URIPathType(exists=True, **_STORAGE_KWARGS),
    required=True,
    help="Results directory whose per-sample H-Plot CSVs are aggregated.",
)
@click.option("--overwrite", is_flag=True, default=False, show_default=True,
              help="Overwrite the aggregated hplot-outputs.csv if present.")
def hplot_finalize_cmd(*, results_dir: URIPath, overwrite: bool) -> None:
    """Aggregate per-sample H-Plot CSVs into hplot-outputs.csv."""
    hplot_finalize(results_dir, overwrite=overwrite)
    click.echo(f"Wrote {results_dir / 'hplot-outputs.csv'}")
