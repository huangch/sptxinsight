"""``sptxinsight cci``: cell-cell interaction (ligand-receptor) scoring.

Operates on a results directory whose ``model-outputs-csv/`` was populated by
``sptxinsight ingest`` (or ``run``) and carries per-cell ``expr_*`` transcript
columns. For each sample it builds (or reuses) the cached Delaunay cell graph,
prunes it to a physical distance cutoff, distance-weights the neighbour edges,
and writes per-cell ligand-receptor signalling scores.

Outputs written to ``<results-dir>/``::

    cci-outputs-csv/cells/<id>.csv   per-cell scores
    cci-outputs.csv                  per-sample x per-pair summary

Each LR pair yields four per-cell columns::

    cci_<LIG>_<REC>_out_mean   sender:   ligand_self * mean(neighbour receptor)
    cci_<LIG>_<REC>_out_sum    sender:   ligand_self *  sum(neighbour receptor)
    cci_<LIG>_<REC>_in_mean    receiver: receptor_self * mean(neighbour ligand)
    cci_<LIG>_<REC>_in_sum     receiver: receptor_self *  sum(neighbour ligand)

Neighbours are 1-hop Delaunay cells within ``--d-max`` microns, weighted by a
distance-decay kernel; range is set by ``--d-max`` / ``--lambda`` (no k-hop).
"""

from __future__ import annotations

import click

from ..insightlib.cci_generation import cci_generation
from ..uri_path import URIPath
from ..uri_path import URIPathType
from ._common import _STORAGE_KWARGS
from ._common import csv_to_list


@click.command()
@click.option(
    "-o",
    "--results-dir",
    type=URIPathType(exists=True, **_STORAGE_KWARGS),
    required=True,
    help="Results directory containing model-outputs-csv/ from a prior ingest.",
)
@click.option(
    "--lr-pairs",
    "lr_pairs_path",
    type=URIPathType(exists=True, dir_okay=False, **_STORAGE_KWARGS),
    default=None,
    help="CSV/TSV of ligand-receptor pairs (columns ligand/receptor or "
    "ligand_gene_symbol/receptor_gene_symbol). Defaults to the bundled "
    "human_lr_pair.csv.",
)
@click.option(
    "--genes",
    "restrict_genes",
    callback=csv_to_list,
    default=None,
    help="Comma/space separated gene list; keep only LR pairs whose ligand AND "
    "receptor are in this list (intersected with the sample panel).",
)
@click.option(
    "--d-max",
    "d_max_um",
    default=25.0,
    show_default=True,
    type=click.FloatRange(min=0, min_open=True),
    help="Maximal neighbour distance (microns); prunes the Delaunay graph.",
)
@click.option(
    "--kernel",
    default="exponential",
    show_default=True,
    type=click.Choice(["exponential", "gaussian", "binary"]),
    help="Distance-decay weighting of neighbour edges. 'binary' = no "
    "decay (all neighbours within --d-max count equally).",
)
@click.option(
    "--lambda",
    "lam_um",
    default=25.0,
    show_default=True,
    type=click.FloatRange(min=0, min_open=True),
    help="Decay length (microns) for the exponential/gaussian kernel. "
    "Ignored for --kernel binary.",
)
@click.option(
    "--num-workers",
    default=4,
    show_default=True,
    type=click.IntRange(min=1),
    help="Number of samples to process concurrently (threads).",
)
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    show_default=True,
    help="Recompute and overwrite existing per-sample CCI CSVs.",
)
def cci(
    *,
    results_dir: URIPath,
    lr_pairs_path: URIPath | None,
    restrict_genes: list[str] | None,
    d_max_um: float,
    kernel: str,
    lam_um: float,
    num_workers: int,
    overwrite: bool,
) -> None:
    """Score per-cell ligand-receptor cell-cell interactions (CCI)."""
    stems = (
        sorted(
            p.stem
            for p in (results_dir / "model-outputs-csv").iterdir()
            if p.suffix.lower() == ".csv"
        )
        if (results_dir / "model-outputs-csv").exists()
        else []
    )
    failed = cci_generation(
        results_dir,
        lr_pairs_path=lr_pairs_path,
        restrict_genes=restrict_genes,
        d_max_um=d_max_um,
        kernel=kernel,
        lam_um=lam_um,
        num_workers=num_workers,
        overwrite=overwrite,
    )
    ok = len(stems) - len(failed)
    msg = f"CCI complete: {ok}/{len(stems)} sample(s)."
    if failed:
        msg += f" ({len(failed)} failed: {failed})"
    click.echo(msg)
