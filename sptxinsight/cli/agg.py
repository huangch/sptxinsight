"""``sptxinsight agg``: cell-type aggregate detection over ingested samples.

``sptxinsight agg`` detects connected, density-gated aggregates of a chosen set
of cell types (``--agg-types``) over the cached Delaunay graph, contracts them
into a quotient graph, and names the product via ``--agg-name`` (e.g. ``tls`` for
a T+B aggregate).  The name namespaces every artifact, so multiple ``agg`` runs
on the same sample are additive rather than colliding.

Operates on a results directory whose ``model-outputs-csv/`` was populated by
``sptxinsight ingest`` (or ``run``).  No spatial samples are re-read; sample ids
and the unit micron spacing are derived from the existing CSV stems.
"""

from __future__ import annotations

import math
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List

import click
import pandas as pd
from tqdm import tqdm

from ..insightlib.agg_generation import agg_generation, membership_column
from ..uri_path import URIPath, URIPathType
from ._common import _STORAGE_KWARGS


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


def _coerce_number(token: str) -> int | float | str:
    """Convert CLI token to int/float when possible, otherwise lowercase text."""
    t = token.strip()
    if re.fullmatch(r"[+-]?\d+", t):
        try:
            return int(t)
        except ValueError:
            pass
    try:
        x = float(t)
        if math.isfinite(x):
            return x
    except ValueError:
        pass
    return t.lower()


def _csv_to_list(
    _: click.Context, __: click.Parameter, value: str | list[str] | None
) -> list[int | float | str]:
    if value is None:
        return []
    tokens = (
        value
        if isinstance(value, list)
        else [x for x in re.split(r"[,\s]+", str(value).strip()) if x]
    )
    return [_coerce_number(str(x)) for x in tokens]


def _read_header(csv_path: URIPath) -> set[str]:
    try:
        with csv_path.open("r", encoding="utf-8") as fp:
            return set(pd.read_csv(fp, nrows=0).columns)
    except Exception:
        return set()


def _collect_headers(model_output_dir: URIPath, num_workers: int) -> set[str]:
    """Union of column names across all model-output CSVs."""
    csv_paths = [p for p in model_output_dir.iterdir() if p.name.endswith(".csv")]
    if not csv_paths:
        return set()
    cols: set[str] = set()
    with ThreadPoolExecutor(max_workers=max(1, num_workers)) as ex:
        futures = [ex.submit(_read_header, p) for p in csv_paths]
        for fut in tqdm(
            as_completed(futures),
            total=len(futures),
            desc="Scanning headers",
            unit="csv",
            dynamic_ncols=True,
        ):
            cols.update(fut.result())
    return cols


_NAME_RE = re.compile(r"^[a-z0-9_]+$")


@click.command()
@click.option(
    "-o",
    "--results-dir",
    type=URIPathType(exists=True, **_STORAGE_KWARGS),
    required=True,
    help=(
        "Results directory containing model-outputs-csv/ from a prior ingest. "
        "Must contain a ``model-outputs-csv/`` subdirectory produced by "
        "``sptxinsight ingest`` or ``sptxinsight run``."
    ),
)
@click.option(
    "--agg-name",
    required=True,
    help=(
        "Product label for this aggregate (lower-case [a-z0-9_]+), e.g. 'tls'.  "
        "Namespaces every artifact (object_<name>_prob_<name> column, "
        "agg-<name>-outputs-csv/, agg/<name>/ in the graph cache) and is usable "
        "as a type in `sptxinsight hplot --target-by aggregate`."
    ),
)
@click.option(
    "--agg-types",
    callback=_csv_to_list,
    required=True,
    help="Comma-separated ingredient cell types that may join the aggregate, e.g. 't_cell,b_cell'.",
)
@click.option(
    "--agg-max-neighbor-distance",
    default=25.0,
    type=click.FloatRange(min=0),
    show_default=True,
    help="Maximum distance (um) between neighboring cells in the Delaunay graph.",
)
@click.option(
    "--agg-k",
    default=2,
    type=click.IntRange(min=1),
    show_default=True,
    help="k-hop neighborhood radius for the density gate.",
)
@click.option(
    "--agg-n",
    default=8,
    type=click.IntRange(min=0),
    show_default=True,
    help="Minimum neighborhood size for a cell to be inside an aggregate region.",
)
@click.option(
    "--agg-r",
    default=0.5,
    type=click.FloatRange(min=0, max=1),
    show_default=True,
    help="Minimum fraction of ingredient cells in the neighborhood (density gate).",
)
@click.option(
    "--agg-min-size",
    default=10,
    type=click.IntRange(min=1),
    show_default=True,
    help="Drop aggregates with fewer than this many cells.",
)
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    show_default=True,
    help="Recompute and overwrite existing per-sample agg outputs for this name.",
)
@click.option(
    "--num-workers",
    default=8,
    show_default=True,
    type=click.IntRange(min=1),
    help="Number of samples to process concurrently.",
)
def agg(
    *,
    results_dir: URIPath,
    agg_name: str,
    agg_types: List,
    agg_max_neighbor_distance: float = 25.0,
    agg_k: int = 2,
    agg_n: int = 8,
    agg_r: float = 0.5,
    agg_min_size: int = 10,
    overwrite: bool = False,
    num_workers: int = 8,
) -> None:
    """Detect cell-type aggregates (e.g. TLS) and write namespaced outputs.

    \b
    Outputs written to <results-dir>/:
      model-outputs-csv/<id>.csv         gains object_<name>_prob_<name> column
      agg-<name>-outputs-csv/<id>.csv     one row per aggregate (hplot-consumable)
      graphs/<id>.h5                      gains an agg/<name>/ quotient-graph group
    """
    name = str(agg_name).strip().lower()
    if not _NAME_RE.match(name):
        raise click.ClickException(
            f"--agg-name must match [a-z0-9_]+ (got {agg_name!r}); "
            "lower-case letters, digits and underscores only."
        )

    agg_type_list = [str(t).strip().lower() for t in agg_types if str(t).strip()]
    if not agg_type_list:
        raise click.ClickException("--agg-types must list at least one cell type.")

    slide_paths, mpp_lookup = _slide_paths_from_results(results_dir)

    model_output_dir = results_dir / "model-outputs-csv"

    # --- name-collision guard: cell-type labels, CME ids, region_/object_ tags ---
    headers = _collect_headers(model_output_dir, num_workers)
    cell_types = {
        c.lower().removeprefix("prob_") for c in headers if c.lower().startswith("prob_")
    }
    if name in cell_types:
        raise click.ClickException(
            f"--agg-name '{name}' collides with an existing cell-type label. "
            "Choose a name distinct from the model's cell types."
        )
    own_col = membership_column(name).lower()
    colliding = [
        c for c in headers
        if c.lower().startswith((f"object_{name}_", f"region_{name}_", f"cme_{name}"))
        and c.lower() != own_col
    ]
    if colliding:
        raise click.ClickException(
            f"--agg-name '{name}' collides with existing columns {sorted(colliding)}. "
            "Choose a different name (or pass --overwrite to recompute this aggregate)."
        )

    click.secho(f"\nRunning aggregate (agg) analysis for '{name}'.\n", fg="green")

    failed = agg_generation(
        wsi_dir=None,
        slide_paths=slide_paths,
        results_dir=results_dir,
        name=name,
        agg_types=agg_type_list,
        max_neighbor_distance_um=agg_max_neighbor_distance,
        k=agg_k,
        N=agg_n,
        R=agg_r,
        min_size=agg_min_size,
        num_workers=num_workers,
        slide_mpp_lookup=mpp_lookup,
        overwrite=overwrite,
    )

    if failed:
        click.secho(f"\nagg failed for {len(failed)} sample(s):", fg="yellow")
        click.secho("\n".join(failed), fg="yellow")
    else:
        click.secho("\nagg completed successfully.\n", fg="green")
