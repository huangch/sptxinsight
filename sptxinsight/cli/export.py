"""``sptxinsight export``: report the location of aggregated H-Plot outputs."""

from __future__ import annotations

import click

from ..uri_path import URIPath, URIPathType
from ._common import _STORAGE_KWARGS


@click.command()
@click.option(
    "-o",
    "--results-dir",
    type=URIPathType(exists=True, **_STORAGE_KWARGS),
    required=True,
    help="Results directory produced by `sptxinsight run`.",
)
def export(*, results_dir: URIPath) -> None:
    """Print the path to the aggregated H-Plot table if present."""
    hplot_csv = results_dir / "hplot-outputs.csv"
    if not hplot_csv.exists():
        raise click.ClickException(
            f"{hplot_csv} not found; run `sptxinsight run` or `hplot-finalize` first."
        )
    click.echo(str(hplot_csv))
