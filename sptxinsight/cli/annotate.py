"""``sptxinsight annotate``: validate / report per-sample cell typing.

MVP behaviour is pass-through: samples are expected to already carry a cell-type
label in ``adata.obs``. This command verifies the label exists and reports the
per-type counts, providing a hook for a future learned cell-typer.

Verbosity (``-v`` / ``--verbose``, repeatable) controls the aggregate rosters
printed after the per-sample lines: ``-v`` lists the available cell types and
``-vv`` additionally lists the gene panel (for ``--*-type-by gene``).
"""

from __future__ import annotations

from collections import Counter

import click

from ..io import read_sample
from ..uri_path import URIPath, URIPathType
from ._common import _STORAGE_KWARGS, enumerate_sample_uris


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
    "--cell-type-key",
    default="cell_type",
    show_default=True,
    help="Column in adata.obs expected to hold the per-cell type label.",
)
@click.option(
    "-v",
    "--verbose",
    count=True,
    help="Increase output detail. -v lists the available cell types; "
    "-vv additionally lists the gene panel (use with --*-type-by gene).",
)
def annotate(*, sptx_dir: URIPath, cell_type_key: str, verbose: int) -> None:
    """Verify each sample is cell-typed and report per-type counts."""
    samples = enumerate_sample_uris(sptx_dir)
    if not samples:
        raise click.ClickException(f"No .h5ad/.zarr samples found under {sptx_dir}")

    missing: list[str] = []
    all_counts: Counter[str] = Counter()
    all_genes: set[str] = set()
    for uri in samples:
        adata = read_sample(uri)
        all_genes.update(str(g) for g in adata.var_names)
        if cell_type_key not in adata.obs:
            missing.append(uri.sample_id)
            click.echo(f"{uri.sample_id}: MISSING obs[{cell_type_key!r}]")
            continue
        counts = adata.obs[cell_type_key].astype(str).value_counts()
        all_counts.update({str(t): int(n) for t, n in counts.items()})
        types_str = ", ".join(f"{t} ({int(n)})" for t, n in counts.items())
        click.echo(f"{uri.sample_id}: {adata.n_obs} cells, {len(counts)} types: {types_str}")

    n_ok = len(samples) - len(missing)
    if verbose >= 1 and all_counts:
        click.echo("")
        click.echo(
            f"Available cell types across {n_ok} sample(s) "
            f"(use with --base-type / --target-type):"
        )
        for t, n in all_counts.most_common():
            click.echo(f"  {t}: {n} cells")

    if verbose >= 2 and all_genes:
        click.echo("")
        click.echo(
            f"Available genes across {len(samples)} sample(s) "
            f"({len(all_genes)} total; use with --base-type-by / --target-type-by gene):"
        )
        click.echo("  " + ", ".join(sorted(all_genes)))

    if missing:
        raise click.ClickException(
            f"{len(missing)} sample(s) lack obs[{cell_type_key!r}]: {missing}. "
            "Cell-typing is required before H-Plot (no learned typer in MVP)."
        )
