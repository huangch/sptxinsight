"""``sptxinsight verify``: validate / report per-sample cell typing.

MVP behaviour is pass-through: samples are expected to already carry a cell-type
label in ``adata.obs``. This command verifies the label exists and reports the
per-type counts, providing a hook for a future learned cell-typer.

Verbosity (``-v`` / ``--verbose``, repeatable) controls the aggregate rosters
printed after the per-sample lines: ``-v`` lists the available cell types and
``-vv`` additionally lists the gene panel (for ``--*-type-by gene``).
"""

from __future__ import annotations

from collections import Counter
from collections import defaultdict

import click

from ..io import read_sample
from ..uri_path import URIPath
from ..uri_path import URIPathType
from ._common import _STORAGE_KWARGS
from ._common import enumerate_sample_uris


def _suggest_obs_keys(obs_columns: list[str]) -> list[str]:
    """Return likely cell-type annotation columns from adata.obs."""
    tokens = ("cell", "type", "annot", "class", "label", "cluster")
    hits = [c for c in obs_columns if any(t in str(c).lower() for t in tokens)]
    return sorted(hits)


@click.command(name="verify")
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
def verify(*, sptx_dir: URIPath, cell_type_key: str, verbose: int) -> None:
    """Verify each sample is cell-typed and report per-type counts."""
    try:
        samples = enumerate_sample_uris(sptx_dir)
    except ValueError as exc:
        raise click.UsageError(str(exc)) from None
    if not samples:
        raise click.ClickException(f"No .h5ad/.zarr samples found under {sptx_dir}")

    missing: list[str] = []
    all_counts: Counter[str] = Counter()
    all_genes: set[str] = set()

    # Manifests with repeated stems (e.g. many "annotated.h5ad") are common;
    # disambiguate output by including full URI when sample ids collide.
    id_counts = Counter(uri.sample_id for uri in samples)
    if any(n > 1 for n in id_counts.values()):
        dup = sorted([sid for sid, n in id_counts.items() if n > 1])
        click.echo(
            "Warning: duplicate sample ids detected from input paths "
            f"{dup}. Output lines will include full paths; for stable ids, "
            "use sptx-list:///... with an explicit TAB/comma sample-id column."
        )

    for uri in samples:
        adata = read_sample(uri)
        all_genes.update(str(g) for g in adata.var_names)
        label = (
            uri.sample_id
            if id_counts[uri.sample_id] == 1
            else f"{uri.sample_id} ({uri})"
        )
        if cell_type_key not in adata.obs:
            obs_cols = [str(c) for c in adata.obs.columns]
            suggestions = _suggest_obs_keys(obs_cols)
            if suggestions:
                hint = f" candidates={suggestions}"
            else:
                hint = f" obs_columns_first30={obs_cols[:30]}"
            missing.append(label)
            click.echo(f"{label}: MISSING obs[{cell_type_key!r}]{hint}")
            continue
        counts = adata.obs[cell_type_key].astype(str).value_counts()
        all_counts.update({str(t): int(n) for t, n in counts.items()})
        types_str = ", ".join(f"{t} ({int(n)})" for t, n in counts.items())
        click.echo(
            f"{label}: {adata.n_obs} cells, {len(counts)} types: {types_str}"
        )

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
            "Cell-typing is required before H-Plot (no learned typer in MVP). "
            "Pass --cell-type-key <column> if labels already exist under a different obs column."
        )
