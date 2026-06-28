"""``sptxinsight ingest``: AnnData samples → H-Plot CSV contract."""

from __future__ import annotations

import click

from ..adapt import anndata_to_contract
from ..io import read_sample
from ..uri_path import URIPath
from ..uri_path import URIPathType
from ._common import _STORAGE_KWARGS
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
    help="Directory to store results. Per-sample CSVs are written under "
    "model-outputs-csv/.",
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
def ingest(
    *,
    sptx_dir: URIPath,
    results_dir: URIPath,
    cell_type_key: str,
    spatial_key: str,
) -> None:
    """Read spatial samples and write the per-sample H-Plot CSV contract."""
    results_dir.mkdir(parents=True, exist_ok=True)
    try:
        samples = enumerate_sample_uris(sptx_dir)
    except ValueError as exc:
        raise click.UsageError(str(exc)) from None
    if not samples:
        raise click.ClickException(f"No .h5ad/.zarr samples found under {sptx_dir}")

    for uri in samples:
        slide_id = uri.sample_id
        adata = read_sample(uri)
        types = anndata_to_contract(
            adata,
            slide_id,
            results_dir,
            cell_type_key=cell_type_key,
            spatial_key=spatial_key,
        )
        click.echo(f"{slide_id}: {adata.n_obs} cells, {len(types)} types")

    click.echo(f"Wrote {len(samples)} CSV(s) to {results_dir / 'model-outputs-csv'}")
