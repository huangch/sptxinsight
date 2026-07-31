"""``sptxinsight annotate``: run KurtoRank cell-type annotation."""

from __future__ import annotations

from pathlib import Path

import click

from ..integrations.kurtorank_adapter import AnnotateConfig
from ..integrations.kurtorank_adapter import run_kurtorank_annotate_input
from ..uri_path import URIPath
from ..uri_path import URIPathType
from ._common import _STORAGE_KWARGS


def _apply_options(*option_decorators):
    def _decorator(func):
        for option in reversed(option_decorators):
            func = option(func)
        return func

    return _decorator


_input_output_options = _apply_options(
    click.option(
        "-i",
        "--input-path",
        type=URIPathType(exists=True, **_STORAGE_KWARGS),
        required=True,
        help="Input dataset root. Supports local Xenium directories containing "
        "the output files directly, and AnnData inputs "
        "(.h5ad/.zarr directory or sptx-list manifest).",
    ),
    click.option(
        "-o",
        "--output-dir",
        type=click.Path(path_type=Path, exists=False, file_okay=False, dir_okay=True),
        required=True,
        help="Output directory for annotation artifacts (non-destructive).",
    ),
    click.option(
        "--tissue-type",
        required=True,
        help="KurtoRank tissue type (e.g., bladder, lung, breast).",
    ),
)


_annndata_options = _apply_options(
    click.option(
        "--spatial-key",
        default="spatial",
        show_default=True,
        help="Coordinates key in adata.obsm for AnnData mode.",
    ),
    click.option(
        "--cell-type-key",
        default="cell_type",
        show_default=True,
        help="Output cell-type column in adata.obs for AnnData mode.",
    ),
)


_kurtorank_options = _apply_options(
    click.option(
        "--markers-csv",
        type=click.Path(path_type=Path, exists=True, file_okay=True, dir_okay=False),
        default=None,
        help="Optional custom KurtoRank markers CSV.",
    ),
    click.option("--common-only/--no-common-only", default=True, show_default=True),
    click.option("--normal-only/--include-cancer", default=False, show_default=True),
    click.option("--use-graphclust/--use-leiden", default=True, show_default=True),
    click.option("--chosen-leiden-res", default=0.5, show_default=True, type=float),
    click.option("--min-genes", default=10, show_default=True, type=int),
    click.option("--min-cells", default=5, show_default=True, type=int),
    click.option("--lower-percentile", default=1.0, show_default=True, type=float),
    click.option("--upper-percentile", default=99.0, show_default=True, type=float),
    click.option("--n-top-genes", default=3500, show_default=True, type=int),
    click.option("--n-perm", default=1000, show_default=True, type=int),
    click.option("--n-jobs", default=32, show_default=True, type=int),
    click.option(
        "--allow-cuda-parallel/--no-allow-cuda-parallel",
        default=False,
        show_default=True,
    ),
    click.option(
        "--method",
        "methods",
        multiple=True,
        help="Repeatable KurtoRank method selector.",
    ),
    click.option("--generate-plots/--no-generate-plots", default=True, show_default=True),
    click.option("--overwrite/--no-overwrite", default=False, show_default=True),
    click.option(
        "--regenerate-plots/--no-regenerate-plots",
        default=False,
        show_default=True,
    ),
    click.option(
        "--plot-format",
        type=click.Choice(["png", "svg"], case_sensitive=False),
        default="png",
        show_default=True,
    ),
    click.option("--plot-dpi", default=600, show_default=True, type=int),
    click.option("--use-top-k-markers", default=None, type=int),
    click.option("--verbose/--quiet", default=False, show_default=True),
)


@_input_output_options
@_annndata_options
@_kurtorank_options
@click.command(name="annotate")
def annotate(
    *,
    input_path: URIPath,
    output_dir: Path,
    tissue_type: str,
    spatial_key: str,
    cell_type_key: str,
    markers_csv: Path | None,
    common_only: bool,
    normal_only: bool,
    use_graphclust: bool,
    chosen_leiden_res: float,
    min_genes: int,
    min_cells: int,
    lower_percentile: float,
    upper_percentile: float,
    n_top_genes: int,
    n_perm: int,
    n_jobs: int,
    allow_cuda_parallel: bool,
    methods: tuple[str, ...],
    generate_plots: bool,
    overwrite: bool,
    regenerate_plots: bool,
    plot_format: str,
    plot_dpi: int,
    use_top_k_markers: int | None,
    verbose: bool,
) -> None:
    """Run KurtoRank annotation for supported inputs."""
    config = AnnotateConfig(
        tissue_type=tissue_type,
        spatial_key=spatial_key,
        cell_type_key=cell_type_key,
        markers_csv=markers_csv,
        common_only=common_only,
        normal_only=normal_only,
        use_graphclust=use_graphclust,
        chosen_leiden_res=chosen_leiden_res,
        min_genes=min_genes,
        min_cells=min_cells,
        lower_percentile=lower_percentile,
        upper_percentile=upper_percentile,
        n_top_genes=n_top_genes,
        n_perm=n_perm,
        n_jobs=n_jobs,
        allow_cuda_parallel=allow_cuda_parallel,
        methods=methods,
        generate_plots=generate_plots,
        overwrite=overwrite,
        regenerate_plots=regenerate_plots,
        plot_format=plot_format,
        plot_dpi=plot_dpi,
        use_top_k_markers=use_top_k_markers,
        verbose=verbose,
    )

    run_kurtorank_annotate_input(
        input_path=input_path,
        output_dir=output_dir,
        config=config,
    )

    click.echo(f"KurtoRank annotation completed. Outputs written to {output_dir}")
