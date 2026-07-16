"""`sptxinsight marker-init`: build a skeleton marker panel via KurtoRank."""

from __future__ import annotations

from pathlib import Path

import click


@click.command(name="marker-init")
@click.option(
    "--atlases",
    "atlases_csv",
    default=None,
    help="Comma-separated DISCO atlas slugs (e.g. blood,lung,adipose_cell).",
)
@click.option(
    "--all-atlases",
    is_flag=True,
    default=False,
    help="Fetch every atlas matching the type filter (default: type==tissue).",
)
@click.option(
    "--list-atlases",
    is_flag=True,
    default=False,
    help="Print the DISCO atlas catalog and exit (no download).",
)
@click.option(
    "--include-disease",
    is_flag=True,
    default=False,
    help="Also include atlases where type=='disease' (default: excluded).",
)
@click.option(
    "--include-celltype",
    is_flag=True,
    default=False,
    help="Also include atlases where type=='cell type' (default: excluded).",
)
@click.option(
    "--logfc-min",
    type=float,
    default=1.0,
    show_default=True,
    help="Minimum log fold-change to keep a DEG row as a marker.",
)
@click.option(
    "--pct1-min",
    type=float,
    default=0.25,
    show_default=True,
    help="Minimum fraction of cells in the target cell type that must express the gene.",
)
@click.option(
    "--max-markers",
    type=int,
    default=None,
    help="Optional per-cell-type cap on marker count (top-N by logfc).",
)
@click.option(
    "--output",
    "output_path",
    "-o",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Output CSV path. Required unless --list-atlases is given.",
)
@click.option(
    "--cache-dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=Path.home() / ".kurtorank" / "disco",
    show_default=True,
    help="Directory for the DISCO HTTP response cache.",
)
@click.option(
    "--skip-missing",
    is_flag=True,
    default=False,
    help="Skip (warn) unknown atlas names instead of failing.",
)
@click.option(
    "--timeout",
    type=int,
    default=30,
    show_default=True,
    help="HTTP timeout in seconds per request.",
)
def marker_init(
    atlases_csv: str | None,
    all_atlases: bool,
    list_atlases: bool,
    include_disease: bool,
    include_celltype: bool,
    logfc_min: float,
    pct1_min: float,
    max_markers: int | None,
    output_path: Path | None,
    cache_dir: Path,
    skip_missing: bool,
    timeout: int,
) -> None:
    """Build a skeleton marker panel (KurtoRank build-panel wrapper)."""
    try:
        from kurtorank.seed.main import build_panel_cmd
    except ImportError as exc:
        raise click.ClickException(
            "KurtoRank marker-init support is unavailable. Install with: "
            "pip install 'sptxinsight[kurtorank]'"
        ) from exc

    argv: list[str] = []
    if atlases_csv:
        argv.extend(["--atlases", atlases_csv])
    if all_atlases:
        argv.append("--all-atlases")
    if list_atlases:
        argv.append("--list-atlases")
    if include_disease:
        argv.append("--include-disease")
    if include_celltype:
        argv.append("--include-celltype")
    argv.extend(["--logfc-min", str(logfc_min)])
    argv.extend(["--pct1-min", str(pct1_min)])
    if max_markers is not None:
        argv.extend(["--max-markers", str(max_markers)])
    if output_path is not None:
        argv.extend(["--output", str(output_path)])
    argv.extend(["--cache-dir", str(cache_dir)])
    if skip_missing:
        argv.append("--skip-missing")
    argv.extend(["--timeout", str(timeout)])

    try:
        build_panel_cmd.main(args=argv, prog_name="kurtorank build-panel", standalone_mode=False)
    except click.ClickException:
        raise
    except Exception as exc:
        raise click.ClickException(f"KurtoRank marker-init failed: {exc}") from None
