"""`sptxinsight marker-rerank`: rerank marker panels via KurtoRank."""

from __future__ import annotations

from pathlib import Path

import click


@click.command(name="marker-rerank")
@click.option(
    "--input",
    "input_csv",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Input marker CSV. Defaults to KurtoRank bundled markers CSV.",
)
@click.option(
    "--output",
    "output_csv",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Output marker CSV path. Defaults to in-place overwrite of input.",
)
@click.option(
    "--qc-output",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Optional per-gene QC output CSV path.",
)
@click.option(
    "--tissue",
    default=None,
    help="Limit to one tissue_type (debug).",
)
@click.option(
    "--tissues",
    default=None,
    help="Comma-separated tissue_types to process (takes precedence over --tissue).",
)
@click.option("--dry-run", is_flag=True, default=False, help="Compute but do not write outputs.")
@click.option("--seed", type=int, default=1234, show_default=True, help="RNG seed.")
@click.option(
    "--census-version",
    default=None,
    help="Census release tag; default = latest.",
)
@click.option(
    "--census-uri",
    default=None,
    help="Path to a local Census SOMA for faster reranking.",
)
@click.option(
    "--parallel",
    type=int,
    default=4,
    show_default=True,
    help="Number of parallel tissue workers (0 or 1 = sequential).",
)
@click.option("--verbose", is_flag=True, default=False, help="Show per-query Census timing.")
@click.option(
    "--checkpoint",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Write intermediate checkpoint CSV after each tissue.",
)
@click.option(
    "--log-file",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Append logs to this file.",
)
def marker_rerank(
    input_csv: Path | None,
    output_csv: Path | None,
    qc_output: Path | None,
    tissue: str | None,
    tissues: str | None,
    dry_run: bool,
    seed: int,
    census_version: str | None,
    census_uri: str | None,
    parallel: int,
    verbose: bool,
    checkpoint: Path | None,
    log_file: Path | None,
) -> None:
    """Rerank marker genes against CELLxGENE Census (KurtoRank wrapper)."""
    try:
        from kurtorank import rerank_markers
        from kurtorank.markers import default_markers_csv
    except ImportError as exc:
        raise click.ClickException(
            "KurtoRank marker-rerank support is unavailable. Install with: "
            "pip install 'sptxinsight[kurtorank]'"
        ) from exc

    try:
        rerank_markers(
            input_csv=input_csv if input_csv is not None else default_markers_csv(),
            output_csv=output_csv,
            qc_output=qc_output,
            tissue=tissue,
            tissues=tissues,
            dry_run=dry_run,
            seed=seed,
            census_version=census_version,
            census_uri=census_uri,
            parallel=parallel,
            verbose=verbose,
            checkpoint=checkpoint,
            log_file=log_file,
            configure_logging=True,
        )
    except click.ClickException:
        raise
    except Exception as exc:
        raise click.ClickException(f"KurtoRank marker-rerank failed: {exc}") from None
