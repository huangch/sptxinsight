"""Executable entry point: configure multiprocessing then dispatch the CLI."""

from __future__ import annotations

import multiprocessing as mp
import os

import click

from .cancel import install_sigint_handler
from .cli.cli import cli


def main() -> None:
    """Initialize runtime knobs and invoke the Click CLI."""
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
    mp.set_start_method("spawn", force=True)
    install_sigint_handler()

    try:
        cli()
    except (click.Abort, KeyboardInterrupt):
        click.secho("\nsptxinsight: aborted by user.", fg="yellow", err=True)
        raise SystemExit(130) from None


if __name__ == "__main__":
    main()
