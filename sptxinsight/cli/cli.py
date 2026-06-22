"""Top-level Click group wiring sptxinsight's ingest, annotate, and run commands."""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path
from typing import Any
from typing import Literal

import click

from ..io import set_backend
from .agg import agg
from .annotate import annotate
from .cci import cci
from .cme import cme
from .cme import cme_profile_cmd
from .export import export
from .hplot import hplot
from .hplot import hplot_finalize_cmd
from .ingest import ingest
from .run import run

_logging_levels = ["debug", "info", "warning", "error", "critical"]

# Subcommands hidden unless the user opts into experimental features by setting
# the SPTXINSIGHT_EXPERIMENTAL environment variable. Mirrors WSInsight's gating.
_EXPERIMENTAL_COMMANDS = ("hplot", "hplot-finalize", "cci", "agg")


def _experimental_enabled() -> bool:
    v = os.environ.get("SPTXINSIGHT_EXPERIMENTAL", "").strip().lower()
    return v in {"1", "true", "yes", "on"}


@click.group()
@click.option(
    "--backend",
    default=None,
    help="Backend for loading spatial samples.",
    type=click.Choice(["anndata", "zarr", "spatialdata"]),
)
@click.option(
    "--log-level",
    default="info",
    type=click.Choice(_logging_levels),
    help="Set the loudness of logging.",
)
@click.version_option()
@click.pass_context
def cli(
    ctx: click.Context,
    backend: Literal["anndata"] | Literal["zarr"] | Literal["spatialdata"] | None,
    log_level: str,
) -> None:
    """Run cell-typing and spatial-heterogeneity analysis on spatial transcriptomics."""

    levels = {level: getattr(logging, level.upper()) for level in _logging_levels}
    logging.basicConfig(level=levels[log_level])

    if backend is not None:
        set_backend(backend)

    sub = ctx.invoked_subcommand
    if sub in _EXPERIMENTAL_COMMANDS and not _experimental_enabled():
        raise click.UsageError(
            f"'{sub}' is an experimental sptxinsight subcommand. "
            "Set SPTXINSIGHT_EXPERIMENTAL=1 to enable it."
        )


cli.add_command(run)
cli.add_command(ingest)
cli.add_command(annotate)
cli.add_command(export)
cli.add_command(hplot)
cli.add_command(hplot_finalize_cmd)
cli.add_command(cme)
cli.add_command(cme_profile_cmd)
cli.add_command(cci)
cli.add_command(agg)

# Hide experimental commands from --help unless SPTXINSIGHT_EXPERIMENTAL is set.
# They remain registered so `describe` can emit the full schema; invocation is
# blocked in the group callback above.
if not _experimental_enabled():
    for _name in _EXPERIMENTAL_COMMANDS:
        _cmd = cli.commands.get(_name)
        if _cmd is not None:
            _cmd.hidden = True


def _describe_param(param: click.Parameter) -> dict[str, Any]:
    """Serialise one Click option/argument into a JSON-friendly dict."""
    kind: str
    choices: list[str] = []
    t = param.type
    if isinstance(t, click.Choice):
        kind = "choice"
        choices = list(t.choices)
    elif isinstance(t, click.Path):
        kind = "path"
    elif isinstance(t, click.types.BoolParamType):
        kind = "bool"
    elif isinstance(t, click.types.IntParamType):
        kind = "int"
    elif isinstance(t, click.types.FloatParamType):
        kind = "float"
    else:
        kind = "string"

    default = param.default
    if callable(default):
        default = None
    if isinstance(default, Path):
        default = str(default)
    try:
        json.dumps(default)
    except TypeError:
        default = None

    entry: dict[str, Any] = {
        "name": param.name,
        "kind": kind,
        "required": bool(param.required),
        "default": default,
        "help": (param.help if isinstance(param, click.Option) else "") or "",
        "multiple": bool(getattr(param, "multiple", False)),
        "is_flag": bool(getattr(param, "is_flag", False)),
    }
    if isinstance(param, click.Option):
        entry["param_type"] = "option"
        entry["flags"] = list(param.opts) + list(param.secondary_opts)
    else:
        entry["param_type"] = "argument"
        entry["flags"] = []
    if choices:
        entry["choices"] = choices
    if kind == "path":
        entry["path_file_okay"] = bool(getattr(t, "file_okay", True))
        entry["path_dir_okay"] = bool(getattr(t, "dir_okay", True))
        entry["path_exists"] = bool(getattr(t, "exists", False))
    return entry


def _describe_command(name: str, cmd: click.Command) -> dict[str, Any]:
    ctx = click.Context(cmd, info_name=name)
    params: list[dict[str, Any]] = []
    for p in cmd.get_params(ctx):
        if p.name == "help":
            continue
        params.append(_describe_param(p))
    return {
        "name": name,
        "help": (cmd.help or cmd.short_help or "").strip(),
        "params": params,
    }


@cli.command(name="describe")
@click.option(
    "--output",
    "output_path",
    type=click.Path(dir_okay=False, writable=True),
    default=None,
    help="Write the schema JSON to this file instead of stdout.",
)
def describe_cmd(output_path: str | None) -> None:
    """Emit a machine-readable JSON schema of every sptxinsight subcommand."""
    schema: dict[str, Any] = {"schema_version": 1, "commands": {}}
    for name, cmd in cli.commands.items():
        if name == "describe":
            continue
        schema["commands"][name] = _describe_command(name, cmd)
    payload = json.dumps(schema, indent=2, sort_keys=True)
    if output_path:
        Path(output_path).write_text(payload + "\n", encoding="utf-8")
    else:
        click.echo(payload)
