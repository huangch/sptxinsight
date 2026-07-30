"""Guards for the CLI option-naming convention.

A standalone subcommand must not repeat its own name in its option flags:
it is ``sptxinsight niche --clusters``, never ``sptxinsight niche --niche-clusters``.
Only ``sptxinsight run`` prefixes options with a stage name, because it
orchestrates several subcommands and has to keep their namespaces apart.

This convention has silently regressed before (option prefixes were removed and
later reappeared), so it is pinned here rather than left to review.
"""

from __future__ import annotations

import click
import pytest

from sptxinsight.cli.cli import cli


def _subcommands() -> list[tuple[str, click.Command]]:
    """Every registered subcommand except the ``run`` orchestrator."""
    return [(name, cmd) for name, cmd in cli.commands.items() if name != "run"]


def _flags(cmd: click.Command) -> list[str]:
    """All long-form option strings declared on a command."""
    return [opt for param in cmd.params for opt in getattr(param, "opts", [])
            if opt.startswith("--")]


@pytest.mark.parametrize("name,cmd", _subcommands(), ids=lambda v: v if isinstance(v, str) else "")
def test_standalone_command_does_not_prefix_its_own_options(name, cmd):
    offenders = [flag for flag in _flags(cmd) if flag.startswith(f"--{name}-")]
    assert not offenders, (
        f"`sptxinsight {name}` declares self-prefixed option(s) {offenders}. "
        f"Standalone subcommands take unprefixed flags (e.g. --k, --clusters); "
        f"only `sptxinsight run` carries the stage prefix (--{name}-...)."
    )
