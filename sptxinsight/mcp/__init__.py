"""sptxinsight MCP (Model Context Protocol) server.

This subpackage exposes the sptxinsight CLI as MCP tools so that AI agents
(e.g. Claude Desktop, VS Code Copilot, the ClawSight Hermes plugin) can
invoke sptxinsight subcommands through the same surface as human users.

Entry point: ``sptxinsight-mcp`` (see :mod:`sptxinsight.mcp.__main__`).
"""

from __future__ import annotations

__all__ = ["build_server"]


def build_server(*args, **kwargs):  # pragma: no cover - thin re-export
    from sptxinsight.mcp.server import build_server as _b

    return _b(*args, **kwargs)
