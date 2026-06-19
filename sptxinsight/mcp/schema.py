"""Schema-driven generation of MCP tool definitions from the sptxinsight CLI.

The single source of truth is the live Click command tree (the same data
``sptxinsight describe`` serialises). This module builds it in-process,
classifies commands as stable / experimental and short / long-running, and
converts each Click parameter into a JSON-schema property suitable for use as
an MCP tool input schema. Building it live (rather than reading a committed
JSON file) means the MCP surface can never drift from the CLI.
"""

from __future__ import annotations

from typing import Any

# Subcommands exposed by default. Mirrors the gating in
# sptxinsight/cli/cli.py (experimental commands are hidden unless the
# SPTXINSIGHT_EXPERIMENTAL env var is set).
STABLE_COMMANDS: frozenset[str] = frozenset(
    {"run", "ingest", "annotate", "export", "cme", "cme-profile"}
)
EXPERIMENTAL_COMMANDS: frozenset[str] = frozenset(
    {"hplot", "hplot-finalize", "cci"}
)

# Commands that may run for many minutes or hours. These are exposed as
# background-job tools in the MCP server (the tool returns a job_id and the
# agent polls job_status / job_logs / cancel_job). All other stable commands
# run synchronously.
LONG_RUNNING_COMMANDS: frozenset[str] = frozenset(
    {"run", "ingest", "annotate", "cme", "hplot", "cci"}
)

_KIND_TO_JSON_TYPE: dict[str, str] = {
    "string": "string",
    "int": "integer",
    "integer": "integer",
    "float": "number",
    "number": "number",
    "bool": "boolean",
    "boolean": "boolean",
    "path": "string",
    "choice": "string",
}


def load_schema() -> dict[str, Any]:
    """Build and return the sptxinsight CLI schema live from the Click app."""
    from sptxinsight.cli.cli import _describe_command, cli

    commands: dict[str, Any] = {}
    for name, cmd in cli.commands.items():
        if name == "describe":
            continue
        commands[name] = _describe_command(name, cmd)
    return {"schema_version": 1, "commands": commands}


def _param_to_json_property(param: dict[str, Any]) -> dict[str, Any]:
    """Convert one Click parameter entry into a JSON-schema property dict."""
    kind = str(param.get("kind", "string")).lower()
    json_type = _KIND_TO_JSON_TYPE.get(kind, "string")
    prop: dict[str, Any] = {"type": json_type}
    help_text = param.get("help")
    if help_text:
        prop["description"] = " ".join(str(help_text).split())
    if param.get("choices"):
        prop["enum"] = list(param["choices"])
    if param.get("multiple"):
        prop = {"type": "array", "items": prop}
        if help_text:
            prop["description"] = " ".join(str(help_text).split())
    default = param.get("default")
    if (
        default is not None
        and not param.get("required", False)
        and not param.get("multiple", False)
    ):
        prop["default"] = default
    return prop


def command_to_input_schema(command: dict[str, Any]) -> dict[str, Any]:
    """Build a JSON-schema ``object`` describing one command's parameters."""
    properties: dict[str, Any] = {}
    required: list[str] = []
    for param in command.get("params", []):
        name = str(param["name"])
        properties[name] = _param_to_json_property(param)
        if param.get("required"):
            required.append(name)
    schema: dict[str, Any] = {
        "type": "object",
        "properties": properties,
        "additionalProperties": False,
    }
    if required:
        schema["required"] = required
    return schema


def discover_commands(experimental: bool = False) -> dict[str, dict[str, Any]]:
    """Return ``{name: command_dict}`` for the commands the server should expose."""
    raw = load_schema().get("commands", {})
    allowed = set(STABLE_COMMANDS)
    if experimental:
        allowed |= set(EXPERIMENTAL_COMMANDS)
    return {name: cmd for name, cmd in raw.items() if name in allowed}


def is_long_running(name: str) -> bool:
    """Return True if the named command should be exposed as a background job."""
    return name in LONG_RUNNING_COMMANDS


__all__ = [
    "STABLE_COMMANDS",
    "EXPERIMENTAL_COMMANDS",
    "LONG_RUNNING_COMMANDS",
    "load_schema",
    "command_to_input_schema",
    "discover_commands",
    "is_long_running",
]
