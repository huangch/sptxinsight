#!/usr/bin/env bash
# conda-setup.sh — create and populate the standalone sptxinsight conda environment.
#
# Usage:  sh ./conda-setup.sh [-n ENV_NAME] [-r|--reset] [-m|--mcp]
#
#   -n | --name  ENV_NAME   Conda environment to use (default: current active env).
#   -r | --reset            Deactivate, remove, recreate, and activate the env.
#                           Without this flag the script skips env creation and
#                           only (re-)installs packages into the existing env.
#   -m | --mcp              Also install fastmcp (MCP server support).
#                           Not installed by default to avoid entangling fastmcp's
#                           jaraco.* dep chain with the main resolution.
#
# NOTE: sptxinsight is also co-installable inside the shared wsinsight env
# via:  conda activate wsinsight && pip install --no-deps -e .
# This script creates a *separate* sptxinsight environment instead.
#
# NOTE: spatialdata / squidpy require numpy>=2 and are INCOMPATIBLE with this
# environment (pinned numpy<2 by pyproject.toml). scanpy/anndata ARE installed
# here, pinned to the last numpy<2-compatible line: anndata>=0.12,<0.13 and
# scanpy<1.11. anndata 0.12 is required to read 0.12-format `annotated.h5ad`
# (older anndata raises IORegistryError on encoding_type='null'); newer scanpy
# (>=1.11) drags in numpy>=2. Do NOT relax these pins.

set -e   # abort on first error

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# ── Argument parsing ──────────────────────────────────────────────────────────
ENV_NAME="${CONDA_DEFAULT_ENV:-}"   # default = current active env
DO_RESET=0
DO_MCP=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        -n|--name)
            if [[ -z "${2:-}" ]]; then
                echo "Error: -n/--name requires an environment name." >&2
                exit 1
            fi
            ENV_NAME="$2"
            shift 2
            ;;
        -r|--reset)
            DO_RESET=1
            shift
            ;;
        -m|--mcp)
            DO_MCP=1
            shift
            ;;
        *)
            echo "Unknown option: $1" >&2
            echo "Usage: sh ./conda-setup.sh [-n ENV_NAME] [-r|--reset] [-m|--mcp]" >&2
            exit 1
            ;;
    esac
done

if [[ -z "$ENV_NAME" ]]; then
    echo "Error: no conda environment specified and no environment is currently active." >&2
    echo "       Use -n ENV_NAME to specify one." >&2
    exit 1
fi

echo "Target conda environment: ${ENV_NAME}  (reset=${DO_RESET}, mcp=${DO_MCP})"

# ── (Re-)create environment ───────────────────────────────────────────────────
CONDA_BASE="$(conda info --base 2>/dev/null || true)"
if [[ -z "${CONDA_BASE}" ]]; then
    for _base in /opt/conda /opt/anaconda3; do
        if [[ -f "${_base}/etc/profile.d/conda.sh" ]]; then
            CONDA_BASE="${_base}"
            break
        fi
    done
fi
if [[ -z "${CONDA_BASE}" || ! -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
    echo "Error: cannot locate conda.sh. Activate conda first or set CONDA_BASE." >&2
    exit 1
fi
source "${CONDA_BASE}/etc/profile.d/conda.sh"

if [[ "$DO_RESET" -eq 1 ]]; then
    conda deactivate
    conda env remove -n "${ENV_NAME}" -y 2>/dev/null || true
    # Python only — sptxinsight uses geopandas via pip+pyogrio (no GDAL binary needed).
    conda create -n "${ENV_NAME}" python=3.11 "setuptools<67" -c conda-forge -y
fi

conda activate "${ENV_NAME}"
pip install --upgrade pip

# ── Pip cache fix (NAS inode quota) ──────────────────────────────────────────
# Redirect pip's wheel cache to /tmp to bypass NAS inode quotas. Exported before
# any purge: `pip cache purge` obeys this variable, so purging first wiped the
# user's global ~/.cache/pip. Shared dir so the sibling repos reuse wheels.
export PIP_CACHE_DIR="${PIP_CACHE_DIR:-/tmp/pip-cache-wsinsight-stack}"

pip install "numpy<2"

# ── Install all sptxinsight dependencies declared in pyproject.toml ───────────
# Installs torch, torchvision, torch_geometric, scanpy, anndata, geopandas,
# aiobotocore, boto3, fsspec, s3fs, gcsfs, igraph, leidenalg, etc.
# zarr<3 and harmonypy are installed via [zarr,harmony] extras below.
pip install -c "${SCRIPT_DIR}/constraints.txt" -e "${SCRIPT_DIR}[zarr,harmony]"

# ── MCP server (fastmcp) ──────────────────────────────────────────────────────
# Optional (-m/--mcp). Installed separately to avoid entangling fastmcp's
# jaraco.* dep chain with the main sptxinsight resolution. Version pins are
# in constraints.txt.
if [ "${DO_MCP}" -eq 1 ]; then
    pip install fastmcp
fi

# ── Smoke test ────────────────────────────────────────────────────────────────
# Hard checks are fatal: a half-installed env must not look like a success.
# The test suite is reported but does not fail the setup.
echo "---- smoke test ----"
SMOKE_FAIL=0
smoke() {                       # smoke <label> <command...>
    label="$1"; shift
    if "$@" >/dev/null 2>&1; then
        printf '  PASS  %s\n' "$label"
    else
        printf '  FAIL  %s\n' "$label"
        SMOKE_FAIL=$((SMOKE_FAIL + 1))
    fi
}

python -c 'import importlib.metadata as m; print("  numpy", m.version("numpy"), "| anndata", m.version("anndata"), "| scanpy", m.version("scanpy"))' || true

smoke "sptxinsight on PATH"  command -v sptxinsight
smoke "sptxinsight --help"   sptxinsight --help
smoke "import sptxinsight"   python -c 'import sptxinsight'
smoke "import scanpy"        python -c 'import scanpy'
smoke "numpy < 2"            python -c 'import numpy, sys; sys.exit(int(numpy.__version__.split(".")[0]) >= 2)'
# 0.12 is required to read h5ad written in the 0.12 format.
smoke "anndata >= 0.12"      python -c 'import sys, importlib.metadata as m; from packaging.version import Version; sys.exit(Version(m.version("anndata")) < Version("0.12"))'
smoke "scanpy < 1.11"        python -c 'import sys, importlib.metadata as m; from packaging.version import Version; sys.exit(Version(m.version("scanpy")) >= Version("1.11"))'
if [ "${DO_MCP}" -eq 1 ]; then
    smoke "sptxinsight-mcp on PATH" command -v sptxinsight-mcp
    smoke "sptxinsight-mcp --help"  sptxinsight-mcp --help
fi

if [ -d "${SCRIPT_DIR}/tests" ]; then
    if python -c "import pytest" >/dev/null 2>&1; then
        python -m pytest "${SCRIPT_DIR}/tests" -q \
            && echo "  PASS  test suite" \
            || echo "  WARN  test suite did not pass (non-fatal)"
    else
        echo "  SKIP  test suite (pytest not installed; pip install -e '.[dev]')"
    fi
fi

if [ "${SMOKE_FAIL}" -ne 0 ]; then
    echo "smoke test: ${SMOKE_FAIL} check(s) FAILED" >&2
    exit 1
fi
echo "smoke test: all checks passed"
