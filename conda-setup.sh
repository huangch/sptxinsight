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
pip cache purge || true
# Redirect pip's wheel cache to /tmp to bypass NAS inode quotas.
export PIP_CACHE_DIR=/tmp/pip-cache-sptxinsight

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

# ── Safety checks ─────────────────────────────────────────────────────────────
python -c "
import numpy; v = numpy.__version__
assert int(v.split('.')[0]) < 2, f'ERROR: numpy {v} >= 2.0; re-run: pip install \"numpy<2\"'
import importlib.metadata as m
from packaging.version import Version
ad = m.version('anndata'); sp = m.version('scanpy')
assert Version(ad) >= Version('0.12'), f'ERROR: anndata {ad} < 0.12 cannot read 0.12-format h5ad'
assert Version(sp) < Version('1.11'), f'ERROR: scanpy {sp} >= 1.11 requires numpy>=2'
import scanpy  # import must succeed
print(f'numpy {v} | anndata {ad} | scanpy {sp} OK')
"

# ── Smoke test ────────────────────────────────────────────────────────────────
sptxinsight --help
if [ "${DO_MCP}" -eq 1 ]; then
    sptxinsight-mcp --help
fi
