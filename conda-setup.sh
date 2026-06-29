#!/usr/bin/env bash
# conda-setup.sh — create and populate the standalone sptxinsight conda environment.
#
# Usage:  sh ./conda-setup.sh [-n ENV_NAME] [-r|--reset]
#
#   -n | --name  ENV_NAME   Conda environment to use (default: current active env).
#   -r | --reset            Deactivate, remove, recreate, and activate the env.
#                           Without this flag the script skips env creation and
#                           only (re-)installs packages into the existing env.
#
# NOTE: sptxinsight is also co-installable inside the shared wsinsight env
# via:  conda activate wsinsight && pip install --no-deps -e .
# This script creates a *separate* sptxinsight environment instead.
#
# NOTE: spatialdata / squidpy require numpy>=2 and are INCOMPATIBLE with this
# environment (pinned numpy<2 by pyproject.toml). scanpy/anndata ARE installed
# here, but pinned to the last numpy<2-compatible line: anndata>=0.12,<0.13 and
# scanpy<1.11. anndata 0.12 is required to read 0.12-format `annotated.h5ad`
# (older anndata raises IORegistryError on encoding_type='null'); newer scanpy
# (>=1.11) drags in numpy>=2. Do NOT relax these pins.

set -e   # abort on first error

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# ── Argument parsing ──────────────────────────────────────────────────────────
ENV_NAME="${CONDA_DEFAULT_ENV:-}"   # default = current active env
DO_RESET=0

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
        *)
            echo "Unknown option: $1" >&2
            echo "Usage: sh ./conda-setup.sh [-n ENV_NAME] [-r|--reset]" >&2
            exit 1
            ;;
    esac
done

if [[ -z "$ENV_NAME" ]]; then
    echo "Error: no conda environment specified and no environment is currently active." >&2
    echo "       Use -n ENV_NAME to specify one." >&2
    exit 1
fi

echo "Target conda environment: ${ENV_NAME}  (reset=${DO_RESET})"

# ── (Re-)create environment ───────────────────────────────────────────────────
source /opt/anaconda3/etc/profile.d/conda.sh

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

# ── Torch stack (heaviest; download first) ────────────────────────────────────
pip install torch torchvision torch-geometric

# ── Core scientific / bioinformatics stack ────────────────────────────────────
pip install scipy pandas h5py tqdm click
# Pin anndata/scanpy to the last numpy<2-compatible line. anndata>=0.12 reads
# the 0.12-format h5ad written by wsinsight; scanpy<1.11 avoids the numpy>=2
# requirement. zarr<3 keeps anndata zarr-readers compatible with wsinsight.
pip install "anndata>=0.12,<0.13" "scanpy<1.11" "zarr<3"
pip install scikit-learn joblib

# ── Geometry / GIS — pyogrio as OGR backend (no GDAL binary required) ─────────
pip install pyogrio shapely geopandas

# ── Graph clustering ──────────────────────────────────────────────────────────
pip install igraph leidenalg

# ── Cloud I/O (version-capped to stay compatible with wsinsight's zarr<3 stack)
# Pre-install aiobotocore + boto3 with explicit compatible versions to avoid
# pip spending minutes backtracking through 90+ boto3 versions.
pip install "aiobotocore>=2.5.4,<3.0.0" "boto3>=1.41,<1.42"
pip install "fsspec>=2023.1.0,<2026" "s3fs<2026" "gcsfs<2026" \
    requests platformdirs

# ── Optional extras (always install — both CLI entry points should work) ──────
pip install "fastmcp>=2.0"       # sptxinsight-mcp server
pip install "harmonypy>=0.0.9"   # --cme-batch-correct harmony

# ── Install sptxinsight itself ────────────────────────────────────────────────
pip install -e "${SCRIPT_DIR}"

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
sptxinsight-mcp --help
