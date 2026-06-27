#!/usr/bin/env bash
# conda-setup.sh — create and populate the standalone sptxinsight conda environment.
#
# Usage:  sh ./conda-setup.sh          (from /workspace/wsinsight/sptxinsight/)
#
# NOTE: sptxinsight is also co-installable inside the shared wsinsight env
# via:  conda activate wsinsight && pip install --no-deps -e .
# This script creates a *separate* sptxinsight environment instead.
#
# NOTE: scanpy / spatialdata extras require numpy>=2 and are INCOMPATIBLE with
# this environment (pinned numpy<2 by pyproject.toml). Install those extras
# only in a dedicated numpy-2 environment.

set -e   # abort on first error

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# ── Reset environment ─────────────────────────────────────────────────────────
source /opt/anaconda3/etc/profile.d/conda.sh
conda deactivate
conda env remove -n sptxinsight -y 2>/dev/null || true

# Python only — sptxinsight uses geopandas via pip+pyogrio (no GDAL binary needed).
conda create -n sptxinsight python=3.11 "setuptools<67" -c conda-forge -y
conda activate sptxinsight
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
pip install anndata "zarr<3"          # zarr<3 keeps anndata compat with wsinsight
pip install scikit-learn joblib

# ── Geometry / GIS — pyogrio as OGR backend (no GDAL binary required) ─────────
pip install pyogrio shapely geopandas

# ── Graph clustering ──────────────────────────────────────────────────────────
pip install python-igraph leidenalg

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
print(f'numpy {v} OK')
"

# ── Smoke test ────────────────────────────────────────────────────────────────
sptxinsight --help
sptxinsight-mcp --help
