#!/usr/bin/env bash
# conda-setup.sh — create and populate the standalone sptxinsight conda environment.
#
# <<<USAGE_START>>>
# Usage:  sh ./conda-setup.sh ENV_NAME [-r|--reset] [-m|--mcp] [-d|--dev] [-h|--help]
#
#   ENV_NAME                (positional, REQUIRED) Conda environment to use/create.
#                           There is NO fallback to the currently-activated conda env:
#                           the name is mandatory so `-r` can never accidentally
#                           destroy a different active environment.
#   -r | --reset            Deactivate, remove, recreate, and activate ENV_NAME.
#                           Without this flag the script skips env creation and
#                           only (re-)installs packages into the existing env.
#   -m | --mcp              Also install fastmcp (MCP server support).
#                           Not installed by default to avoid entangling fastmcp's
#                           jaraco.* dep chain with the main resolution.
#   -d | --dev              Also install the [dev] extra (pytest, pytest-cov,
#                           pre_commit) so the post-install smoke test can run the
#                           real test suite. Without -d the suite is SKIPped if
#                           pytest is missing; with -d it FAILS (you asked for it).
#                           The package itself is always installed editable (-e).
#   -h | --help             Print this help message and exit.
# <<<USAGE_END>>>
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
# ENV_NAME is the FIRST POSITIONAL argument and is REQUIRED. We deliberately do
# NOT fall back to $CONDA_DEFAULT_ENV: with `-r` in play, a hidden dependency on
# whatever env happens to be active is a footgun (it would silently destroy an
# unrelated env). Make the caller name the env explicitly, every time.
DO_RESET=0
DO_MCP=0
DO_DEV=0

print_usage() {
    awk '
        /<<<USAGE_START>>>/ {capture=1; next}
        /<<<USAGE_END>>>/   {capture=0}
        capture            {sub(/^# ?/, ""); print}
    ' "$0"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            print_usage
            exit 0
            ;;
        -r|--reset)
            DO_RESET=1
            shift
            ;;
        -m|--mcp)
            DO_MCP=1
            shift
            ;;
        -d|--dev)
            DO_DEV=1
            shift
            ;;
        -*)
            echo "Unknown option: $1" >&2
            echo "Run '${0##*/} --help' for usage." >&2
            exit 1
            ;;
        *)
            # First non-flag token is ENV_NAME. Reject a second positional
            # argument (we only have one positional slot).
            if [[ -n "${ENV_NAME:-}" ]]; then
                echo "Error: only one positional argument (ENV_NAME) is accepted; got '$ENV_NAME' and '$1'." >&2
                echo "Run '${0##*/} --help' for usage." >&2
                exit 1
            fi
            ENV_NAME="$1"
            shift
            ;;
    esac
done

if [[ -z "${ENV_NAME:-}" ]]; then
    echo "Error: ENV_NAME is required." >&2
    echo "       Got: $0 $*" >&2
    echo "       Run '${0##*/} --help' for usage." >&2
    exit 1
fi

echo "Target conda environment: ${ENV_NAME}  (reset=${DO_RESET}, mcp=${DO_MCP}, dev=${DO_DEV})"

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
# aiobotocore, boto3, fsspec, s3fs, gcsfs, igraph, leidenalg, zarr, harmonypy, etc.
# With -d/--dev, also install the [dev] extra (pytest, pytest-cov, ruff, pre_commit)
# so the smoke test can run the suite; without -d, the suite is SKIPped if
# pytest is missing and only WARN-ed if it fails.
if [[ "${DO_DEV}" -eq 1 ]]; then
    pip install -c "${SCRIPT_DIR}/constraints.txt" -e "${SCRIPT_DIR}[dev]"
else
    pip install -c "${SCRIPT_DIR}/constraints.txt" -e "${SCRIPT_DIR}"
fi

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
    elif [ "${DO_DEV}" -eq 1 ]; then
        # User asked for the [dev] extra: pytest should be present. FAIL loudly
        # instead of silently SKIPping, or the install is misconfigured.
        echo "  FAIL  test suite: pytest missing but -d/--dev was requested" >&2
        smoke "pytest importable (dev)" python -c "import pytest"
    else
        echo "  SKIP  test suite (pytest not installed; rerun with -d/--dev)"
    fi
fi

if [ "${SMOKE_FAIL}" -ne 0 ]; then
    echo "smoke test: ${SMOKE_FAIL} check(s) FAILED" >&2
    exit 1
fi
echo "smoke test: all checks passed"
