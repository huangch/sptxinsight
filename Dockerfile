# ====================================
# sptxinsight — standalone GPU image
# CUDA 12.8 + cuDNN + Ubuntu 22.04
# ====================================
# Independent of the wsinsight image. Ships the sptxinsight CLI and the
# `sptxinsight-mcp` MCP server. Runs cell-typing / niche niche discovery / H-Plot
# / CCI over spatial-transcriptomics data on the GPU.
FROM nvidia/cuda:12.8.0-cudnn-devel-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

# ------------------------------------
# Basic system dependencies
# (libgdal for geopandas/shapely region merge; the rest are common wheels deps)
# ------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    git curl wget unzip vim ca-certificates build-essential pkg-config \
    libglib2.0-0 && \
    rm -rf /var/lib/apt/lists/*

# ------------------------------------
# Install Miniconda (Python 3.11 base) — gdal from conda-forge so geopandas
# resolves cleanly against numpy<2.
# ------------------------------------
ENV CONDA_DIR=/opt/conda
RUN curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/mc.sh \
 && bash /tmp/mc.sh -b -p "$CONDA_DIR" \
 && rm /tmp/mc.sh
ENV PATH="$CONDA_DIR/bin:$PATH"

RUN conda --version && \
    (conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main -y || true) && \
    (conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r -y || true)

RUN conda update -n base --yes --override-channels -c conda-forge conda && \
    conda create -y --override-channels -n sptxinsight -c conda-forge \
        python=3.11 gdal=3.11.3 pip "setuptools<67" && \
    conda clean -afy

# ------------------------------------
# Global Conda initialization (Docker bash doesn't read /etc/profile.d/*)
# ------------------------------------
RUN echo '. /opt/conda/etc/profile.d/conda.sh' >> /etc/bash.bashrc && \
    echo 'conda activate sptxinsight' >> /etc/bash.bashrc && \
    echo '. /opt/conda/etc/profile.d/conda.sh' >> /etc/skel/.bashrc && \
    echo 'conda activate sptxinsight' >> /etc/skel/.bashrc

ENV CONDA_DEFAULT_ENV=sptxinsight
ENV PATH="$CONDA_DIR/envs/sptxinsight/bin:$PATH"

RUN python -m pip install --upgrade pip

# ------------------------------------
# GPU stack: torch + torch_geometric (niche niche autoencoder training)
# Copied ahead of the source tree so this layer is constrained (unpinned, torch
# resolves to >=2.9 and pulls the CUDA 13 wheels over nvidia/cudnn/lib) and so
# the cache survives source edits.
# ------------------------------------
COPY constraints.txt /app/sptxinsight/constraints.txt
RUN pip install --retries 10 -c /app/sptxinsight/constraints.txt \
        "numpy<2" torch torchvision torch_geometric

# ------------------------------------
# Install sptxinsight with MCP server, zarr (zarr<3) and harmony extras.
# NOTE: scanpy is a core sptxinsight dependency but versions ≥1.11 require
# numpy>=2, which conflicts with the numpy<2 stack here.  Pin scanpy to a
# numpy<2-compatible release and pin numpy<2 after to prevent pip from
# upgrading it during dependency resolution.
# ------------------------------------
WORKDIR /app/sptxinsight
COPY . .
RUN pip install --retries 10 -c /app/sptxinsight/constraints.txt \
        "numpy<2" ".[mcp,zarr,harmony]" && \
    pip install --retries 10 "numpy<2"

# ------------------------------------
# Sanity check (build-time)
# ------------------------------------
RUN python - <<'PY'
import torch
import sptxinsight
from sptxinsight.mcp.server import build_server
m = build_server(experimental=True)
print("sptxinsight:", getattr(sptxinsight, "__version__", "?"),
      "Torch:", torch.__version__, "CUDA:", torch.version.cuda,
      "GPU?", torch.cuda.is_available(), "MCP:", type(m).__name__)
PY

# Fail the build if the console scripts were not installed. Without this the
# image can ship with an importable package but no `sptxinsight` on PATH.
RUN command -v sptxinsight && \
    command -v sptxinsight-mcp && \
    sptxinsight --help > /dev/null

# ------------------------------------
# Non-root user
# ------------------------------------
# The container starts as root; the runtime entrypoint (docker-entrypoint.sh)
# remaps the pre-created ``user`` account to the owner of the mounted /workspace
# (or to $HOST_UID/$HOST_GID) and then drops privileges via setpriv. The uid/gid
# baked here is only a throwaway placeholder, immediately overwritten at run
# time, so it is hard-coded (1000) rather than exposed as a build ARG — the
# shared entrypoint looks the account up by the name ``user``.
RUN groupadd -g 1000 user && \
    useradd -m -u 1000 -g 1000 -s /bin/bash user && \
    bash -lc 'echo ". /opt/conda/etc/profile.d/conda.sh" >> /home/user/.bashrc' && \
    bash -lc 'echo "conda activate sptxinsight" >> /home/user/.bashrc' && \
    chown -R 1000:1000 /home/user

# Install the runtime uid/gid-remapping entrypoint.
RUN install -m 0755 /app/sptxinsight/docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh

WORKDIR /workspace
RUN chown -R 1000:1000 /workspace
# NOTE: no ``USER`` here on purpose — the container starts as root so the
# entrypoint can remap ``user`` to the mount owner, then drops privileges via
# setpriv. Passing ``docker run --user ...`` still works: the entrypoint detects
# a non-root start and execs the command unchanged.

SHELL ["/bin/bash","-lc"]
ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]
CMD ["bash"]
