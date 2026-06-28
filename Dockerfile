# ====================================
# sptxinsight — standalone GPU image
# CUDA 12.8 + cuDNN + Ubuntu 22.04
# ====================================
# Independent of the wsinsight image. Ships the sptxinsight CLI and the
# `sptxinsight-mcp` MCP server. Runs cell-typing / CME niche discovery / H-Plot
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
    git curl wget unzip ca-certificates build-essential pkg-config \
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
# GPU stack: torch + torch_geometric (CME niche autoencoder training)
# ------------------------------------
RUN pip install --retries 10 "numpy<2" torch torchvision torch_geometric

# ------------------------------------
# Install sptxinsight with MCP server, zarr (zarr<3) and harmony extras.
# NOTE: scanpy is a core sptxinsight dependency but versions ≥1.11 require
# numpy>=2, which conflicts with the numpy<2 stack here.  Pin scanpy to a
# numpy<2-compatible release and pin numpy<2 after to prevent pip from
# upgrading it during dependency resolution.
# ------------------------------------
WORKDIR /app/sptxinsight
COPY . .
RUN pip install --retries 10 "scanpy<1.11" "numpy<2" ".[mcp,zarr,harmony]" && \
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

# ------------------------------------
# Non-root user
# ------------------------------------
ARG USERNAME=user
ARG UID=1000
ARG GID=1000
RUN groupadd -g ${GID} ${USERNAME} && \
    useradd -m -u ${UID} -g ${GID} -s /bin/bash ${USERNAME} && \
    bash -lc 'echo ". /opt/conda/etc/profile.d/conda.sh" >> /home/'"${USERNAME}"'/.bashrc' && \
    bash -lc 'echo "conda activate sptxinsight" >> /home/'"${USERNAME}"'/.bashrc' && \
    chown -R ${UID}:${GID} /home/${USERNAME}
WORKDIR /workspace
RUN chown -R ${UID}:${GID} /workspace
USER ${USERNAME}

SHELL ["/bin/bash","-lc"]
CMD ["bash"]
