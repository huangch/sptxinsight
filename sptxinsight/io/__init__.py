"""Spatial-sample I/O for sptxinsight.

Mirrors WSInsight's ``wsi`` package role (``set_backend`` + readers) but for
spatial-omics objects. The default backend reads AnnData ``.h5ad`` files.
"""

from __future__ import annotations

from typing import Literal

from .anndata_reader import read_sample

_BACKEND: str = "anndata"
_VALID_BACKENDS = ("anndata", "zarr", "spatialdata")


def set_backend(name: Literal["anndata", "zarr", "spatialdata"]) -> None:
    """Select the reader backend used by :func:`read_sample`."""
    global _BACKEND
    if name not in _VALID_BACKENDS:
        raise ValueError(f"Unknown backend {name!r}; choose from {_VALID_BACKENDS}.")
    _BACKEND = name


def get_backend() -> str:
    """Return the currently selected reader backend."""
    return _BACKEND


__all__ = ["set_backend", "get_backend", "read_sample"]
