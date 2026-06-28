"""Read spatial-omics samples into AnnData.

The reader materializes remote URIs (``s3://``, ``gs://`` …) to a local path via
:class:`~sptxinsight.uri_path.URIPath`, then loads them with ``scanpy``, falling
back to ``anndata`` when scanpy is unavailable. Only ``.h5ad`` is supported in
the MVP; ``.zarr`` and SpatialData are stubbed for later.
"""

from __future__ import annotations

import os
from typing import TYPE_CHECKING

from ..uri_path import URIPath

if TYPE_CHECKING:  # pragma: no cover - typing only
    from anndata import AnnData


def read_sample(uri: "str | URIPath") -> "AnnData":
    """Load a single spatial sample into an :class:`AnnData`.

    Parameters
    ----------
    uri:
        Local path or fsspec URI to an ``.h5ad`` (or ``.zarr``) sample.
    """
    import scanpy as sc

    path = uri if isinstance(uri, URIPath) else URIPath(str(uri))
    local = os.fspath(path)  # downloads + caches remote URIs; no-op for local
    suffix = str(local).lower()
    if suffix.endswith(".h5ad"):
        return sc.read_h5ad(local)
    if suffix.endswith(".zarr"):
        import anndata as ad
        return ad.read_zarr(local)
    raise ValueError(
        f"Unsupported spatial sample format: {path}. Expected .h5ad or .zarr."
    )
