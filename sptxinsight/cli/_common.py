"""Shared helpers for sptxinsight CLI subcommands."""

from __future__ import annotations

import re
from typing import List

from ..uri_path import _LIST_SCHEMES
from ..uri_path import URIPath
from ._paths import default_storage_kwargs

_STORAGE_KWARGS = default_storage_kwargs()

_SAMPLE_SUFFIXES = (".h5ad", ".zarr")


def csv_to_list(ctx, param, value):  # noqa: ANN001 - Click callback signature
    """Click callback: split a comma/space separated string into a list."""
    if value is None:
        return None
    tokens = [x for x in re.split(r"[,\s]+", str(value).strip()) if x]
    return tokens or None


def enumerate_sample_uris(sptx_dir: URIPath) -> List[URIPath]:
    """List spatial-sample URIs under ``sptx_dir``.

    Accepts a directory of ``.h5ad`` / ``.zarr`` samples or a
    ``sptx-list:///path/to/list.txt`` manifest (one sample path per line, with
    an optional TAB/comma-separated 2nd column giving an explicit sample id;
    ``image-list://`` is accepted as a deprecated alias).
    """
    sptx_dir = sptx_dir.coerce_sample_list()
    if not sptx_dir.exists():
        raise FileNotFoundError(f"Spatial sample directory not found: {sptx_dir}")

    samples = sorted(
        p
        for p in sptx_dir.iterdir()
        if sptx_dir.scheme in _LIST_SCHEMES or p.suffix.lower() in _SAMPLE_SUFFIXES
    )
    return samples
