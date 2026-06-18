"""Adapt an AnnData spatial sample into the WSInsight H-Plot CSV contract.

The vendored :func:`sptxinsight.insightlib.hplot_generation.hplot_generation`
reads one CSV per sample from ``<results_dir>/model-outputs-csv/<slide_id>.csv``
with:

- ``center_x`` / ``center_y``: cell-center coordinates **in microns**,
- ``prob_<type>`` columns: per-cell soft membership; the engine takes ``idxmax``
  to get the hard label at runtime, and
- ``expr_<gene>`` columns (optional): per-cell expression of selected genes,
  written only when a gene-based base/target is requested. The engine averages
  them per layer to build expression H-Plots.

Xenium/Visium coordinates in ``adata.obsm["spatial"]`` are already microns, so
the matching ``slide_mpp_lookup`` value is ``1.0`` (see :mod:`sptxinsight.pipeline`).
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING, List, Mapping, Sequence

import numpy as np
import pandas as pd

from .uri_path import URIPath

if TYPE_CHECKING:  # pragma: no cover - typing only
    from anndata import AnnData

_MODEL_OUTPUTS_SUBDIR = "model-outputs-csv"


def _sanitize_type(name: str) -> str:
    """Normalize a cell-type or gene label into a ``prob_``/``expr_`` suffix."""
    slug = re.sub(r"\s+", "_", str(name).strip().lower())
    slug = re.sub(r"[^0-9a-z_]+", "_", slug)
    return re.sub(r"_+", "_", slug).strip("_")


def _resolve_expression_matrix(adata: "AnnData", which: str):
    """Return ``(matrix, var_names)`` for the requested expression source.

    ``which`` accepts ``"X"`` (default), ``"raw"`` / ``"raw.X"`` for
    ``adata.raw.X``, or the name of an entry in ``adata.layers``. A leading
    ``"adata."`` is tolerated (e.g. ``"adata.raw.X"``).
    """
    token = (which or "X").strip()
    if token.lower().startswith("adata."):
        token = token[len("adata."):]
    low = token.lower()
    if low in ("", "x"):
        return adata.X, list(adata.var_names)
    if low in ("raw", "raw.x"):
        if adata.raw is None:
            raise ValueError(
                "adata.raw is None but --expression-matrix raw was requested."
            )
        return adata.raw.X, list(adata.raw.var_names)
    if token in adata.layers:
        return adata.layers[token], list(adata.var_names)
    raise ValueError(
        f"Unknown expression matrix {token!r}; expected 'X', 'raw', or a layer "
        f"in {list(adata.layers)}."
    )


def _extract_gene_columns(
    adata: "AnnData", genes: Sequence[str], expression_matrix: str
) -> Mapping[str, np.ndarray]:
    """Extract dense per-cell vectors for ``genes`` from the chosen matrix.

    Gene symbols are matched exactly first, then case-insensitively. Raises
    ``KeyError`` listing any symbols absent from the matrix.
    """
    from scipy import sparse

    mat, var_names = _resolve_expression_matrix(adata, expression_matrix)
    name_to_idx = {n: i for i, n in enumerate(var_names)}
    lower_to_idx = {str(n).lower(): i for i, n in enumerate(var_names)}

    out: dict[str, np.ndarray] = {}
    missing: list[str] = []
    for gene in genes:
        idx = name_to_idx.get(gene)
        if idx is None:
            idx = lower_to_idx.get(str(gene).lower())
        if idx is None:
            missing.append(gene)
            continue
        col = mat[:, idx]
        if sparse.issparse(col):
            col = np.asarray(col.todense()).ravel()
        else:
            col = np.asarray(col).ravel()
        out[gene] = col.astype(float, copy=False)
    if missing:
        raise KeyError(
            f"Gene(s) not found in expression matrix '{expression_matrix}': "
            f"{missing}"
        )
    return out


def anndata_to_contract(
    adata: "AnnData",
    slide_id: str,
    results_dir: "str | URIPath",
    *,
    cell_type_key: str = "cell_type",
    spatial_key: str = "spatial",
    genes: Sequence[str] = (),
    expression_matrix: str = "X",
) -> List[str]:
    """Write ``model-outputs-csv/<slide_id>.csv`` for one AnnData sample.

    Parameters
    ----------
    adata:
        Annotated sample. Requires ``obsm[spatial_key]`` (N x 2 micron coords)
        and a categorical/string ``obs[cell_type_key]``.
    slide_id:
        Stem used for the output CSV and for ``slide_mpp_lookup`` keys.
    results_dir:
        H-Plot results directory; the CSV is written under its
        ``model-outputs-csv/`` subdirectory (created if absent).
    cell_type_key:
        Column in ``adata.obs`` holding the per-cell type label.
    spatial_key:
        Key in ``adata.obsm`` holding the spatial coordinates (microns).
    genes:
        Optional gene symbols to also write as ``expr_<gene>`` columns (used
        for gene-based base/target H-Plots). Empty by default.
    expression_matrix:
        Source matrix for ``genes`` extraction: ``"X"`` (default), ``"raw"``,
        or a layer name in ``adata.layers``.

    Returns
    -------
    list[str]
        The sanitized cell-type vocabulary (without the ``prob_`` prefix),
        suitable for ``--base-type`` / ``--target-type`` resolution.
    """
    if spatial_key not in adata.obsm:
        raise KeyError(
            f"adata.obsm[{spatial_key!r}] not found; cannot build H-Plot contract."
        )
    if cell_type_key not in adata.obs:
        raise KeyError(
            f"adata.obs[{cell_type_key!r}] not found; sample must be cell-typed."
        )

    xy = np.asarray(adata.obsm[spatial_key], dtype=float)
    if xy.ndim != 2 or xy.shape[1] < 2:
        raise ValueError(
            f"adata.obsm[{spatial_key!r}] must be N x 2 (got shape {xy.shape})."
        )

    df = pd.DataFrame({"center_x": xy[:, 0], "center_y": xy[:, 1]})

    labels = adata.obs[cell_type_key].astype(str).to_numpy()
    onehot = pd.get_dummies(labels).astype(float)
    onehot.columns = [f"prob_{_sanitize_type(c)}" for c in onehot.columns]
    onehot.index = df.index

    frames = [df, onehot]
    if genes:
        gene_cols = _extract_gene_columns(adata, genes, expression_matrix)
        expr = pd.DataFrame(
            {f"expr_{_sanitize_type(g)}": v for g, v in gene_cols.items()}
        )
        expr.index = df.index
        frames.append(expr)

    out = pd.concat(frames, axis=1)

    results_dir = (
        results_dir if isinstance(results_dir, URIPath) else URIPath(str(results_dir))
    )
    out_dir = results_dir / _MODEL_OUTPUTS_SUBDIR
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / f"{slide_id}.csv"
    with csv_path.open("w") as fp:
        out.to_csv(fp, index=False)

    return sorted(c[len("prob_"):] for c in onehot.columns)
