"""Cell-cell interaction (CCI) scoring for sptxinsight.

Given a set of ligand-receptor (LR) pairs and an already-ingested results
directory (``model-outputs-csv/<id>.csv`` with ``expr_*`` columns), this module
computes, for every cell, a per-pair signalling score in both directions:

* ``cci_<LIG>_<REC>_out_*`` -- the cell as a **sender**: its own ligand level
  times the (weighted) receptor level of its spatial neighbours.
* ``cci_<LIG>_<REC>_in_*``  -- the cell as a **receiver**: its own receptor level
  times the (weighted) ligand level of its spatial neighbours.

Each direction is aggregated over neighbours two ways: ``_mean`` (density-
invariant average over the neighbourhood, the recommended default for distance/
border H-Plots) and ``_sum`` (total signalling load, scales with neighbour
count).

Neighbours are the **1-hop** cells of the cached Delaunay graph, pruned to a
physical distance cutoff ``d_max`` (microns) and weighted by a distance-decay
kernel (exponential ``exp(-d/lambda)`` by default). This is contact / short-
range signalling, so no k-hop expansion is used; range is controlled purely by
``d_max`` / ``lambda``.

Outputs (under ``<results-dir>/``)::

    cci-outputs-csv/cells/<id>.csv   per-cell scores (center_x, center_y,
                                     prob_*, cci_<pair>_{in,out}_{mean,sum})
    cci-outputs.csv                  per-sample x per-pair summary

The per-cell ``cci_*`` columns are continuous, like ``expr_*``, so they can be
fed to distance/niche analyses downstream.
"""

from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import as_completed
from pathlib import Path
from typing import List
from typing import Mapping
from typing import Optional
from typing import Sequence
from typing import Tuple

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse import diags

from ..uri_path import URIPath
from .graph_cache import get_or_build_delaunay
from .insight_helpers import compute_cell_center_points

_logger = logging.getLogger(__name__)

# Path to the bundled human ligand-receptor reference (CellTalkDB-style).
_BUNDLED_LR_CSV = Path(__file__).resolve().parents[1] / "human_lr_pair.csv"

_KERNELS = ("exponential", "gaussian", "binary")


# ---------------------------------------------------------------------------
# LR pair resolution
# ---------------------------------------------------------------------------


def load_lr_pairs(
    lr_pairs_path: Optional[str | Path | URIPath] = None,
) -> List[Tuple[str, str]]:
    """Return a list of ``(ligand, receptor)`` symbol tuples.

    Reads ``lr_pairs_path`` when given (a CSV/TSV with ligand/receptor columns),
    otherwise the bundled ``human_lr_pair.csv``. Recognised column names are
    ``ligand_gene_symbol``/``receptor_gene_symbol`` or ``ligand``/``receptor``.
    """
    if lr_pairs_path is None:
        path = _BUNDLED_LR_CSV
        df = pd.read_csv(path)
    else:
        path = lr_pairs_path
        opener = getattr(path, "open", None)
        sep = "\t" if str(path).lower().endswith((".tsv", ".txt")) else ","
        if callable(opener) and not isinstance(path, (str, Path)):
            with path.open("r", encoding="utf-8") as fp:
                df = pd.read_csv(fp, sep=sep)
        else:
            df = pd.read_csv(path, sep=sep)

    lower = {c.lower(): c for c in df.columns}
    lig_col = lower.get("ligand_gene_symbol") or lower.get("ligand")
    rec_col = lower.get("receptor_gene_symbol") or lower.get("receptor")
    if lig_col is None or rec_col is None:
        raise ValueError(
            f"LR table {path!s} must have ligand/receptor columns "
            "(ligand_gene_symbol/receptor_gene_symbol or ligand/receptor); "
            f"got {list(df.columns)}."
        )

    pairs: List[Tuple[str, str]] = []
    seen: set[Tuple[str, str]] = set()
    for lig, rec in zip(df[lig_col].astype(str), df[rec_col].astype(str), strict=False):
        lig, rec = lig.strip(), rec.strip()
        if not lig or not rec or lig.lower() == "nan" or rec.lower() == "nan":
            continue
        key = (lig.upper(), rec.upper())
        if key in seen:
            continue
        seen.add(key)
        pairs.append((lig, rec))
    return pairs


def filter_pairs_to_panel(
    pairs: Sequence[Tuple[str, str]],
    panel_genes: Sequence[str],
    restrict_genes: Optional[Sequence[str]] = None,
) -> List[Tuple[str, str]]:
    """Keep only pairs whose ligand AND receptor are in ``panel_genes``.

    Matching is case-insensitive. ``restrict_genes`` (when given) further limits
    pairs to those whose ligand and receptor are both in that list.
    """
    panel = {g.lower() for g in panel_genes}
    restrict = {g.lower() for g in restrict_genes} if restrict_genes else None
    kept: List[Tuple[str, str]] = []
    for lig, rec in pairs:
        ll, rl = lig.lower(), rec.lower()
        if ll not in panel or rl not in panel:
            continue
        if restrict is not None and (ll not in restrict or rl not in restrict):
            continue
        kept.append((lig, rec))
    return kept


# ---------------------------------------------------------------------------
# Weighted neighbour graph
# ---------------------------------------------------------------------------


def _decay_weights(length_px: np.ndarray, kernel: str, lam: float) -> np.ndarray:
    """Distance-decay edge weights for a 1-hop neighbour graph."""
    d = np.asarray(length_px, dtype=np.float64)
    if kernel == "binary":
        return np.ones_like(d)
    if kernel == "gaussian":
        return np.exp(-(d * d) / (2.0 * lam * lam))
    # exponential (default)
    return np.exp(-d / lam)


def build_weight_matrices(
    edges_df: pd.DataFrame,
    n_cells: int,
    kernel: str,
    lam: float,
):
    """Build symmetric neighbour-weight matrices ``(W_sum, W_mean)``.

    ``W_sum`` carries the raw distance-decay weights; ``W_mean`` is its row-
    normalised counterpart (each row sums to 1 where the cell has neighbours).
    Self-loops are absent (Delaunay edges never connect a cell to itself).
    """
    src = edges_df["source"].to_numpy()
    dst = edges_df["target"].to_numpy()
    w = _decay_weights(edges_df["length"].to_numpy(), kernel, lam)

    # Symmetrise: each undirected edge contributes (i->j) and (j->i).
    rows = np.concatenate([src, dst])
    cols = np.concatenate([dst, src])
    data = np.concatenate([w, w])
    W_sum = coo_matrix((data, (rows, cols)), shape=(n_cells, n_cells)).tocsr()

    rowsum = np.asarray(W_sum.sum(axis=1)).ravel()
    inv = np.divide(1.0, rowsum, out=np.zeros_like(rowsum), where=rowsum > 0)
    W_mean = diags(inv) @ W_sum
    return W_sum, W_mean


# ---------------------------------------------------------------------------
# Per-sample scoring
# ---------------------------------------------------------------------------


def _sanitize(sym: str) -> str:
    """Make a gene symbol safe for a column name (keep alnum, collapse rest)."""
    return "".join(ch if ch.isalnum() else "" for ch in sym).upper()


def score_sample(
    df: pd.DataFrame,
    pairs: Sequence[Tuple[str, str]],
    *,
    slide_id: str,
    mpp_um_per_px: float,
    d_max_um: float,
    kernel: str,
    lam_um: float,
    graph_cache_dir: Optional[Path] = None,
) -> pd.DataFrame:
    """Compute per-cell CCI scores for one sample.

    Returns a DataFrame aligned to ``df`` rows with ``center_x``/``center_y``,
    any ``prob_*`` columns (copied through for downstream celltype anchoring),
    and four ``cci_<pair>_{in,out}_{mean,sum}`` columns per pair.
    """
    df = compute_cell_center_points(df.copy())
    n = len(df)
    _cx = "center_x_um" if "center_x_um" in df.columns else "center_x"
    _cy = "center_y_um" if "center_y_um" in df.columns else "center_y"
    centers = df[[_cx, _cy]].to_numpy(dtype=np.float32)

    max_edge_len_px = float(d_max_um) / float(mpp_um_per_px)
    lam_px = float(lam_um) / float(mpp_um_per_px)
    if graph_cache_dir is not None:
        centers_int = np.asarray(centers, dtype=np.int32)
        edges_df = get_or_build_delaunay(
            graph_cache_dir, slide_id, centers_int, mpp_um_per_px, max_edge_len_px
        )
    else:
        from .insight_helpers import delaunay_triangulation

        edges_df = delaunay_triangulation(centers, max_edge_len_px)

    W_sum, W_mean = build_weight_matrices(edges_df, n, kernel, lam_px)

    # Resolve expression columns case-insensitively.
    expr_lower = {
        c[len("expr_") :].lower(): c for c in df.columns if c.startswith("expr_")
    }
    used_genes = sorted({g.lower() for pr in pairs for g in pr})
    used_genes = [g for g in used_genes if g in expr_lower]
    gidx = {g: i for i, g in enumerate(used_genes)}
    if used_genes:
        E = df[[expr_lower[g] for g in used_genes]].to_numpy(dtype=np.float64)
        E = np.nan_to_num(E)
        NS = W_sum @ E  # neighbour weighted sum, per gene  (n, G)
        NM = W_mean @ E  # neighbour weighted mean, per gene (n, G)
    else:
        E = NS = NM = np.zeros((n, 0))

    out = {"center_x_um": df[_cx].to_numpy(), "center_y_um": df[_cy].to_numpy()}
    for c in df.columns:
        if c.startswith("prob_"):
            out[c] = df[c].to_numpy()

    for lig, rec in pairs:
        ll, rl = lig.lower(), rec.lower()
        if ll not in gidx or rl not in gidx:
            continue
        li, ri = gidx[ll], gidx[rl]
        L = E[:, li]
        R = E[:, ri]
        tag = f"cci_{_sanitize(lig)}_{_sanitize(rec)}"
        out[f"{tag}_out_mean"] = L * NM[:, ri]
        out[f"{tag}_out_sum"] = L * NS[:, ri]
        out[f"{tag}_in_mean"] = R * NM[:, li]
        out[f"{tag}_in_sum"] = R * NS[:, li]

    return pd.DataFrame(out)


# ---------------------------------------------------------------------------
# Cohort driver
# ---------------------------------------------------------------------------


def cci_generation(
    results_dir: URIPath,
    *,
    lr_pairs_path: Optional[str | Path | URIPath] = None,
    restrict_genes: Optional[Sequence[str]] = None,
    d_max_um: float = 25.0,
    kernel: str = "exponential",
    lam_um: float = 25.0,
    slide_mpp_lookup: Optional[Mapping[str, float]] = None,
    num_workers: int = 4,
    overwrite: bool = False,
) -> List[str]:
    """Compute per-cell CCI scores for every ingested sample.

    Reads ``<results-dir>/model-outputs-csv/<id>.csv`` and writes
    ``cci-outputs-csv/cells/<id>.csv`` per sample plus an aggregated
    ``cci-outputs.csv``. Returns the list of sample ids that failed.
    """
    if kernel not in _KERNELS:
        raise ValueError(f"kernel must be one of {_KERNELS}, got {kernel!r}.")

    model_dir = results_dir / "model-outputs-csv"
    if not model_dir.exists():
        raise FileNotFoundError(
            f"{model_dir} not found; run `sptxinsight ingest` first."
        )
    stems = sorted(p.stem for p in model_dir.iterdir() if p.suffix.lower() == ".csv")
    if not stems:
        raise FileNotFoundError(f"No model-output CSVs under {model_dir}.")

    # Resolve pairs against the panel of the first sample (shared across cohort).
    with (model_dir / f"{stems[0]}.csv").open("r", encoding="utf-8") as fp:
        head = pd.read_csv(fp, nrows=1)
    panel_genes = [c[len("expr_") :] for c in head.columns if c.startswith("expr_")]
    if not panel_genes:
        raise ValueError(
            f"Sample {stems[0]!r} has no expr_ columns; CCI needs transcript data."
        )
    all_pairs = load_lr_pairs(lr_pairs_path)
    pairs = filter_pairs_to_panel(all_pairs, panel_genes, restrict_genes)
    if not pairs:
        raise ValueError(
            "No ligand-receptor pairs have both genes in the sample panel "
            f"({len(panel_genes)} genes). Provide --lr-pairs/--genes that match."
        )
    _logger.info(
        "CCI: %d LR pair(s) within the %d-gene panel.", len(pairs), len(panel_genes)
    )

    out_dir = results_dir / "cci-outputs-csv"
    out_dir.mkdir(exist_ok=True)
    cells_dir = out_dir / "cells"
    cells_dir.mkdir(exist_ok=True)
    graph_cache_dir = Path(str(results_dir)) / "graphs"
    graph_cache_dir.mkdir(parents=True, exist_ok=True)

    def _process(stem: str):
        cells_csv = cells_dir / f"{stem}.csv"
        if not overwrite and cells_csv.exists():
            with cells_csv.open("r", encoding="utf-8") as fp:
                return stem, pd.read_csv(fp)
        with (model_dir / f"{stem}.csv").open("r", encoding="utf-8") as fp:
            df = pd.read_csv(fp)
        mpp = 1.0
        if slide_mpp_lookup:
            mpp = slide_mpp_lookup.get(stem) or 1.0
        scored = score_sample(
            df,
            pairs,
            slide_id=stem,
            mpp_um_per_px=mpp,
            d_max_um=d_max_um,
            kernel=kernel,
            lam_um=lam_um,
            graph_cache_dir=graph_cache_dir,
        )
        with cells_csv.open("w", encoding="utf-8") as fp:
            scored.to_csv(fp, index=False)
        return stem, scored

    failed: List[str] = []
    summaries: List[pd.DataFrame] = []

    def _summarise(stem: str, scored: pd.DataFrame) -> pd.DataFrame:
        score_cols = [c for c in scored.columns if c.startswith("cci_")]
        rows = []
        for col in score_cols:
            # col = cci_<LIG>_<REC>_<dir>_<agg>
            body, agg = col.rsplit("_", 1)
            body, direction = body.rsplit("_", 1)
            _, lig, rec = body.split("_", 2)
            v = scored[col].to_numpy()
            rows.append(
                {
                    "sample": stem,
                    "ligand": lig,
                    "receptor": rec,
                    "direction": direction,
                    "aggregation": agg,
                    "n_cells": len(v),
                    "mean": float(np.mean(v)),
                    "frac_pos": float(np.mean(v > 0)),
                }
            )
        return pd.DataFrame(rows)

    if num_workers and num_workers > 1:
        with ThreadPoolExecutor(max_workers=num_workers) as ex:
            futs = {ex.submit(_process, s): s for s in stems}
            for fut in as_completed(futs):
                stem = futs[fut]
                try:
                    sid, scored = fut.result()
                    summaries.append(_summarise(sid, scored))
                except Exception as exc:  # noqa: BLE001
                    _logger.warning("CCI failed for %s: %s", stem, exc)
                    failed.append(stem)
    else:
        for stem in stems:
            try:
                sid, scored = _process(stem)
                summaries.append(_summarise(sid, scored))
            except Exception as exc:  # noqa: BLE001
                _logger.warning("CCI failed for %s: %s", stem, exc)
                failed.append(stem)

    if summaries:
        agg = pd.concat(summaries, ignore_index=True)
        with (results_dir / "cci-outputs.csv").open("w", encoding="utf-8") as fp:
            agg.to_csv(fp, index=False)

    return failed
