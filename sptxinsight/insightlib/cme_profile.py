"""Profile discovered CMEs (niches) so each ``cme_n`` can be named.

Reads ``cme-outputs-csv/cells/<id>.csv`` produced by ``sptxinsight cme`` and, for
each CME cluster, summarises:

* the mean cell-type composition (from ``prob_*`` columns), and
* the most enriched marker genes (from ``expr_*`` columns, if present).

These two fingerprints let a human assign biological names (TLS, tumor nest,
perivascular, ...) to the otherwise arbitrary cluster ids.

This module deliberately depends only on numpy/pandas so it can be imported and
run without the deep-learning ``[cme]`` extra.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd


# Mode -> (cells subdir, one-hot column prefix, output filename suffix). Kept in
# sync with ``_CME_MODE_SPEC`` in cme_generation.py but defined locally so this
# module stays free of the heavy deep-learning ``[cme]`` extra.
_MODE_SPEC = {
    "celltype":   ("cme-outputs-csv",        "cme_",    ""),
    "expression": ("cme-gex-outputs-csv",    "gexcme_", "-gex"),
    "both":       ("cme-hybrid-outputs-csv", "hcme_",   "-hybrid"),
}


def _cme_columns(columns: List[str], prefix: str = "cme_") -> List[str]:
    """Return ``<prefix>*`` columns ordered by their integer index."""
    cols = [c for c in columns if c.startswith(prefix)]

    def _idx(c: str) -> int:
        try:
            return int(c[len(prefix):])
        except (IndexError, ValueError):
            return 1 << 30

    return sorted(cols, key=_idx)


def cme_profile(
    results_dir: str | Path,
    top_genes: int = 10,
    top_types: int = 5,
    write: bool = True,
    mode: str = "celltype",
) -> Tuple[pd.DataFrame, Optional[pd.DataFrame]]:
    """Aggregate per-cell CME labels into per-cluster composition + marker tables.

    ``mode`` selects which niche family to profile (``celltype`` ->
    ``cme-outputs-csv``/``cme_``; ``expression`` -> ``cme-gex-outputs-csv``/
    ``gexcme_``; ``both`` -> ``cme-hybrid-outputs-csv``/``hcme_``).

    Returns ``(composition_df, markers_df)``. ``markers_df`` is ``None`` when the
    cells CSVs carry no ``expr_`` columns. When ``write`` is True the tables are
    saved to ``<results-dir>/cme-profile-composition<suffix>.csv`` and
    ``<results-dir>/cme-profile-markers<suffix>.csv``.
    """
    if mode not in _MODE_SPEC:
        raise ValueError(f"mode must be one of {sorted(_MODE_SPEC)}, got {mode!r}")
    subdir, prefix, suffix = _MODE_SPEC[mode]

    results_dir = Path(str(results_dir))
    cells_dir = results_dir / subdir / "cells"
    csvs = sorted(cells_dir.glob("*.csv"))
    if not csvs:
        raise FileNotFoundError(
            f"No CME cell CSVs found under {cells_dir}. "
            f"Run `sptxinsight cme --cme-mode {mode}` first."
        )

    head = pd.read_csv(csvs[0], nrows=1)
    cme_cols = _cme_columns(list(head.columns), prefix)
    if not cme_cols:
        raise ValueError(
            f"No {prefix} columns in {csvs[0]}; was this produced by "
            f"`sptxinsight cme --cme-mode {mode}`?"
        )
    prob_cols = [c for c in head.columns if c.startswith("prob_")]
    expr_cols = [c for c in head.columns if c.startswith("expr_")]

    K = len(cme_cols)
    counts = np.zeros(K, dtype=np.float64)
    sum_prob = np.zeros((K, len(prob_cols)), dtype=np.float64)
    sum_expr = np.zeros((K, len(expr_cols)), dtype=np.float64) if expr_cols else None

    usecols = cme_cols + prob_cols + expr_cols
    for csv in csvs:
        df = pd.read_csv(csv, usecols=lambda c: c in usecols)
        oh = df[cme_cols].fillna(0).to_numpy()
        assigned = oh.sum(axis=1) > 0
        if not assigned.any():
            continue
        lab = oh[assigned].argmax(axis=1)
        P = df.loc[assigned, prob_cols].fillna(0).to_numpy(dtype=np.float64) if prob_cols else None
        E = df.loc[assigned, expr_cols].fillna(0).to_numpy(dtype=np.float64) if expr_cols else None
        for k in range(K):
            m = lab == k
            nk = int(m.sum())
            if nk == 0:
                continue
            counts[k] += nk
            if P is not None:
                sum_prob[k] += P[m].sum(axis=0)
            if E is not None:
                sum_expr[k] += E[m].sum(axis=0)

    denom = np.maximum(counts, 1.0)[:, None]
    total = max(counts.sum(), 1.0)

    # ---- composition table ----
    comp = pd.DataFrame(
        sum_prob / denom,
        columns=[c[len("prob_"):] for c in prob_cols],
        index=[f"{prefix}{k}" for k in range(K)],
    )
    comp.insert(0, "n_cells", counts.astype(int))
    comp.insert(1, "frac", counts / total)
    if prob_cols:
        type_names = [c[len("prob_"):] for c in prob_cols]
        comp["top_types"] = [
            ", ".join(
                f"{type_names[j]}={comp.iloc[k][type_names[j]]:.2f}"
                for j in np.argsort(-(sum_prob[k] / denom[k]))[:top_types]
            )
            for k in range(K)
        ]

    # ---- marker table ----
    markers = None
    if expr_cols and sum_expr is not None:
        mean_expr = sum_expr / denom
        global_mean = sum_expr.sum(axis=0) / total
        eps = 1e-6
        enr = np.log2((mean_expr + eps) / (global_mean[None, :] + eps))
        genes = [c[len("expr_"):] for c in expr_cols]
        rows = []
        for k in range(K):
            order = np.argsort(-enr[k])[:top_genes]
            for rank, gi in enumerate(order, start=1):
                rows.append({
                    "cme": f"{prefix}{k}",
                    "rank": rank,
                    "gene": genes[gi],
                    "mean_expr": float(mean_expr[k, gi]),
                    "log2_enrichment": float(enr[k, gi]),
                })
        markers = pd.DataFrame(rows)

    if write:
        comp.to_csv(results_dir / f"cme-profile-composition{suffix}.csv")
        if markers is not None:
            markers.to_csv(results_dir / f"cme-profile-markers{suffix}.csv", index=False)

    return comp, markers


def _nmi(a: np.ndarray, b: np.ndarray) -> float:
    """Normalized mutual information (arithmetic mean), numpy-only."""
    a = np.asarray(a)
    b = np.asarray(b)
    n = a.size
    if n == 0:
        return float("nan")
    ct = pd.crosstab(pd.Series(a), pd.Series(b)).to_numpy(dtype=np.float64)
    pij = ct / n
    pi = pij.sum(axis=1, keepdims=True)
    pj = pij.sum(axis=0, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        mi = np.nansum(pij * np.log(pij / (pi @ pj)))
        hi = -np.nansum(pi * np.log(pi))
        hj = -np.nansum(pj * np.log(pj))
    denom = (hi + hj) / 2.0
    return float(mi / denom) if denom > 0 else 0.0


def cme_agreement(
    results_dir: str | Path,
    write: bool = True,
) -> Optional[Tuple[float, pd.DataFrame]]:
    """Compare cell-type niches against gene-expression niches on the same cells.

    Requires both ``cme-outputs-csv/cells/`` (``cme_*``) and
    ``cme-gex-outputs-csv/cells/`` (``gexcme_*``). For each shared sample the two
    label sets are aligned by row index (both derive from the same
    model-output CSV, so rows correspond), restricted to cells assigned in both,
    argmax-decoded, and pooled. Returns ``(nmi, crosstab)`` or ``None`` when one
    family is missing. With ``write`` the count crosstab is saved to
    ``<results-dir>/cme-agreement.csv``.
    """
    results_dir = Path(str(results_dir))
    ct_dir = results_dir / "cme-outputs-csv" / "cells"
    gx_dir = results_dir / "cme-gex-outputs-csv" / "cells"
    if not ct_dir.exists() or not gx_dir.exists():
        return None

    ct_csvs = {p.name: p for p in ct_dir.glob("*.csv")}
    gx_csvs = {p.name: p for p in gx_dir.glob("*.csv")}
    common = sorted(set(ct_csvs) & set(gx_csvs))
    if not common:
        return None

    ct_chunks: List[np.ndarray] = []
    gx_chunks: List[np.ndarray] = []
    for name in common:
        ct = pd.read_csv(ct_csvs[name], usecols=lambda c: c.startswith("cme_"))
        gx = pd.read_csv(gx_csvs[name], usecols=lambda c: c.startswith("gexcme_"))
        ct_cols = _cme_columns(list(ct.columns), "cme_")
        gx_cols = _cme_columns(list(gx.columns), "gexcme_")
        if not ct_cols or not gx_cols:
            continue
        ct_oh = ct[ct_cols].fillna(0).to_numpy()
        gx_oh = gx[gx_cols].fillna(0).to_numpy()
        n = min(len(ct_oh), len(gx_oh))
        ct_oh, gx_oh = ct_oh[:n], gx_oh[:n]
        assigned = (ct_oh.sum(axis=1) > 0) & (gx_oh.sum(axis=1) > 0)
        if not assigned.any():
            continue
        ct_chunks.append(ct_oh[assigned].argmax(axis=1))
        gx_chunks.append(gx_oh[assigned].argmax(axis=1))

    if not ct_chunks:
        return None

    ct_lab = np.concatenate(ct_chunks)
    gx_lab = np.concatenate(gx_chunks)
    nmi = _nmi(ct_lab, gx_lab)
    crosstab = pd.crosstab(
        pd.Series(ct_lab, name="celltype_niche"),
        pd.Series(gx_lab, name="gene_niche"),
    )
    crosstab.index = [f"cme_{i}" for i in crosstab.index]
    crosstab.columns = [f"gexcme_{j}" for j in crosstab.columns]

    if write:
        crosstab.to_csv(results_dir / "cme-agreement.csv")

    return nmi, crosstab

