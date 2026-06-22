"""Cellular microenvironment (CME) discovery for sptxinsight.

Vendored and adapted from WSInsight's ``insightlib.cme_generation``. The pipeline
builds a per-sample Delaunay cell graph, computes k-hop cell-type composition
features, trains one shared Deep Graph Infomax (DGI) encoder across all samples,
clusters the embeddings into recurring microenvironments, and writes per-cell CME
labels (and, optionally, merged annotation-level regions).

Differences from the WSInsight original:

* No whole-slide images. Spatial coordinates in ``model-outputs-csv/<id>.csv``
  are already in microns, so ``mpp = 1.0`` (resolved through ``slide_mpp_lookup``;
  there is nothing to read from a slide).
* The H-Optimus image-morphology branch is removed (sptxinsight has no pixels).
  Features are purely k-hop cell-type composition.
* Slide graphs are built sequentially so the inner k-hop process pool is not
  nested inside an outer one.

Heavy dependencies (``torch``, ``torch_geometric``, ``scikit-learn``,
``python-igraph``, ``leidenalg``) are imported at module load, so this module is
only imported lazily from the ``cme`` CLI command. The annotation-level region
merge additionally needs ``geopandas``/``shapely`` and is imported on demand.
"""

from __future__ import annotations

import math
import os
from collections import deque
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import as_completed
from pathlib import Path
from typing import Any
from typing import Dict
from typing import Iterable
from typing import List
from typing import Mapping
from typing import Optional
from typing import Sequence
from typing import Tuple

import click
import igraph as ig
import joblib
import leidenalg as la
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.cluster import KMeans
from sklearn.metrics import normalized_mutual_info_score
from sklearn.metrics import silhouette_score
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import StandardScaler
from torch_geometric.data import Data
from torch_geometric.loader import DataListLoader
from torch_geometric.loader import DataLoader as GeoDataLoader
from torch_geometric.nn import DataParallel as GeoDataParallel
from torch_geometric.nn import GCNConv
from torch_geometric.nn.models import DeepGraphInfomax
from tqdm import tqdm

from .. import errors
from ..cancel import critical_section
from ..cancel import raise_if_cancelled
from ..io._wsi_stub import _validate_wsi_directory
from ..io._wsi_stub import get_avg_mpp
from ..uri_path import URIPath
from .graph_cache import get_or_build_delaunay
from .insight_helpers import compute_cell_center_points
from .insight_helpers import create_adjacency_list_fast
from .insight_helpers import delaunay_triangulation

# =============================================================================
# Worker-count helpers (replacements for WSInsight's num_worker_optimizer)
# =============================================================================


def pick_workers_safe(max_workers: int | None = None, min_workers: int = 1) -> int:
    """Pick a process-pool size in ``[1, cpu_count]``.

    Targets ``max_workers`` (capped to the CPU count) but never drops below
    ``min_workers`` (itself clamped to the CPU count).
    """
    cpu = os.cpu_count() or 1
    hi = cpu if (not max_workers or max_workers < 1) else min(cpu, int(max_workers))
    lo = max(1, min(int(min_workers), cpu))
    return min(cpu, max(lo, hi))


def throttle_when_busy() -> None:
    """No-op placeholder for WSInsight's load-aware throttle."""
    return None


# =============================================================================
# Utilities: probabilities, edges, isolation
# =============================================================================


def probs_from_df(
    df: pd.DataFrame, class_order: Optional[List[str]] = None
) -> Tuple[np.ndarray, List[str]]:
    """Extract [N,C] soft probabilities from columns like 'prob_*'."""
    cols = [c for c in df.columns if c.startswith("prob_")]
    if class_order is not None:
        want = [f"prob_{k}" for k in class_order]
        missing = [c for c in want if c not in cols]
        if missing:
            raise ValueError(f"Missing probability columns: {missing}")
        cols = want
        classes = class_order
    else:
        classes = [c[len("prob_") :] for c in cols]

    P = df[cols].to_numpy(dtype=np.float32)  # [N,C]
    s = P.sum(axis=1, keepdims=True) + 1e-8
    P = P / s
    return P, classes


def to_edge_index(
    edges_df: pd.DataFrame,
    src_col: str = "source",
    dst_col: str = "target",
    undirected: bool = True,
    drop_self_loops: bool = True,
) -> np.ndarray:
    """DataFrame -> edge_index [2,E]. Assumes 0-based indices already capped."""
    u = edges_df[src_col].to_numpy()
    v = edges_df[dst_col].to_numpy()
    if drop_self_loops:
        keep = u != v
        u, v = u[keep], v[keep]
    if undirected:
        ei = np.r_[u, v]
        ej = np.r_[v, u]
    else:
        ei, ej = u, v
    return np.vstack([ei, ej]).astype(np.int64)


def drop_isolated(edge_index: np.ndarray, N: int) -> Tuple[np.ndarray, np.ndarray]:
    """Remove nodes with degree 0. Returns (edge_index_kept, kept_indices)."""
    if edge_index.size == 0:
        return edge_index, np.array([], dtype=np.int64)
    ei, ej = edge_index
    deg = np.bincount(np.r_[ei, ej], minlength=N)
    kept = np.where(deg > 0)[0]
    if len(kept) == N:
        return edge_index, kept

    # remap
    map_old2new = -np.ones(N, dtype=np.int64)
    map_old2new[kept] = np.arange(len(kept), dtype=np.int64)
    ei_m = map_old2new[ei]
    ej_m = map_old2new[ej]
    mask = (ei_m >= 0) & (ej_m >= 0)
    edge_index_new = np.vstack([ei_m[mask], ej_m[mask]]).astype(np.int64)
    return edge_index_new, kept


# =============================================================================
# k-hop composition (EXACT hop bins)
# =============================================================================


def _khop_rows_worker(
    start: int,
    end: int,
    k: int,
    alpha: float,
    P: np.ndarray,
    adj: dict,
    mode: str = "soft",
    labels: np.ndarray | None = None,
) -> np.ndarray:
    """Compute X rows [start:end) using EXACT-hop BFS on 'adj'.

    Returns X_block with shape [(end-start), (k+1)*C].

    mode="soft":
      - 0-hop: P[i]
      - h>=1 : Laplace-smoothed mean of neighbors' P at EXACT hop h
               out = (mean + alpha/C) / (1+alpha)

    mode="hard":
      - 0-hop: one-hot of argmax(P[i])
      - h>=1 : histogram proportions of argmax labels among EXACT hop h nodes,
               Dirichlet-smoothed with alpha.
    """
    _, C = P.shape
    H = end - start
    Xblk = np.zeros((H, (k + 1) * C), dtype=np.float32)

    if mode == "hard":
        if labels is None:
            labels = np.asarray(P.argmax(axis=1), dtype=np.int64)

    # 0-hop block
    if mode == "soft":
        Xblk[:, :C] = P[start:end]
    else:  # hard
        oh = np.zeros((H, C), dtype=np.float32)
        oh[np.arange(H), labels[start:end]] = 1.0
        Xblk[:, :C] = oh

    for row, i in enumerate(range(start, end)):
        # EXACT-hop BFS bins
        seen = {i}
        q = deque([(i, 0)])
        bins = [list() for _ in range(k + 1)]
        bins[0].append(i)
        while q:
            u, d = q.popleft()
            if d == k:
                continue
            for v in adj.get(u, []):
                if v in seen:
                    continue
                seen.add(v)
                nh = d + 1
                bins[nh].append(v)
                q.append((v, nh))

        # aggregate per hop
        for h in range(1, k + 1):
            idx = bins[h]
            off = h * C
            if not idx:
                Xblk[row, off : off + C] = 1.0 / C
                continue

            if mode == "soft":
                mean_prob = P[idx].mean(axis=0)
                Xblk[row, off : off + C] = (mean_prob + (alpha / C)) / (1.0 + alpha)
            else:
                counts = np.bincount(labels[idx], minlength=C).astype(np.float32)
                props = counts / counts.sum()
                Xblk[row, off : off + C] = (props + (alpha / C)) / (1.0 + alpha)

    return Xblk


def khop_features(
    P: np.ndarray,
    edge_index: np.ndarray,
    N: int,
    k: int = 2,
    alpha: float = 1.0,
    mode: str = "soft",
) -> np.ndarray:
    """Build k-hop feature blocks X of shape [N, (k+1)*C].

    mode="soft":
      0-hop: P[i]; h>=1: Laplace-smoothed mean of neighbors' P at EXACT hop h.
    mode="hard":
      0-hop: one-hot of argmax(P[i]); h>=1: histogram proportions at EXACT hop h.

    EXACT-hop rings (not <=h). Empty hop rings are filled with uniform 1/C.
    """
    N_nodes, C = P.shape
    assert N_nodes == N, "P and N mismatch"

    # No edges -> baseline blocks
    if edge_index.size == 0:
        X = np.zeros((N, (k + 1) * C), dtype=np.float32)
        if mode == "soft":
            X[:, :C] = P
        else:
            labels = P.argmax(axis=1)
            X[np.arange(N), labels] = 1.0  # 0-hop one-hot
        for h in range(1, k + 1):
            X[:, h * C : (h + 1) * C] = 1.0 / C
        return X

    # Build undirected unique edge list -> adjacency
    ei, ej = edge_index
    a = np.minimum(ei, ej)
    b = np.maximum(ei, ej)
    pairs = np.unique(np.stack([a, b], axis=1), axis=0)
    edges_df = pd.DataFrame({"source": pairs[:, 0], "target": pairs[:, 1]})
    adj = create_adjacency_list_fast(
        edges_df, dedup_neighbors=True, sort_neighbors=False
    )

    # Output buffer and 0-hop block
    X = np.zeros((N, (k + 1) * C), dtype=np.float32)
    if mode == "soft":
        X[:, :C] = P
        labels = None
    else:
        labels = P.argmax(axis=1).astype(np.int64)
        oh = np.zeros((N, C), dtype=np.float32)
        oh[np.arange(N), labels] = 1.0
        X[:, :C] = oh

    # Decide workers and chunking
    max_workers = pick_workers_safe(
        max_workers=(os.cpu_count() or 1) - 8, min_workers=8
    )
    chunk_size = max(1, math.ceil(N / max_workers))
    ranges = [(s, min(s + chunk_size, N)) for s in range(0, N, chunk_size)]

    if max_workers == 1 or len(ranges) == 1:
        X[:, :] = _khop_rows_worker(0, N, k, alpha, P, adj, mode=mode, labels=labels)
        return X

    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        futures = {
            ex.submit(_khop_rows_worker, s, e, k, alpha, P, adj, mode, labels): (s, e)
            for (s, e) in ranges
        }
        for fut in as_completed(futures):
            throttle_when_busy()
            s, e = futures[fut]
            X[s:e, :] = fut.result()

    return X


def khop_mean_features(
    V: np.ndarray, edge_index: np.ndarray, N: int, k: int = 2
) -> np.ndarray:
    """k-hop mean of node values (e.g. gene expression).

    Returns X of shape [N, (k+1)*G]:
      - 0-hop: V[i]
      - h>=1 : mean of V over the EXACT-hop-h neighbor ring (empty ring -> 0).

    Unlike :func:`khop_features` (composition probabilities) there is no
    simplex/Laplace smoothing; raw means are kept and rescaled later by the
    global ``StandardScaler``.
    """
    N_nodes, G = V.shape
    assert N_nodes == N, "V and N mismatch"
    X = np.zeros((N, (k + 1) * G), dtype=np.float32)
    X[:, :G] = V
    if edge_index.size == 0:
        return X

    ei, ej = edge_index
    a = np.minimum(ei, ej)
    b = np.maximum(ei, ej)
    pairs = np.unique(np.stack([a, b], axis=1), axis=0)
    edges_df = pd.DataFrame({"source": pairs[:, 0], "target": pairs[:, 1]})
    adj = create_adjacency_list_fast(
        edges_df, dedup_neighbors=True, sort_neighbors=False
    )

    for i in range(N):
        seen = {i}
        q = deque([(i, 0)])
        bins = [[i]] + [list() for _ in range(k)]
        while q:
            u, d = q.popleft()
            if d == k:
                continue
            for v in adj.get(u, []):
                if v in seen:
                    continue
                seen.add(v)
                bins[d + 1].append(v)
                q.append((v, d + 1))
        for h in range(1, k + 1):
            idx = bins[h]
            if idx:
                X[i, h * G : (h + 1) * G] = V[idx].mean(axis=0)
    return X


# =============================================================================
# PyG: GCN + DGI (shared across samples)
# =============================================================================


class GCLEncoder(nn.Module):
    def __init__(self, in_dim, hidden=64, out_dim=32, dropout=0.2):
        super().__init__()
        self.conv1 = GCNConv(in_dim, hidden)
        self.conv2 = GCNConv(hidden, out_dim)
        self.act = nn.ReLU()
        self.drop = nn.Dropout(dropout)

    def forward(self, x, edge_index):
        h = self.drop(self.act(self.conv1(x, edge_index)))
        z = self.conv2(h, edge_index)
        return z


class DGIModule(nn.Module):
    """Wrap DGI so it can accept a PyG Data object; encoder behavior unchanged."""

    def __init__(self, encoder: nn.Module):
        super().__init__()

        if not hasattr(encoder, "conv2") or not hasattr(encoder.conv2, "out_channels"):
            raise ValueError("Encoder must expose conv2.out_channels")
        enc_out_dim = int(encoder.conv2.out_channels)

        def summary(z, *args, **kwargs):
            return torch.sigmoid(z.mean(dim=0))

        def corruption(x_in, edge_in):
            perm = torch.randperm(x_in.size(0), device=x_in.device)
            return x_in[perm], edge_in

        self.dgi = DeepGraphInfomax(
            hidden_channels=enc_out_dim,
            encoder=encoder,
            summary=summary,
            corruption=corruption,
        )

    def forward(self, data: Data):
        return self.dgi(data.x, data.edge_index)

    def loss(self, pos_z, neg_z, s):
        hd = self.dgi.hidden_channels
        if s.ndim != 1 or s.numel() != hd:
            s = s.reshape(-1, hd).mean(dim=0)
        return self.dgi.loss(pos_z, neg_z, s)


def train_dgi_multi(slides, hidden=64, out_dim=32, epochs=300, lr=1e-3, wd=1e-4):
    """Train a shared DGI encoder across sample graphs and return embeddings."""
    ngpu = torch.cuda.device_count() if torch.cuda.is_available() else 0
    primary = torch.device("cuda:0" if ngpu > 0 else "cpu")

    in_dim = slides[0]["X"].shape[1]
    enc = GCLEncoder(in_dim, hidden, out_dim).to(primary)
    model = DGIModule(enc)

    data_list = [
        Data(
            x=torch.from_numpy(s["X"]).float(),
            edge_index=torch.from_numpy(s["edge_index"]).long(),
        )
        for s in slides
    ]

    if ngpu > 1:
        per_gpu_graphs = 1
        for cand in range(4, 0, -1):  # try 4,3,2,1
            try:
                test_bs = cand * max(1, ngpu)
                _ = (
                    DataListLoader(data_list[:test_bs], batch_size=test_bs)
                    .__iter__()
                    .__next__()
                )
                per_gpu_graphs = cand
                break
            except RuntimeError as e:
                if "out of memory" in str(e).lower():
                    torch.cuda.empty_cache()
                    continue
                else:
                    raise
        batch_size = per_gpu_graphs * max(1, ngpu)
        loader = DataListLoader(data_list, batch_size=batch_size, shuffle=True)
        model = GeoDataParallel(model, device_ids=list(range(ngpu))).to(primary)
    else:
        loader = GeoDataLoader(data_list, batch_size=1, shuffle=True)
        model = model.to(primary)

    enc_out = enc.conv2.out_channels
    if ngpu > 1:
        print(
            f"[DGI check] encoder_out_dim={enc_out}, dgi_hidden={model.module.dgi.hidden_channels}"
        )
        assert model.module.dgi.hidden_channels == enc_out
    else:
        print(
            f"[DGI check] encoder_out_dim={enc_out}, dgi_hidden={model.dgi.hidden_channels}"
        )
        assert model.dgi.hidden_channels == enc_out

    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=wd)

    for _ in tqdm(range(epochs)):
        for batch in loader:
            opt.zero_grad()
            if ngpu > 1:
                pos_z, neg_z, s = model(batch)
                loss = model.module.loss(pos_z, neg_z, s)
            else:
                batch = batch.to(primary)
                pos_z, neg_z, s = model(batch)
                loss = model.loss(pos_z, neg_z, s)
            loss.backward()
            opt.step()

    enc_eval = (model.module.dgi.encoder if ngpu > 1 else model.dgi.encoder).to(primary)
    enc_eval.eval()
    Z_list = []
    with torch.no_grad():
        for s in slides:
            x = torch.from_numpy(s["X"]).float().to(primary)
            ei = torch.from_numpy(s["edge_index"]).long().to(primary)
            Z_list.append(enc_eval(x, ei).cpu().numpy().astype(np.float32))
    return enc_eval, Z_list


# =============================================================================
# Sample-graph construction (k-hop composition only; no image features)
# =============================================================================


def prepare_slide_graph(
    cme_detection_df: pd.DataFrame,
    mpp_um_per_px: float,
    max_edge_len_um: float,
    class_order: Optional[List[str]] = None,
    k_hops: int = 2,
    alpha: float = 1.0,
    graph_cache_dir: Optional[Path] = None,
    slide_id: Optional[str] = None,
    mode: str = "hard",
    use_expression: bool = False,
    feature_source: str = "celltype",
    expr_pca=None,
    expr_pca_cols: Optional[List[str]] = None,
) -> Dict[str, np.ndarray]:
    """Build one sample graph.

    - centers from bbox (shared helper)
    - Delaunay + distance cap (shared helper); for sptxinsight ``mpp_um_per_px``
      is ``1.0`` and coordinates are already in microns.
    - drop isolated cells
    - features selected by ``feature_source``:
        * ``"celltype"``   -> k-hop cell-type composition (default; unchanged).
        * ``"expression"`` -> k-hop mean gene expression only (needs ``expr_``).
        * ``"both"``       -> composition concatenated with mean expression.
      The legacy ``use_expression`` flag is treated as ``feature_source="both"``.
      When ``expr_pca`` (a fitted PCA) is given, the per-cell expression is first
      projected onto that shared basis (using ``expr_pca_cols`` for a fixed gene
      order) so the expression feature block carries PC scores instead of raw
      genes.

    Returns: {'X','edge_index','kept_idx','classes','genes','feature_source','edges_df'}.
    """
    if use_expression and feature_source == "celltype":
        feature_source = "both"
    if feature_source not in ("celltype", "expression", "both"):
        raise ValueError(
            f"feature_source must be celltype/expression/both, got {feature_source!r}"
        )
    want_celltype = feature_source in ("celltype", "both")
    want_expression = feature_source in ("expression", "both")

    df = compute_cell_center_points(cme_detection_df.copy())
    centers_px = df[["center_x", "center_y"]].to_numpy(dtype=np.float32)
    N = len(df)

    max_edge_len_px = float(max_edge_len_um) / float(mpp_um_per_px)
    if graph_cache_dir is not None and slide_id is not None:
        centers_int = np.asarray(centers_px, dtype=np.int32)
        edges_df = get_or_build_delaunay(
            graph_cache_dir, slide_id, centers_int, mpp_um_per_px, max_edge_len_px
        )
    else:
        edges_df = delaunay_triangulation(centers_px, max_edge_len_px)

    edge_index = to_edge_index(
        edges_df,
        src_col="source",
        dst_col="target",
        undirected=True,
        drop_self_loops=True,
    )
    edge_index, kept_idx = drop_isolated(edge_index, N)
    if kept_idx.size == 0:
        raise ValueError("All nodes are isolated after distance cap; nothing to train.")
    N_kept = len(kept_idx)

    blocks: List[np.ndarray] = []
    classes: List[str] = []
    genes: List[str] = []

    if want_celltype:
        P_all, classes = probs_from_df(df, class_order=class_order)  # [N,C]
        P = P_all[kept_idx]  # [N_kept,C]
        X_khop = khop_features(
            P=P, edge_index=edge_index, N=N_kept, k=k_hops, alpha=alpha, mode=mode
        )
        blocks.append(X_khop.astype(np.float32))

    if want_expression:
        expr_cols = [c for c in df.columns if c.startswith("expr_")]
        if not expr_cols:
            if feature_source == "expression":
                raise errors.WsinferException(
                    f"--cme-mode expression needs expr_ columns, but sample "
                    f"{slide_id!r} has none. Re-ingest with transcript data, or "
                    f"use --cme-mode celltype."
                )
            # 'both' on a sample without expression: fall back to composition.
        elif expr_pca is not None:
            # Project per-cell expression onto the shared PCA basis BEFORE the
            # k-hop aggregation. ``expr_pca_cols`` fixes the gene order so every
            # sample is projected with the same components, keeping niches
            # comparable across the cohort. The feature block then carries PC
            # scores (pc0..) instead of raw genes; the model df keeps the
            # interpretable expr_ columns for cme-profile markers.
            use_cols = expr_pca_cols if expr_pca_cols is not None else expr_cols
            V = df.reindex(columns=use_cols).to_numpy(dtype=np.float32)
            V = np.nan_to_num(V)[kept_idx]
            V = expr_pca.transform(V).astype(np.float32)
            genes = [f"pc{i}" for i in range(V.shape[1])]
            X_expr = khop_mean_features(V, edge_index=edge_index, N=N_kept, k=k_hops)
            blocks.append(X_expr.astype(np.float32))
        else:
            genes = [c[len("expr_") :] for c in expr_cols]
            V = df[expr_cols].to_numpy(dtype=np.float32)
            V = np.nan_to_num(V)[kept_idx]
            X_expr = khop_mean_features(V, edge_index=edge_index, N=N_kept, k=k_hops)
            blocks.append(X_expr.astype(np.float32))

    if not blocks:
        raise ValueError(
            f"No features built for sample {slide_id!r} (mode={feature_source})."
        )

    X = np.hstack(blocks).astype(np.float32)

    return {
        "X": X,
        "edge_index": edge_index.astype(np.int64),
        "kept_idx": kept_idx.astype(np.int64),
        "classes": classes,
        "genes": genes,
        "feature_source": feature_source,
        "edges_df": edges_df,
    }


# =============================================================================
# Clustering: Leiden sweep / KMeans on DGI embeddings
# =============================================================================


def _knn_graph_connectivity(Z: np.ndarray, k_nn: int = 15):
    A = kneighbors_graph(Z, n_neighbors=k_nn, mode="connectivity", include_self=False)
    A = A.maximum(A.T).tocsr()  # symmetrize
    return A


def _igraph_from_sparse(A) -> ig.Graph:
    """Convert a scipy sparse adjacency matrix to an undirected igraph graph."""
    A = A.tocoo()
    g = ig.Graph(
        n=A.shape[0],
        edges=list(zip(A.row.tolist(), A.col.tolist(), strict=False)),
        directed=False,
    )
    g.simplify(combine_edges="ignore")
    return g


def _leiden_worker(
    n_nodes: int, edges: np.ndarray, resolution: float
) -> Tuple[np.ndarray, float, float]:
    """Run a single Leiden clustering pass and return labels plus modularity."""
    g_local = ig.Graph(n=n_nodes, edges=edges.tolist(), directed=False)
    g_local.simplify(combine_edges="ignore")
    part = la.find_partition(
        g_local,
        la.RBConfigurationVertexPartition,
        resolution_parameter=float(resolution),
    )
    labels = np.asarray(part.membership, dtype=int)
    return labels, float(part.modularity), float(resolution)


def _reduce_resolution_worker(args):
    """Summarize repeated Leiden runs for one resolution value."""
    r, runs, Z = args
    best_labels, best_mod = max(runs, key=lambda x: x[1])

    nmis = []
    if len(np.unique(best_labels)) > 1:
        for lab, _ in runs:
            if len(np.unique(lab)) > 1:
                nmis.append(normalized_mutual_info_score(lab, best_labels))
    stability = float(np.mean(nmis)) if nmis else 0.0

    if len(np.unique(best_labels)) > 1:
        sil = float(
            silhouette_score(
                Z, best_labels, sample_size=np.min([len(Z), 10000]), metric="euclidean"
            )
        )
    else:
        sil = -1.0

    counts = np.bincount(best_labels)
    min_frac = float(counts.min() / counts.sum()) if counts.size else 0.0

    return {
        "resolution": float(r),
        "n_clusters": int(len(np.unique(best_labels))),
        "modularity": float(best_mod),
        "stability": stability,
        "silhouette": sil,
        "min_frac": min_frac,
        "labels": best_labels,
    }


def _leiden_sweep_on_graph(
    Z: np.ndarray,
    g: ig.Graph,
    cme_clustering_resolutions: Iterable[float] = np.arange(0.2, 2.05, 0.1),
    n_repeats: int = 5,
) -> Dict[str, Any]:
    """Parallel sweep over (resolution, repeat) and parallel reduction per resolution.

    Returns {"winner": {...}, "all": [ per-resolution dicts ... ]}.
    """
    n_nodes = g.vcount()
    el = np.array(g.get_edgelist(), dtype=np.int64)
    if el.size == 0:
        labels = np.zeros(n_nodes, dtype=int)
        out = {
            "resolution": float(next(iter(cme_clustering_resolutions), 1.0)),
            "n_clusters": 1,
            "modularity": 0.0,
            "stability": 1.0,
            "silhouette": -1.0,
            "min_frac": 1.0,
            "labels": labels,
        }
        return {"winner": out, "all": [out]}

    # ---- Phase A: parallel Leiden runs over (resolution, repeat) ----
    tasks = []
    for r in cme_clustering_resolutions:
        for _ in range(n_repeats):
            tasks.append((n_nodes, el, float(r)))

    n_jobs = pick_workers_safe(max_workers=(os.cpu_count() or 1) - 2, min_workers=2)
    results_by_r: Dict[float, list] = {}
    with ProcessPoolExecutor(max_workers=n_jobs) as ex:
        futs = [ex.submit(_leiden_worker, *t) for t in tasks]
        for fut in as_completed(futs):
            throttle_when_busy()
            labels, modularity, r = fut.result()
            results_by_r.setdefault(r, []).append((labels, modularity))

    # ---- Phase B: parallel reduction per resolution ----
    logs = []
    with ProcessPoolExecutor(max_workers=n_jobs) as ex:
        futs = [
            ex.submit(_reduce_resolution_worker, (r, results_by_r[r], Z))
            for r in results_by_r.keys()
        ]
        for fut in as_completed(futs):
            throttle_when_busy()
            logs.append(fut.result())

    logs.sort(key=lambda d: d["resolution"])

    filtered = [d for d in logs if d["min_frac"] >= 0.005] or logs
    winner = sorted(
        filtered,
        key=lambda d: (d["stability"], d["modularity"], d["silhouette"]),
        reverse=True,
    )[0]

    return {"winner": winner, "all": logs}


def estimate_cmes_from_Z_list(
    Z_list: List[np.ndarray],
    mode: str = "global",  # "global" (recommended) or "per_slide"
    k_nn: int = 15,
    cme_clustering_resolutions=np.arange(0.2, 2.05, 0.1),  # noqa: B008
    n_repeats: int = 5,
) -> Dict[str, Any]:
    """Auto-estimate the number of CMEs via a Leiden community-detection sweep."""
    if mode == "global":
        offsets = np.cumsum([0] + [Z.shape[0] for Z in Z_list[:-1]])
        Z_all = np.vstack(Z_list)
        A = _knn_graph_connectivity(Z_all, k_nn=k_nn)
        g = _igraph_from_sparse(A)
        sweep = _leiden_sweep_on_graph(
            Z_all,
            g,
            cme_clustering_resolutions=cme_clustering_resolutions,
            n_repeats=n_repeats,
        )
        w = sweep["winner"]
        labels_all = w["labels"]
        labels_list = []
        for off, Z in zip(offsets, Z_list, strict=False):
            labels_list.append(labels_all[off : off + len(Z)])
        return {
            "clusters_k": w["n_clusters"],
            "labels_list": labels_list,
            "winner": w,
            "all_results": sweep["all"],
        }

    elif mode == "per_slide":
        labels_list = []
        winners = []
        all_logs = []
        n_clusters_list = []
        for Z in Z_list:
            A = _knn_graph_connectivity(Z, k_nn=k_nn)
            g = _igraph_from_sparse(A)
            sweep = _leiden_sweep_on_graph(
                Z,
                g,
                cme_clustering_resolutions=cme_clustering_resolutions,
                n_repeats=n_repeats,
            )
            w = sweep["winner"]
            labels_list.append(w["labels"])
            winners.append(w)
            all_logs.append(sweep["all"])
            n_clusters_list.append(w["n_clusters"])
        return {
            "clusters_k": int(np.median(n_clusters_list)),
            "labels_list": labels_list,
            "winner": winners,
            "all_results": all_logs,
        }
    else:
        raise ValueError("mode must be 'global' or 'per_slide'")


def _mpp_for(wsi_path, slide_mpp_lookup: Mapping[str, float] | None) -> float:
    """Resolve microns-per-pixel for a sample.

    For sptxinsight ``slide_mpp_lookup`` supplies ``1.0`` (coordinates already in
    microns). Falls back to the WSI stub only if no lookup entry exists, which
    raises a clear error.
    """
    slide_id = Path(str(wsi_path)).stem
    if slide_mpp_lookup:
        mpp = slide_mpp_lookup.get(slide_id) or slide_mpp_lookup.get(str(wsi_path))
        if mpp:
            return float(mpp)
    return float(get_avg_mpp(wsi_path))


# =============================================================================
# CME mode namespacing + cross-sample batch correction
# =============================================================================

# Each mode writes to its own output dir / one-hot column prefix / checkpoint
# suffix so cell-type, gene-expression, and hybrid niches can coexist on the
# same cells without clobbering each other. ``celltype`` keeps the original
# (unsuffixed) paths and ``cme_`` prefix for backward compatibility.
_CME_MODE_SPEC: Dict[str, Dict[str, str]] = {
    "celltype": {"subdir": "cme-outputs-csv", "prefix": "cme", "ckpt": ""},
    "expression": {"subdir": "cme-gex-outputs-csv", "prefix": "gexcme", "ckpt": "-gex"},
    "both": {"subdir": "cme-hybrid-outputs-csv", "prefix": "hcme", "ckpt": "-hybrid"},
}


def center_per_sample(Z_list: List[np.ndarray]) -> List[np.ndarray]:
    """Per-sample mean-centering of DGI embeddings (native batch correction).

    Subtracts each sample's mean embedding and adds back the cohort grand mean,
    so per-sample location shifts (a common technical batch effect) are removed
    while the shared embedding geometry is preserved. Zero dependencies.
    """
    if not Z_list:
        return Z_list
    grand = np.vstack(Z_list).mean(axis=0)
    return [(Z - Z.mean(axis=0) + grand).astype(np.float32) for Z in Z_list]


def _harmony_correct(
    Z_list: List[np.ndarray], sample_ids: Sequence[str]
) -> List[np.ndarray]:
    """Harmony batch correction over the pooled DGI embeddings (optional dep)."""
    try:
        from harmonypy import run_harmony
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise RuntimeError(
            "Harmony batch correction needs the optional 'harmonypy' package. "
            "Install it with: pip install harmonypy  "
            "(or pip install 'sptxinsight[harmony]')."
        ) from exc
    lengths = [len(Z) for Z in Z_list]
    Z_all = np.vstack(Z_list).astype(np.float64)
    batch = np.repeat([str(s) for s in sample_ids], lengths)
    meta = pd.DataFrame({"sample": batch})
    ho = run_harmony(Z_all, meta, ["sample"])
    Z_corr = np.asarray(ho.Z_corr).T.astype(np.float32)  # [N, d]
    out: List[np.ndarray] = []
    off = 0
    for n in lengths:
        out.append(Z_corr[off : off + n])
        off += n
    return out


def _apply_batch_correction(
    Z_list: List[np.ndarray], batch_correct: str, sample_ids: Sequence[str]
) -> List[np.ndarray]:
    """Dispatch ``none`` / ``center`` / ``harmony`` correction on ``Z_list``."""
    if batch_correct in (None, "none"):
        return Z_list
    if batch_correct == "center":
        return center_per_sample(Z_list)
    if batch_correct == "harmony":
        return _harmony_correct(Z_list, sample_ids)
    raise ValueError(
        f"batch_correct must be none/center/harmony, got {batch_correct!r}"
    )


def _fit_expression_pca(model_output_paths: Sequence[Path], n_components: int):
    """Fit one shared PCA basis on the pooled per-cell expression of a cohort.

    Gene expression (``expr_`` columns) is high-dimensional and highly
    correlated, so feeding all genes (times the hop count) into the GCN is
    redundant and noisy. This reduces every cell's expression to a small set of
    principal components BEFORE the k-hop aggregation. The basis is fit once on
    the pooled cells of all samples (``IncrementalPCA`` so the full cohort is
    never held in memory at once) and then applied identically to each sample,
    which is what keeps the resulting niches comparable across the cohort.

    Returns ``(ipca, expr_cols)`` with ``expr_cols`` fixing the gene order, or
    ``None`` when PCA is disabled or no sample carries ``expr_`` columns.
    """
    if not n_components or n_components < 1:
        return None
    from sklearn.decomposition import IncrementalPCA

    expr_cols: Optional[List[str]] = None
    for p in model_output_paths:
        cols = sorted(
            c for c in pd.read_csv(p, nrows=0).columns if c.startswith("expr_")
        )
        if cols:
            expr_cols = cols
            break
    if not expr_cols:
        return None

    n_comp = min(int(n_components), len(expr_cols))
    ipca = IncrementalPCA(n_components=n_comp)
    fitted_any = False
    for p in model_output_paths:
        head = pd.read_csv(p, nrows=0).columns
        if not any(c.startswith("expr_") for c in head):
            continue
        V = (
            pd.read_csv(p, usecols=expr_cols)
            .reindex(columns=expr_cols)
            .to_numpy(np.float32)
        )
        V = np.nan_to_num(V)
        if V.shape[0] < n_comp:  # IncrementalPCA needs >= n_components rows per batch
            continue
        ipca.partial_fit(V)
        fitted_any = True
    if not fitted_any:
        return None
    return ipca, expr_cols


# =============================================================================
# End-to-end multi-sample training + clustering
# =============================================================================


def cme_generation(
    wsi_dir: str | URIPath | None,
    wsi_paths: Sequence[Path | URIPath] | None,
    results_dir: str | Path,
    max_edge_len_um: float,
    max_cell_radius_um: float,
    class_order: Optional[List[str]] = None,
    k_hops: int = 2,
    alpha: float = 1.0,
    # encoder
    hidden: int = 64,
    out_dim: int = 32,
    epochs: int = 300,
    # clustering
    cme_cellular: bool = False,
    cme_annotation: bool = False,
    cme_clustering_k: int | None = 10,
    cme_clustering_resolutions: List[float] = [0.5, 1.0, 2.0],  # noqa: B006
    cme_soft_mode: bool = False,
    use_expression: bool = False,
    cme_mode: str = "celltype",
    batch_correct: str = "none",
    expression_pca: int = 50,
    overwrite: bool = False,
    slide_mpp_lookup: Mapping[str, float] | None = None,
) -> None:
    """Discover CMEs across a cohort of spatial-transcriptomics samples.

    Builds per-sample Delaunay cell graphs, trains one shared DGI encoder, clusters
    the embeddings, and writes per-cell CME labels (and optional region merges).

    ``cme_mode`` selects which features drive the niches and namespaces all
    outputs/checkpoints (see :data:`_CME_MODE_SPEC`): ``celltype`` (default,
    byte-identical to the original behaviour), ``expression`` (gene niches), or
    ``both`` (fused). The legacy ``use_expression`` flag maps to ``both``.
    ``batch_correct`` (``none``/``center``/``harmony``) is applied to the DGI
    embeddings before clustering to remove cross-sample technical shifts.
    ``expression_pca`` (>0, expression/both modes) reduces the per-cell gene
    panel to that many shared principal components before the k-hop aggregation;
    set it to 0 to feed all genes in.
    """
    if use_expression and cme_mode == "celltype":
        cme_mode = "both"
    if cme_mode not in _CME_MODE_SPEC:
        raise ValueError(
            f"cme_mode must be one of {sorted(_CME_MODE_SPEC)}, got {cme_mode!r}"
        )
    mode_spec = _CME_MODE_SPEC[cme_mode]

    if os.getenv("SPTXINSIGHT_FORCE_CPU", "0").lower() not in {"0", "f", "false"}:
        device = torch.device("cpu")
    elif torch.cuda.is_available():
        device = torch.device("cuda")
    elif torch.backends.mps.is_available() and torch.backends.mps.is_built():
        device = torch.device("mps")
    else:
        device = torch.device("cpu")
    print(f'Using device "{device}"')

    wsi_dir_path = URIPath(wsi_dir) if wsi_dir is not None else None
    if wsi_dir_path is not None and not wsi_dir_path.exists():
        raise errors.WholeSlideImageDirectoryNotFound(
            f"directory not found: {wsi_dir_path}"
        )

    if wsi_paths is not None:
        slide_paths = [p if isinstance(p, URIPath) else URIPath(p) for p in wsi_paths]
    elif wsi_dir_path is not None:
        slide_paths = [p for p in wsi_dir_path.iterdir() if p.is_file()]
    else:
        raise ValueError("wsi_paths must be provided when wsi_dir is None")

    if not slide_paths:
        context = wsi_dir_path or "provided sample paths"
        raise errors.WholeSlideImagesNotFound(context)

    results_dir = Path(results_dir)
    if not results_dir.exists():
        raise errors.ResultsDirectoryNotFound(results_dir)

    if wsi_dir_path is not None:
        _validate_wsi_directory(wsi_dir_path)
    else:
        stems = [p.stem for p in slide_paths]
        if len(stems) != len(set(stems)):
            raise errors.DuplicateFilePrefixesFound(
                "A sample with the same prefix but different extensions has been found"
            )

    model_output_dir = results_dir / "model-outputs-csv"
    if not model_output_dir.exists():
        raise errors.ResultsDirectoryNotFound(
            "The 'model-outputs-csv' directory was not found in results directory."
        )
    model_output_paths = [
        model_output_dir / p.with_suffix(".csv").name for p in slide_paths
    ]

    if len(model_output_paths) != len(slide_paths):
        raise errors.ResultsDirectoryNotFound(
            "The 'model-outputs-csv' and sample directory were mismatched."
        )

    cme_output_dir = results_dir / mode_spec["subdir"]
    cme_output_dir.mkdir(exist_ok=True)
    cme_cells_output_dir = cme_output_dir / "cells"
    cme_cells_output_dir.mkdir(exist_ok=True)
    cme_cmes_output_dir = cme_output_dir / "cmes"
    cme_cmes_output_dir.mkdir(exist_ok=True)
    # PCA only reduces expression features; tag the checkpoints so toggling
    # --disable-pca never silently reuses an incompatible cached graph/embedding.
    pca_on = cme_mode in ("expression", "both") and bool(expression_pca)
    pca_tag = f"-pca{int(expression_pca)}" if pca_on else ""
    cme_slide_graph_file = (
        results_dir / f"slide-graphs{mode_spec['ckpt']}{pca_tag}.joblib"
    )
    cme_dgi_embeddings_file = (
        results_dir / f"dgi-embeddings{mode_spec['ckpt']}{pca_tag}.joblib"
    )

    if overwrite:
        for _ckpt in (cme_slide_graph_file, cme_dgi_embeddings_file):
            if _ckpt.exists():
                _ckpt.unlink()

    # ---- Phase 1/5: build sample graphs ------------------------------------
    slides: list = []
    classes = None

    if cme_slide_graph_file.exists():
        click.secho(
            "\nPhase 1/5: Build sample graphs for CMEGCN.\n"
            f"Load existing graph file: {cme_slide_graph_file}\n",
            fg="green",
        )
        slides = joblib.load(cme_slide_graph_file)
    else:
        click.secho("\nPhase 1/5: build sample graphs for CMEGCN.\n", fg="green")

        graph_cache_dir = Path(str(results_dir)) / "graphs"
        graph_cache_dir.mkdir(parents=True, exist_ok=True)

        # Fit one shared PCA basis on the pooled cohort expression (expression /
        # both modes only). Done once, up front, so every sample is projected
        # with the same components before k-hop aggregation.
        expr_pca = None
        expr_pca_cols = None
        if pca_on:
            click.secho(
                f"Fitting shared expression PCA ({int(expression_pca)} comps) "
                f"on pooled cohort...",
                fg="green",
            )
            fitted = _fit_expression_pca(model_output_paths, int(expression_pca))
            if fitted is None:
                click.secho(
                    "  No expr_ columns found; continuing without PCA.", fg="yellow"
                )
            else:
                expr_pca, expr_pca_cols = fitted
                evr = float(np.sum(expr_pca.explained_variance_ratio_))
                click.secho(
                    f"  PCA basis: {expr_pca.n_components_} comps over "
                    f"{len(expr_pca_cols)} genes (cum. EVR={evr:.3f}).",
                    fg="green",
                )

        # Sequential over samples so khop_features' inner process pool is not
        # nested inside an outer one.
        for wsi_path, csv_path in tqdm(
            list(zip(slide_paths, model_output_paths, strict=False)),
            total=len(slide_paths),
        ):
            raise_if_cancelled()
            model_output_df = pd.read_csv(csv_path)
            slide_id = Path(str(wsi_path)).stem
            mpp = _mpp_for(wsi_path, slide_mpp_lookup)
            s = prepare_slide_graph(
                model_output_df,
                mpp_um_per_px=mpp,
                max_edge_len_um=max_edge_len_um,
                class_order=class_order,
                k_hops=k_hops,
                alpha=alpha,
                graph_cache_dir=graph_cache_dir,
                slide_id=slide_id,
                mode="soft" if cme_soft_mode else "hard",
                feature_source=cme_mode,
                expr_pca=expr_pca,
                expr_pca_cols=expr_pca_cols,
            )
            slides.append(s)
            if classes is None:
                classes = s["classes"]

        # Global z-score over concatenated features (consistent scale across samples)
        X_all = np.vstack([s["X"] for s in slides]).astype(np.float32)
        scaler = StandardScaler(with_mean=True, with_std=True).fit(X_all)
        for s in slides:
            s["X_normalized"] = scaler.transform(s["X"]).astype(np.float32)

        joblib.dump(slides, cme_slide_graph_file, compress=3)

    # ---- Phase 2/5: shared DGI encoder + embeddings ------------------------
    if cme_dgi_embeddings_file.exists():
        click.secho(
            "\nPhase 2/5: Train shared DGI encoder and get embeddings per sample.\n"
            f"Load existing embeddings file: {cme_dgi_embeddings_file}\n",
            fg="green",
        )
        Z_list = joblib.load(cme_dgi_embeddings_file)
    else:
        click.secho(
            "\nPhase 2/5: Train shared DGI encoder and get embeddings per sample.\n",
            fg="green",
        )
        _, Z_list = train_dgi_multi(
            slides, hidden=hidden, out_dim=out_dim, epochs=epochs
        )
        joblib.dump(Z_list, cme_dgi_embeddings_file, compress=3)

    # ---- Phase 2b/5: cross-sample batch correction on embeddings -----------
    # Applied AFTER (re)loading the raw embeddings so the checkpoint stays
    # correction-agnostic and the method can be changed between runs.
    if batch_correct not in (None, "none"):
        click.secho(
            f"\nPhase 2b/5: Apply '{batch_correct}' cross-sample batch correction.\n",
            fg="green",
        )
        sample_ids = [Path(str(p)).stem for p in slide_paths]
        Z_list = _apply_batch_correction(Z_list, batch_correct, sample_ids)

    # ---- Phase 3/5: choose cluster count -----------------------------------
    if not cme_clustering_k:
        click.secho("\nPhase 3/5: Estimate CME cluster number.\n", fg="green")
        est = estimate_cmes_from_Z_list(
            Z_list,
            mode="global",
            cme_clustering_resolutions=cme_clustering_resolutions,
            k_nn=15,
        )
        cme_clustering_k = est["winner"]["n_clusters"]
        labels_list = est["labels_list"]
    else:
        click.secho(
            f"\nPhase 3/5: Use predefined CME cluster number: cme_clustering_k={cme_clustering_k}.\n",
            fg="green",
        )
        labels_list = [
            KMeans(n_clusters=cme_clustering_k, n_init="auto")
            .fit_predict(Z)
            .astype(np.int32)
            for Z in Z_list
        ]

    # ---- Phase 4/5: cellular-level CME labels ------------------------------
    click.secho(
        "\nPhase 4/5: Perform cellular-level CME analysis per sample.\n", fg="green"
    )
    if cme_cellular:
        for i, (wsi_path, model_output_csv) in tqdm(
            enumerate(zip(slide_paths, model_output_paths, strict=False)),
            total=len(slide_paths),
        ):
            raise_if_cancelled()
            cme_csv_name = Path(str(wsi_path)).with_suffix(".csv").name
            cell_csv = cme_cells_output_dir / cme_csv_name

            if not overwrite and cell_csv.exists():
                continue

            cme_detection_df = pd.read_csv(model_output_csv)

            slide_classes = slides[i]["classes"]
            slide_genes = slides[i].get("genes", [])
            feature_normalized_cols = [
                f"feature_normalized_k{k}_{c.replace('prob_', '')}"
                for k in range(k_hops + 1)
                for c in slide_classes
            ] + [
                f"feature_normalized_k{k}_expr_{g}"
                for k in range(k_hops + 1)
                for g in slide_genes
            ]
            feature_cols = [
                f"feature_raw_k{k}_{c.replace('prob_', '')}"
                for k in range(k_hops + 1)
                for c in slide_classes
            ] + [
                f"feature_raw_k{k}_expr_{g}"
                for k in range(k_hops + 1)
                for g in slide_genes
            ]
            cme_detection_df.loc[
                slides[i]["kept_idx"], feature_normalized_cols
            ] = slides[i]["X_normalized"]
            cme_detection_df.loc[slides[i]["kept_idx"], feature_cols] = slides[i]["X"]
            cme_cols = [
                f"{mode_spec['prefix']}_{idx}" for idx in range(cme_clustering_k)
            ]
            label_one_hot = np.eye(cme_clustering_k, dtype=np.float32)[labels_list[i]]
            cme_detection_df.loc[slides[i]["kept_idx"], cme_cols] = label_one_hot

            with critical_section(
                f"saving CME cell output for {Path(str(wsi_path)).stem}"
            ):
                cme_detection_df.to_csv(cell_csv, index=False)

    # ---- Phase 5/5: annotation-level CME regions ---------------------------
    click.secho(
        "\nPhase 5/5: Perform annotation-level CME analysis per sample.\n", fg="green"
    )
    if cme_annotation:
        try:
            from .vorononi_cme_region_helper import (
                merge_same_label_by_shared_edges_iterative,
            )
            from .vorononi_cme_region_helper import remap_edges_to_valid_indices
        except ImportError as exc:
            raise RuntimeError(
                "Annotation-level CME regions need geopandas and shapely. "
                "Install them with: pip install geopandas shapely."
            ) from exc

        for i, (wsi_path, _model_output_csv) in tqdm(
            enumerate(zip(slide_paths, model_output_paths, strict=False)),
            total=len(slide_paths),
        ):
            raise_if_cancelled()
            cme_csv_name = Path(str(wsi_path)).with_suffix(".csv").name
            cell_csv = cme_cells_output_dir / cme_csv_name
            cme_csv = cme_cmes_output_dir / cme_csv_name

            if not overwrite and cme_csv.exists():
                continue

            mpp = _mpp_for(wsi_path, slide_mpp_lookup)
            cme_detection_df = pd.read_csv(cell_csv)
            valid_mask = np.zeros(len(cme_detection_df), dtype=bool)
            valid_mask[np.asarray(slides[i]["kept_idx"], dtype=int)] = True
            edges_df = remap_edges_to_valid_indices(slides[i]["edges_df"], valid_mask)

            cme_annotation_df = merge_same_label_by_shared_edges_iterative(
                cme_detection_df,
                edges_df,
                cme_clustering_k=cme_clustering_k,
                mpp=mpp,
                max_radius_um=max_cell_radius_um,
                cme_prefix=f"{mode_spec['prefix']}_",
            )

            with critical_section(
                f"saving CME annotation output for {Path(str(wsi_path)).stem}"
            ):
                cme_annotation_df.to_csv(cme_csv, index=False)
