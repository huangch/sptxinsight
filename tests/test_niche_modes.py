"""Tests for the parallel niche modes (celltype / expression / both) and the
cross-sample batch correction + agreement helpers added to ``sptxinsight niche``.

The niche_profile helpers are numpy/pandas-only; the batch-correction helpers live
in niche_generation (which imports torch) and are tested behind a skip guard.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from sptxinsight.insightlib.niche_profile import _MODE_SPEC
from sptxinsight.insightlib.niche_profile import _niche_columns
from sptxinsight.insightlib.niche_profile import _nmi
from sptxinsight.insightlib.niche_profile import niche_agreement
from sptxinsight.insightlib.niche_profile import niche_profile

# --------------------------------------------------------------------------- #
# Mode-spec namespacing
# --------------------------------------------------------------------------- #


def test_mode_spec_is_distinct_and_celltype_unchanged():
    subdirs = {m: s[0] for m, s in _MODE_SPEC.items()}
    prefixes = {m: s[1] for m, s in _MODE_SPEC.items()}
    # celltype keeps the original (unsuffixed) namespace for backward compat.
    assert _MODE_SPEC["celltype"] == ("niche-outputs-csv", "niche_", "")
    # Every mode writes to a distinct folder + column prefix.
    assert len(set(subdirs.values())) == len(subdirs)
    assert len(set(prefixes.values())) == len(prefixes)


def test_niche_columns_orders_by_index_for_each_prefix():
    cols = ["x", "gexniche_10", "gexniche_2", "gexniche_1", "niche_0"]
    assert _niche_columns(cols, "gexniche_") == ["gexniche_1", "gexniche_2", "gexniche_10"]
    assert _niche_columns(cols, "niche_") == ["niche_0"]


# --------------------------------------------------------------------------- #
# NMI
# --------------------------------------------------------------------------- #


def test_nmi_identical_labels_is_one():
    a = np.array([0, 0, 1, 1, 2, 2])
    assert _nmi(a, a.copy()) == pytest.approx(1.0, abs=1e-9)


def test_nmi_independent_labels_is_near_zero():
    rng = np.random.default_rng(0)
    a = rng.integers(0, 3, size=5000)
    b = rng.integers(0, 3, size=5000)
    assert _nmi(a, b) < 0.05


# --------------------------------------------------------------------------- #
# niche_profile + niche_agreement on a tiny synthetic results dir
# --------------------------------------------------------------------------- #


def _write_cells(path, prefix, labels, n_types=2):
    """Write a minimal cells CSV: prob_ columns + one-hot <prefix> columns."""
    labels = np.asarray(labels)
    k = int(labels.max()) + 1
    df = pd.DataFrame()
    for t in range(n_types):
        df[f"prob_type{t}"] = (labels % n_types == t).astype(float)
    oh = np.eye(k)[labels]
    for j in range(k):
        df[f"{prefix}{j}"] = oh[:, j]
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def test_niche_profile_reads_gene_mode_namespace(tmp_path):
    labels = [0, 0, 1, 1, 2, 0]
    _write_cells(
        tmp_path / "niche-gex-outputs-csv" / "cells" / "s1.csv", "gexniche_", labels
    )
    comp, markers = niche_profile(tmp_path, mode="expression", write=True)
    assert list(comp.index) == ["gexniche_0", "gexniche_1", "gexniche_2"]
    assert comp.loc["gexniche_0", "n_cells"] == 3
    assert (tmp_path / "niche-profile-composition-gex.csv").exists()


def test_niche_agreement_writes_crosstab(tmp_path):
    # Same six cells, two labelings that agree perfectly -> NMI == 1.
    ct = [0, 0, 1, 1, 2, 2]
    gx = [1, 1, 0, 0, 2, 2]
    _write_cells(tmp_path / "niche-outputs-csv" / "cells" / "s1.csv", "niche_", ct)
    _write_cells(tmp_path / "niche-gex-outputs-csv" / "cells" / "s1.csv", "gexniche_", gx)
    result = niche_agreement(tmp_path, write=True)
    assert result is not None
    nmi, crosstab = result
    assert nmi == pytest.approx(1.0, abs=1e-9)
    assert (tmp_path / "niche-agreement.csv").exists()
    assert list(crosstab.index) == ["niche_0", "niche_1", "niche_2"]


def test_niche_agreement_returns_none_without_both_families(tmp_path):
    _write_cells(tmp_path / "niche-outputs-csv" / "cells" / "s1.csv", "niche_", [0, 1, 0])
    assert niche_agreement(tmp_path, write=False) is None


# --------------------------------------------------------------------------- #
# Batch correction (torch-backed module; skip if unavailable)
# --------------------------------------------------------------------------- #


def test_center_per_sample_equalizes_means_preserves_grand_mean():
    pytest.importorskip("torch")
    from sptxinsight.insightlib.niche_generation import center_per_sample

    rng = np.random.default_rng(1)
    Z_list = [
        rng.normal(loc=5.0, size=(40, 3)).astype(np.float32),
        rng.normal(loc=-2.0, size=(60, 3)).astype(np.float32),
    ]
    grand_before = np.vstack(Z_list).mean(axis=0)
    out = center_per_sample(Z_list)
    means = [Z.mean(axis=0) for Z in out]
    # All per-sample means collapse onto the shared grand mean.
    assert np.allclose(means[0], means[1], atol=1e-4)
    # Grand mean is preserved.
    assert np.allclose(np.vstack(out).mean(axis=0), grand_before, atol=1e-4)
