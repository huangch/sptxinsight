from __future__ import annotations

import pandas as pd
import pytest

from sptxinsight.insightlib.insight_helpers import calculate_distance_to_border
from sptxinsight.insightlib.insight_helpers import compute_hplot


def test_distance_to_border_exposes_explicit_hop_layers():
    nodes = pd.DataFrame(
        {
            "is_base_region": [True, True, False, False],
            "is_base_border": [False, True, False, False],
            "is_base_type": [True, True, False, False],
            "is_target_type": [False, False, True, True],
        }
    )
    adjacency = {0: [1], 1: [0, 2], 2: [1, 3], 3: [2]}

    out = calculate_distance_to_border(nodes.copy(), adjacency)

    assert out["distance_to_border_hops"].tolist() == [1.0, 0.0, 1.0, 2.0]
    assert out["signed_distance_to_border_hops"].tolist() == [-1.0, -0.0, 1.0, 2.0]
    assert out["hplot_layer"].tolist() == [-1.0, -0.0, 1.0, 2.0]
    assert out["distance_to_border"].equals(out["distance_to_border_hops"])
    assert out["signed_distance_to_border"].equals(out["signed_distance_to_border_hops"])


def test_compute_hplot_separates_layer_from_distance_um():
    nodes = pd.DataFrame(
        {
            "hplot_layer": [-1.0, 0.0, 1.0, 2.0],
            "is_base_type": [True, True, False, False],
            "is_target_type": [False, False, True, True],
        }
    )
    edges = pd.DataFrame(
        {
            "source": [0, 1, 2],
            "target": [1, 2, 3],
            "length": [10.0, 20.0, 30.0],
        }
    )

    hplot = compute_hplot(nodes, edges)

    assert hplot["layer"].tolist() == [-1.0, 0.0, 1.0, 2.0]
    assert hplot["distance_um"].tolist() == pytest.approx([-10.0, 0.0, 20.0, 50.0])
    assert hplot["distance"].tolist() == pytest.approx(hplot["distance_um"].tolist())
    assert hplot.loc[hplot["layer"] == 1.0, "target_type_prop"].item() == pytest.approx(1.0)