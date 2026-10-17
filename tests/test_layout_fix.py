"""Hand-built layouts for the node/edge overlap removal."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from layout_fix import find_overlaps, resolve_overlaps  # noqa: E402

R = np.full(4, 0.1)


def test_node_on_edge_is_moved_off() -> None:
    pos = np.array([[0.0, 0.0], [4.0, 0.0], [2.0, 0.02], [2.0, 3.0]])  # node 2 lies on edge 0-1
    edges = [(0, 1), (2, 3)]
    assert find_overlaps(pos, edges, R, 0.05)[0] == [(2, 0)]
    new, left = resolve_overlaps(pos, edges, R, 0.05)
    assert left == 0 and not any(find_overlaps(new, edges, R, 0.05))
    assert (new[[0, 1, 3]] == pos[[0, 1, 3]]).all()  # only the offending node moved


def test_clean_layout_unchanged_and_incident_edges_ignored() -> None:
    pos = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 1.0]])  # node 1 touches only its own edges
    edges = [(0, 1), (1, 2), (2, 3)]
    new, left = resolve_overlaps(pos, edges, R, 0.05)
    assert left == 0 and (new == pos).all()
