"""Hand-built label layouts for the label overlap removal."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from label_fix import place_labels, rect_hits_segments  # noqa: E402

NO_DISCS = np.zeros((0, 3))
NO_RECTS = np.zeros((0, 4))
BOUNDS = np.array([0.0, 0.0, 200.0, 200.0])


def test_segment_rect_intersection() -> None:
    rect = np.array([[10.0, 10.0, 20.0, 20.0]])
    segs = np.array([[[0.0, 15.0], [30.0, 15.0]],   # through
                     [[0.0, 0.0], [5.0, 30.0]],     # beside (bbox overlaps in y, line misses)
                     [[0.0, 40.0], [30.0, 40.0]]])  # far away
    assert rect_hits_segments(rect, segs).tolist() == [[True, False, False]]


def test_label_on_foreign_edge_is_moved_own_edge_is_fine() -> None:
    segs = np.array([[[0.0, 100.0], [200.0, 100.0]], [[0.0, 50.0], [200.0, 50.0]]])
    size = np.array([[20.0, 10.0], [20.0, 10.0]])
    cands = [np.array([[100.0, 100.0], [100.0, 130.0]]),  # label of edge 1, default sits on edge 0
             np.array([[100.0, 50.0]])]                    # label of edge 1 on its own edge
    pos, before, left = place_labels(size, cands, [1, 1], segs, NO_DISCS, NO_RECTS, BOUNDS)
    assert before == 1 and left == 0 and pos.tolist() == [[100.0, 130.0], [100.0, 50.0]]


def test_unresolvable_overlap_is_reported() -> None:
    size = np.array([[20.0, 10.0], [20.0, 10.0]])
    cands = [np.array([[100.0, 100.0]]), np.array([[105.0, 100.0]])]
    _, before, left = place_labels(size, cands, [-1, -1], np.zeros((0, 2, 2)), NO_DISCS, NO_RECTS, BOUNDS)
    assert before == 1 and left == 1
