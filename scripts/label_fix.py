"""Place text labels so that none overlaps another label, a node disc, a tree edge that is not its own, the legend or leaves the
figure. Works in display pixels. Each label gets a list of candidate centres (first = default); labels in conflict are moved to
the cheapest candidate given all other labels, pass after pass, until nothing is left or nothing improves."""
from __future__ import annotations

import numpy as np
import numpy.typing as npt

Arr = npt.NDArray[np.float64]  # float arrays: (n, 2) centres, (n, 4) rects x0 y0 x1 y1, (m, 2, 2) segments, (k, 3) discs x y r
HARD = 1000.0  # cost of one overlap; distance penalties stay far below it


def rects_at(centres: Arr, size: Arr, margin: float) -> Arr:
    """(n, 4) rectangles of the given size centred on `centres`, grown by `margin` on every side."""
    half = size / 2 + margin
    out: Arr = np.hstack([centres - half, centres + half])
    return out


def rect_hits_rects(a: Arr, b: Arr) -> npt.NDArray[np.bool_]:
    """(len(a), len(b)) True where rectangles a[i] and b[j] overlap."""
    return (a[:, None, 0] < b[None, :, 2]) & (b[None, :, 0] < a[:, None, 2]) & \
           (a[:, None, 1] < b[None, :, 3]) & (b[None, :, 1] < a[:, None, 3])


def rect_hits_discs(rects: Arr, discs: Arr) -> npt.NDArray[np.bool_]:
    """(len(rects), len(discs)) True where the rectangle comes closer to the disc centre than its radius."""
    dx = np.maximum(np.maximum(rects[:, None, 0] - discs[None, :, 0], discs[None, :, 0] - rects[:, None, 2]), 0.0)
    dy = np.maximum(np.maximum(rects[:, None, 1] - discs[None, :, 1], discs[None, :, 1] - rects[:, None, 3]), 0.0)
    return np.hypot(dx, dy) < discs[None, :, 2]


def rect_hits_segments(rects: Arr, segs: Arr) -> npt.NDArray[np.bool_]:
    """(len(rects), len(segs)) True where the segment touches the rectangle (bounding boxes overlap and the line
    does not pass wholly to one side of the rectangle's corners)."""
    p, q = segs[:, 0, :], segs[:, 1, :]
    box = (np.minimum(p[:, 0], q[:, 0])[None] < rects[:, None, 2]) & (rects[:, None, 0] < np.maximum(p[:, 0], q[:, 0])[None]) & \
          (np.minimum(p[:, 1], q[:, 1])[None] < rects[:, None, 3]) & (rects[:, None, 1] < np.maximum(p[:, 1], q[:, 1])[None])
    d = q - p
    side = np.stack([d[None, :, 0] * (cy - p[None, :, 1]) - d[None, :, 1] * (cx - p[None, :, 0])
                     for cx, cy in ((rects[:, None, 0], rects[:, None, 1]), (rects[:, None, 2], rects[:, None, 1]),
                                    (rects[:, None, 0], rects[:, None, 3]), (rects[:, None, 2], rects[:, None, 3]))])
    straddle = (side.min(axis=0) <= 0) & (side.max(axis=0) >= 0)
    hit: npt.NDArray[np.bool_] = box & straddle
    return hit


def place_labels(size: Arr, cands: list[Arr], own_seg: list[int], segs: Arr, discs: Arr, fixed: Arr, bounds: Arr,
                 margin: float = 2.0, max_pass: int = 60, walk: int = 3000) -> tuple[Arr, int, int]:
    """Choose one candidate centre per label. `size` (n, 2) label width/height; `cands[i]` (c_i, 2) candidate centres, the
    first being the default and the starting point; `own_seg[i]` index of the edge label i belongs to (-1 none) which it may
    cover; `fixed` (r, 4) rectangles to keep clear (legend); `bounds` (4,) rectangle all labels must stay inside.
    Returns (centres, overlaps before, overlaps left)."""
    n = len(size)
    static: list[npt.NDArray[np.float64]] = []  # per candidate: HARD * obstacle overlaps + distance penalty
    nhit: list[npt.NDArray[np.int64]] = []
    for i in range(n):
        r = rects_at(cands[i], np.broadcast_to(size[i], (len(cands[i]), 2)), margin)
        seg_hit = rect_hits_segments(r, segs) if len(segs) else np.zeros((len(r), 0), dtype=bool)
        if own_seg[i] >= 0:
            seg_hit[:, own_seg[i]] = False
        hits = seg_hit.sum(axis=1) + rect_hits_discs(r, discs).sum(axis=1) + \
            (rect_hits_rects(r, fixed).sum(axis=1) if len(fixed) else 0) + \
            ((r[:, 0] < bounds[0]) | (r[:, 1] < bounds[1]) | (r[:, 2] > bounds[2]) | (r[:, 3] > bounds[3]))
        nhit.append(np.asarray(hits, dtype=np.int64))
        static.append(HARD * nhit[-1] + 0.01 * np.hypot(*(cands[i] - cands[i][0]).T) + 0.001 * np.arange(len(r)))
    cur = np.zeros(n, dtype=np.int64)

    def centres() -> Arr:
        return np.array([cands[i][cur[i]] for i in range(n)])

    def conflicts() -> npt.NDArray[np.int64]:
        """Per label: obstacle overlaps + overlaps with other labels."""
        rects = rects_at(centres(), size, margin)
        pair = rect_hits_rects(rects, rects)
        np.fill_diagonal(pair, False)
        out: npt.NDArray[np.int64] = np.array([nhit[i][cur[i]] for i in range(n)], dtype=np.int64) + pair.sum(axis=1)
        return out

    def total() -> int:
        rects = rects_at(centres(), size, margin)
        pair = np.triu(rect_hits_rects(rects, rects), 1)
        return int(sum(int(nhit[i][cur[i]]) for i in range(n)) + pair.sum())

    def cost_of(i: int) -> Arr:
        """Cost of every candidate of label i given where all other labels are."""
        others = np.delete(rects_at(centres(), size, margin), i, axis=0)
        r = rects_at(cands[i], np.broadcast_to(size[i], (len(cands[i]), 2)), margin)
        out: Arr = static[i] + HARD * rect_hits_rects(r, others).sum(axis=1)
        return out

    before = total()
    for _ in range(max_pass):  # descent: conflicted labels (worst first) take their cheapest candidate
        bad = conflicts()
        if not bad.any():
            break
        moved = False
        for i in np.argsort(-bad):
            if not bad[i]:
                break
            cost = cost_of(int(i))
            best = int(np.argmin(cost))
            if cost[best] < cost[cur[i]] - 1e-9:
                cur[i] = best
                moved = True
        if not moved:
            break
    if total():  # stuck: min-conflicts random walk (a conflicted label must move, even uphill); keep the best state seen
        rng = np.random.default_rng(0)
        best_cur, best_total = cur.copy(), total()
        for _ in range(walk):
            bad = conflicts()
            if not bad.any():
                break
            i = int(rng.choice(np.nonzero(bad)[0]))
            cost = cost_of(i) + rng.uniform(0.0, 0.5, len(cands[i]))
            cost[cur[i]] = np.inf
            cur[i] = int(np.argmin(cost))
            t = total()
            if t < best_total:
                best_cur, best_total = cur.copy(), t
        cur[:] = best_cur
    return centres(), before, total()
