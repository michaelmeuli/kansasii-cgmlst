"""Remove node/edge overlaps from a tree layout: no node disc may touch an edge it is not an endpoint of, and no two discs
may overlap. Overlapping nodes are nudged perpendicular to the edge (or apart) until the layout is clean."""
from __future__ import annotations

import numpy as np
import numpy.typing as npt

Pos = npt.NDArray[np.float64]  # (n, 2) node centres
Edges = list[tuple[int, int]]


def _closest(p: Pos, a: Pos, b: Pos) -> Pos:
    """Point of segment ab closest to p."""
    ab = b - a
    t = float(np.clip(np.dot(p - a, ab) / max(float(np.dot(ab, ab)), 1e-12), 0.0, 1.0))
    out: Pos = a + t * ab
    return out


def find_overlaps(pos: Pos, edges: Edges, radii: npt.NDArray[np.float64], margin: float) -> tuple[list[tuple[int, int]],
                                                                                         list[tuple[int, int]]]:
    """(node, edge index) pairs of nodes on non-incident edges, and (i, j) pairs of overlapping nodes."""
    n = len(pos)
    on_edge: list[tuple[int, int]] = []
    if edges:
        a = pos[[u for u, _ in edges]]  # (m, 2)
        ab = pos[[v for _, v in edges]] - a
        denom = np.maximum((ab * ab).sum(axis=1), 1e-12)
        rel = pos[:, None, :] - a[None, :, :]  # (n, m, 2)
        t = np.clip((rel * ab[None]).sum(axis=2) / denom[None], 0.0, 1.0)
        dist = np.hypot(*(rel - t[..., None] * ab[None]).transpose(2, 0, 1))  # (n, m)
        ends = np.zeros((n, len(edges)), dtype=bool)
        for ei, (u, v) in enumerate(edges):
            ends[u, ei] = ends[v, ei] = True
        hit = (dist < radii[:, None] + margin) & ~ends
        on_edge = [(int(i), int(ei)) for i, ei in zip(*np.nonzero(hit))]
    dd = np.hypot(*(pos[:, None, :] - pos[None, :, :]).transpose(2, 0, 1))
    close = np.triu(dd < radii[:, None] + radii[None, :] + margin, 1)
    touching = [(int(i), int(j)) for i, j in zip(*np.nonzero(close))]
    return on_edge, touching


def resolve_overlaps(pos: Pos, edges: Edges, radii: npt.NDArray[np.float64], margin: float,
                     max_iter: int = 2000) -> tuple[Pos, int]:
    """Nudge nodes off non-incident edges and off each other; returns (new positions, overlaps left)."""
    pos = pos.copy()
    for _ in range(max_iter):
        on_edge, touching = find_overlaps(pos, edges, radii, margin)
        if not on_edge and not touching:
            return pos, 0
        for n, ei in on_edge:
            u, v = edges[ei]
            c = _closest(pos[n], pos[u], pos[v])
            d = pos[n] - c
            dist = float(np.hypot(*d))
            if dist < 1e-9:  # centred exactly on the edge: step to the side of the edge
                ab = pos[v] - pos[u]
                d = np.array([-ab[1], ab[0]]) / max(float(np.hypot(*ab)), 1e-12)
                dist = 0.0
            else:
                d = d / dist
            pos[n] = pos[n] + d * (radii[n] + margin - dist + 1e-6)
        for i, j in touching:
            d = pos[i] - pos[j]
            dist = float(np.hypot(*d))
            d = d / dist if dist > 1e-9 else np.array([1.0, 0.0])
            push = (radii[i] + radii[j] + margin - dist) / 2 + 1e-6
            pos[i] = pos[i] + d * push
            pos[j] = pos[j] - d * push
    on_edge, touching = find_overlaps(pos, edges, radii, margin)
    return pos, len(on_edge) + len(touching)
