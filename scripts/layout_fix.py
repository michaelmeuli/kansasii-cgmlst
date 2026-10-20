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
                     max_iter: int = 2000, max_scale: int = 25) -> tuple[Pos, int]:
    """Nudge nodes off non-incident edges and off each other; returns (new positions, overlaps left). If nudging alone
    gets stuck in a crowded area, the whole layout is enlarged by 4 % about its centroid (up to max_scale times)."""
    pos = pos.copy()
    for _ in range(max_scale + 1):
        pos, left = _nudge(pos, edges, radii, margin, max_iter)
        if left == 0:
            return pos, 0
        pos = pos.mean(axis=0) + (pos - pos.mean(axis=0)) * 1.04
    return pos, left


def _nudge(pos: Pos, edges: Edges, radii: npt.NDArray[np.float64], margin: float, max_iter: int) -> tuple[Pos, int]:
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


def crossings(pos: Pos, edges: Edges) -> list[tuple[int, int]]:
    """Index pairs of edges that cross each other (edges sharing a node are not counted)."""
    m = len(edges)
    if m < 2:
        return []
    a = pos[[u for u, _ in edges]]
    b = pos[[v for _, v in edges]]
    d = b - a

    def side(p: Pos) -> npt.NDArray[np.float64]:  # (m, m): side of point p[j] relative to the line of edge i
        out: npt.NDArray[np.float64] = d[:, None, 0] * (p[None, :, 1] - a[:, None, 1]) - d[:, None, 1] * (p[None, :, 0] - a[:, None, 0])
        return out

    s1, s2 = side(a), side(b)  # edge j's end points vs edge i's line
    hit = (s1 * s2 < 0) & (s1.T * s2.T < 0)
    ends = np.array(edges)
    share = (ends[:, None, :, None] == ends[None, :, None, :]).any(axis=(2, 3))
    hit &= ~share
    return [(int(i), int(j)) for i, j in zip(*np.nonzero(np.triu(hit, 1)))]


def untangle(pos: Pos, edges: Edges, radii: npt.NDArray[np.float64], margin: float, max_iter: int = 50) -> Pos:
    """Remove edge crossings of a tree layout by rotating (or mirroring) the subtree hanging off one end of a crossing edge about
    its other end; edge lengths are kept. A move is taken when it lowers the number of crossings, ties broken by fewer node overlaps."""
    pos = pos.copy()
    nbr: dict[int, list[int]] = {}
    for u, v in edges:
        nbr.setdefault(u, []).append(v)
        nbr.setdefault(v, []).append(u)

    def subtree(p: int, c: int) -> list[int]:
        seen, stack = {p, c}, [c]
        out = [c]
        while stack:
            for w in nbr[stack.pop()]:
                if w not in seen:
                    seen.add(w)
                    out.append(w)
                    stack.append(w)
        return out

    def score(q: Pos) -> tuple[int, int]:
        return len(crossings(q, edges)), sum(len(x) for x in find_overlaps(q, edges, radii, margin))

    best = score(pos)
    for _ in range(max_iter):
        cross = crossings(pos, edges)
        if not cross:
            break
        improved = False
        for ei in dict.fromkeys(e for pair in cross for e in pair):
            for p, c in (edges[ei], edges[ei][::-1]):
                nodes = subtree(p, c)
                rel = pos[nodes] - pos[p]
                axis = np.arctan2(*(pos[c] - pos[p])[::-1])
                for mirror in (False, True):
                    for deg in range(5, 360, 5):
                        th = np.radians(deg)
                        r = rel
                        if mirror:  # reflect about the line p-c, then rotate
                            r = np.stack([np.cos(2 * axis) * rel[:, 0] + np.sin(2 * axis) * rel[:, 1],
                                          np.sin(2 * axis) * rel[:, 0] - np.cos(2 * axis) * rel[:, 1]], axis=1)
                        rot = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
                        cand = pos.copy()
                        cand[nodes] = pos[p] + r @ rot.T
                        s = score(cand)
                        if s < best:
                            best, pos, improved = s, cand, True
                            rel = pos[nodes] - pos[p]
                            axis = np.arctan2(*(pos[c] - pos[p])[::-1])
        if not improved:
            break
    return pos
