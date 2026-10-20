"""Ridom SeqSphere-style minimum spanning trees from the pairwise-deletion allele distances (one figure per SUBSET).
Identical profiles (distance 0) share a node; node area ~ genome count, nodes with mixed species are pies; edge labels =
allele differences; edge length ~ log(1 + differences) (Ridom's logarithmic scale). MST ties are broken in favour of
edges between nodes holding more genomes (Ridom prefers the better-connected, larger types)."""
from __future__ import annotations

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import numpy.typing as npt
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Wedge
from matplotlib.text import Text

from common import OUT, SPECIES_COLORS, SUBSETS
from label_fix import place_labels
from layout_fix import crossings, find_overlaps, resolve_overlaps, untangle

TITLES = {"complex": "M. kansasii complex", "kansasii": "M. kansasii (sensu stricto)",
          "persicum": "M. persicum", "pseudokansasii": "M. pseudokansasii",
          "persicum_attenuatum_ostraviense_innocens_gastri": "M. persicum, attenuatum, ostraviense, innocens and gastri"}
ap = argparse.ArgumentParser()
ap.add_argument("--subset", choices=list(SUBSETS), default=None, help="only this subset (default: all)")
ap.add_argument("--cluster-threshold", type=int, default=None,
                help="shade nodes joined by edges of <= N allele differences (Ridom cluster alert distance)")
ap.add_argument("--linear", action="store_true", help="linear instead of logarithmic edge lengths")
args = ap.parse_args()
MST = OUT / "mst"


def kruskal(d: np.ndarray, size: np.ndarray) -> list[tuple[int, int, float]]:
    """MST edges; equal distances are taken in order of decreasing combined genome count of the two nodes."""
    n = len(d)
    edges = sorted(((float(d[i, j]), -int(size[i] + size[j]), i, j) for i in range(n) for j in range(i + 1, n)))
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    tree = []
    for w, _, i, j in edges:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj
            tree.append((i, j, w))
    return tree


for name in [args.subset] if args.subset else SUBSETS:
    D = pd.read_csv(MST / f"allele_distance_{name}.csv", index_col=0)
    shared = pd.read_csv(MST / f"shared_loci_{name}.csv", index_col=0)
    samples = pd.read_csv(OUT / "samples.tsv", sep="\t").set_index("sample")
    # nodes = groups of genomes at distance 0 from the group's first member
    groups: list[list[str]] = []
    for s in D.index:
        for g in groups:
            if D.loc[s, g[0]] == 0:
                g.append(s)
                break
        else:
            groups.append([s])
    k = len(groups)
    # distance between nodes = smallest pairwise distance between their members (0 cannot occur between nodes)
    ND = np.zeros((k, k))
    for i in range(k):
        for j in range(i + 1, k):
            ND[i, j] = ND[j, i] = max(float(D.loc[groups[i], groups[j]].to_numpy().min()), 0.5)
    sizes = np.array([len(g) for g in groups])
    tree = kruskal(ND, sizes)
    G = nx.Graph()
    for i, g in enumerate(groups):
        sp = samples.loc[g, "species"].value_counts()
        G.add_node(i, n=len(g), species=sp.to_dict(), ref=bool(samples.loc[g, "is_reference"].any()),
                   ctrl=[str(x) for x in samples.loc[g, "label"].fillna("") if x])
    for i, j, w in tree:
        G.add_edge(i, j, weight=w)
    assert nx.is_connected(G) and G.number_of_edges() == k - 1
    def length(w: float) -> float:
        return w if args.linear else float(np.log1p(w))

    Gl = nx.Graph((u, v, {"weight": length(d["weight"])}) for u, v, d in G.edges(data=True))
    pos = nx.kamada_kawai_layout(Gl, dist=dict(nx.all_pairs_dijkstra_path_length(Gl)))
    xy = np.array([pos[i] for i in range(k)])
    span = float(np.ptp(xy, axis=0).max()) or 1.0
    r0 = 0.012 * span * (1.4 if k < 40 else 1.0)  # radius of a single-genome node (data units)

    def radius(n: int) -> float:
        return r0 * float(np.sqrt(n))

    # no node may sit on an edge it does not belong to (nor on another node): nudge such nodes aside
    edge_list = [(int(u), int(v)) for u, v in G.edges()]
    radii = np.array([radius(G.nodes[i]["n"]) for i in range(k)])
    before = sum(len(x) for x in find_overlaps(xy, edge_list, radii, 0.5 * r0))
    xy, left = resolve_overlaps(xy, edge_list, radii, 0.5 * r0)
    ncross = len(crossings(xy, edge_list))
    for _ in range(3):  # rotate/mirror subtrees until no two edges cross, then clean up node overlaps again
        if not ncross:
            break
        xy = untangle(xy, edge_list, radii, 0.5 * r0)
        print(f"DBG {name}: after untangle crossings {len(crossings(xy, edge_list))} overlaps {sum(len(x) for x in find_overlaps(xy, edge_list, radii, 0.5 * r0))}")
        xy, left = resolve_overlaps(xy, edge_list, radii, 0.5 * r0)
        ncross = len(crossings(xy, edge_list))
        print(f"DBG {name}: after resolve crossings {ncross} left {left}")
    pos = {i: xy[i] for i in range(k)}
    print(f"{name}: node overlaps {before} -> {left}")
    print(f"{name}: edge crossings {ncross}")
    if ncross:
        print(f"WARNING {name}: {ncross} edge crossings could not be resolved")
        print("  crossing edges (node ids):", [(edge_list[i], edge_list[j]) for i, j in crossings(xy, edge_list)])
    if left:
        print(f"WARNING {name}: {left} node/edge overlaps could not be resolved")

    def build(scale: float) -> tuple[Figure, Axes, list[Text], list[Text], int, int]:
        """Draw the tree on a figure `scale` times the base size and place all labels; returns overlaps before/after."""
        fig, ax = plt.subplots(figsize=tuple(scale * v for v in ((16, 13) if k > 60 else (12, 10))))
        ax.set_aspect("equal")
        if args.cluster_threshold is not None:
            for u, v, d in G.edges(data=True):
                if d["weight"] <= args.cluster_threshold:
                    ax.plot(*zip(pos[u], pos[v]), color="#cccccc", lw=22, solid_capstyle="round", alpha=0.6, zorder=0)
        for u, v, d in G.edges(data=True):
            ax.plot(*zip(pos[u], pos[v]), color="#555555", lw=1.0, zorder=1)
        edge_txt = []
        for u, v, d in G.edges(data=True):
            mx, my = (pos[u][0] + pos[v][0]) / 2, (pos[u][1] + pos[v][1]) / 2
            edge_txt.append(ax.text(mx, my, f"{int(round(d['weight']))}", fontsize=6.5, ha="center", va="center", zorder=3,
                                    bbox={"facecolor": "white", "edgecolor": "#999999", "boxstyle": "round,pad=0.15", "lw": 0.4}))
        for i in range(k):
            n, sp = G.nodes[i]["n"], G.nodes[i]["species"]
            x, y = pos[i]
            r = radius(n)
            if len(sp) == 1:
                ax.add_patch(Circle((x, y), r, fc=SPECIES_COLORS.get(next(iter(sp)), "#999999"), ec="black", lw=0.7, zorder=2))
            else:  # pie of species counts
                a0 = 0.0
                for s_name, c in sp.items():
                    a1 = a0 + 360.0 * c / n
                    ax.add_patch(Wedge((x, y), r, a0, a1, fc=SPECIES_COLORS.get(str(s_name), "#999999"), ec="black", lw=0.5, zorder=2))
                    a0 = a1
                ax.add_patch(Circle((x, y), r, fc="none", ec="black", lw=0.7, zorder=2))
            if G.nodes[i]["ref"]:  # public reference genome inside the node
                ax.add_patch(Circle((x, y), r * 1.25, fc="none", ec="black", lw=1.2, ls=(0, (2, 1)), zorder=2))
        # labels: control strains by their LABEL (bold, thick frame), others by sample number(s) ("+n" for big nodes)
        def label(i: int, g: list[str]) -> str:
            names = [s.replace("Mkan329-", "M") for s, c in zip(g, samples.loc[g, "label"].fillna("")) if not c]
            others = ", ".join(names) if len(names) <= 3 else f"{', '.join(names[:2])} +{len(names) - 2}"
            return "\n".join([*G.nodes[i]["ctrl"], *([others] if others else [])])
        texts = [ax.text(pos[i][0], pos[i][1] + radius(G.nodes[i]["n"]) + 0.6 * r0, label(i, g), fontsize=7 if G.nodes[i]["ctrl"] else 6,
                         fontweight="bold" if G.nodes[i]["ctrl"] else "normal", ha="center", va="bottom", zorder=4,
                         bbox={"facecolor": "white", "alpha": 0.95 if G.nodes[i]["ctrl"] else 0.85, "pad": 1.5 if G.nodes[i]["ctrl"] else 0.8,
                               "lw": 1.8 if G.nodes[i]["ctrl"] else 0.8,
                               "edgecolor": SPECIES_COLORS.get(next(iter(G.nodes[i]["species"])), "#999999")})
                 for i, g in enumerate(groups)]
        pad = 0.06 * span
        ax.set_xlim(xy[:, 0].min() - pad, xy[:, 0].max() + pad)
        ax.set_ylim(xy[:, 1].min() - pad, xy[:, 1].max() + pad)
        ax.axis("off")
        pad = 0.06 * span
        ax.set_xlim(xy[:, 0].min() - pad, xy[:, 0].max() + pad)
        ax.set_ylim(xy[:, 1].min() - pad, xy[:, 1].max() + pad)
        ax.axis("off")
        present = sorted({s for _, d in G.nodes(data=True) for s in d["species"]})
        legend = fig.legend(handles=[Line2D([], [], marker="o", ls="", mfc=SPECIES_COLORS.get(s, "#999999"), mec="black", markersize=9, label=s)
                                     for s in present]
                            + ([Line2D([], [], marker="*", ls="", mfc="none", mec="black", markersize=12, label="control strain (ATCC/DSM)")]
                               if any(d["ctrl"] for _, d in G.nodes(data=True)) else []),
                            title="species", loc="lower right", frameon=False)
        fig.suptitle(f"{TITLES[name]}: minimum spanning tree of cgMLST allele profiles\n"
                     f"{len(D)} genomes, {k} types; distance = allele differences, pairwise deletion of missing loci; "
                     f"{'linear' if args.linear else 'logarithmic'} edge scale; dashed ring = public reference genome", fontsize=10)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        fig.canvas.draw()
        # labels are placed in display pixels: none may touch another label, a node (incl. its reference/control ring),
        # an edge that is not its own, or the legend
        for t in [*texts, *edge_txt]:
            t.set_horizontalalignment("center")
            t.set_verticalalignment("center")
        px = ax.transData.transform(xy)
        sc = float(ax.transData.transform((1.0, 0.0))[0] - ax.transData.transform((0.0, 0.0))[0])
        ring = np.array([1.7 if G.nodes[i]["ctrl"] else 1.25 if G.nodes[i]["ref"] else 1.0 for i in range(k)])
        rad_px = radii * sc * ring
        discs = np.column_stack([px, rad_px])
        segs = np.array([[px[u], px[v]] for u, v in edge_list]).reshape(-1, 2, 2)
        sizes = np.array([[b.width, b.height] for b in (t.get_window_extent() for t in [*texts, *edge_txt])])
        cands: list[npt.NDArray[np.float64]] = []
        for i in range(k):
            w, h = sizes[i]
            default = px[i] + np.array([0.0, rad_px[i] + 0.6 * r0 * sc + h / 2])
            ang = np.linspace(0, 2 * np.pi, 48, endpoint=False)
            unit = np.column_stack([np.cos(ang), np.sin(ang)])
            reach = np.abs(unit[:, 0]) * w / 2 + np.abs(unit[:, 1]) * h / 2
            ring_c = [px[i] + unit * (rad_px[i] + gap + reach)[:, None] for gap in (3.0, 8.0, 16.0, 30.0, 50.0, 80.0, 120.0, 200.0, 320.0)]
            cands.append(np.vstack([default, *ring_c]))
        ts = np.concatenate([[0.5], np.linspace(0.15, 0.85, 15)])
        for e, (u, v) in enumerate(edge_list):  # on the edge first, then beside it (either side, further out if crowded; joined by a leader line)
            w, h = sizes[k + e]
            along = px[u] + ts[:, None] * (px[v] - px[u])
            normal = np.array([-(px[v][1] - px[u][1]), px[v][0] - px[u][0]]) / max(float(np.hypot(*(px[v] - px[u]))), 1e-9)
            off = abs(normal[0]) * w / 2 + abs(normal[1]) * h / 2 + 3.0
            beside = [along[::2] + sign * lvl * off * normal for lvl in (1.0, 2.5, 4.0, 6.0, 9.0) for sign in (1.0, -1.0)]
            cands.append(np.vstack([along, *beside]))
        wf, hf = fig.bbox.width, fig.bbox.height
        placed, before, left = place_labels(sizes, cands, [-1] * k + list(range(len(edge_list))), segs, discs,
                                            np.array([legend.get_window_extent().extents]), np.array([0.0, 0.0, wf, 0.95 * hf]))
        for t, c in zip([*texts, *edge_txt], placed):
            t.set_position(tuple(ax.transData.inverted().transform(c)))
        return fig, ax, texts, edge_txt, before, left

    for attempt in range(7):  # not enough room: enlarge the figure (label sizes are fixed in points)
        fig, ax, texts, edge_txt, before, left = build(1.1 ** attempt)
        if left == 0:
            break
        plt.close(fig)
    print(f"{name}: label overlaps {before} -> {left} (figure scale {1.1 ** attempt:.2f})")
    if left:
        print(f"WARNING {name}: {left} label overlaps could not be resolved")
    # leader lines: join every label that is not touching its node to the node edge (labels are moved by label_fix.place_labels)
    fig.canvas.draw()
    inv = ax.transData.inverted()
    for i, txt in enumerate(texts):
        bb = txt.get_window_extent()
        cen = np.asarray(ax.transData.transform(pos[i]), dtype=float)
        rad = float(np.hypot(*(ax.transData.transform((pos[i][0] + radius(G.nodes[i]["n"]), pos[i][1])) - cen)))
        nearest = np.array([min(max(cen[0], bb.x0), bb.x1), min(max(cen[1], bb.y0), bb.y1)])
        gap = float(np.hypot(*(nearest - cen)))
        if G.nodes[i]["ctrl"]:  # star ring around control strains (size from the node radius in pixels)
            ax.scatter([pos[i][0]], [pos[i][1]], marker="*", s=(3.2 * rad * 72 / fig.dpi) ** 2, facecolors="none", edgecolors="black",
                       linewidths=1.3, zorder=2.6)
        if gap > rad + 3:
            start = cen + (nearest - cen) / gap * rad
            (x0, y0), (x1, y1) = inv.transform([tuple(start), tuple(nearest)])
            ax.add_line(Line2D([x0, x1], [y0, y1], color=SPECIES_COLORS.get(next(iter(G.nodes[i]["species"])), "#444444"),
                               lw=0.9, ls=(0, (1.5, 1.5)), zorder=3.5))  # dotted, so it is not mistaken for a tree edge
    for e, (u, v) in enumerate(edge_list):  # edge labels moved off their edge: dotted leader to the edge
        bb = edge_txt[e].get_window_extent()
        a, b = (np.asarray(ax.transData.transform(pos[n]), dtype=float) for n in (u, v))
        c = np.array([(bb.x0 + bb.x1) / 2, (bb.y0 + bb.y1) / 2])
        t = float(np.clip(np.dot(c - a, b - a) / max(float(np.dot(b - a, b - a)), 1e-12), 0.0, 1.0))
        on_edge = a + t * (b - a)
        nearest = np.array([min(max(on_edge[0], bb.x0), bb.x1), min(max(on_edge[1], bb.y0), bb.y1)])
        if float(np.hypot(*(nearest - on_edge))) > 3:
            (x0, y0), (x1, y1) = inv.transform([tuple(on_edge), tuple(nearest)])
            ax.add_line(Line2D([x0, x1], [y0, y1], color="#999999", lw=0.7, ls=(0, (1.5, 1.5)), zorder=3.5))
    for ext in ("png", "pdf"):
        fig.savefig(MST / f"mst_cgmlst_{name}.{ext}", dpi=200)
    plt.close(fig)
    pd.DataFrame({"node": range(k), "n": sizes, "species": [",".join(G.nodes[i]["species"]) for i in range(k)],
                  "control": ["; ".join(G.nodes[i]["ctrl"]) for i in range(k)],
                  "samples": [",".join(g) for g in groups]}).to_csv(MST / f"mst_nodes_{name}.tsv", sep="\t", index=False)
    print(f"{name}: {len(D)} genomes, {k} types, wrote mst_cgmlst_{name}.png/pdf "
          f"(min shared loci {int(shared.to_numpy().min())})")
