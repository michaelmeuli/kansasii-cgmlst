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
import pandas as pd
from adjustText import adjust_text
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Wedge

from common import OUT, SPECIES_COLORS, SUBSETS

TITLES = {"complex": "M. kansasii complex", "kansasii": "M. kansasii (sensu stricto)",
          "persicum": "M. persicum", "pseudokansasii": "M. pseudokansasii"}
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
        G.add_node(i, n=len(g), species=sp.to_dict(), ref=bool(samples.loc[g, "is_reference"].any()))
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

    fig, ax = plt.subplots(figsize=(16, 13) if k > 60 else (12, 10))
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
    # labels: sample number(s) for small nodes, "n=" for big ones; references by accession
    def label(g: list[str]) -> str:
        names = [s.replace("Mkan329-", "") for s in g]
        return ", ".join(names) if len(g) <= 3 else f"{', '.join(names[:2])} +{len(g) - 2}"
    texts = [ax.text(pos[i][0], pos[i][1] + radius(G.nodes[i]["n"]) + 0.6 * r0, label(g), fontsize=6, ha="center", va="bottom", zorder=4)
             for i, g in enumerate(groups)]
    pad = 0.06 * span
    ax.set_xlim(xy[:, 0].min() - pad, xy[:, 0].max() + pad)
    ax.set_ylim(xy[:, 1].min() - pad, xy[:, 1].max() + pad)
    ax.axis("off")
    fig.canvas.draw()
    adjust_text(texts, x=xy[:, 0], y=xy[:, 1], objects=edge_txt, ax=ax, expand=(1.1, 1.3), iter_lim=300)
    present = sorted({s for _, d in G.nodes(data=True) for s in d["species"]})
    fig.legend(handles=[Line2D([], [], marker="o", ls="", mfc=SPECIES_COLORS.get(s, "#999999"), mec="black", markersize=9, label=s)
                       for s in present], title="species", loc="lower right", frameon=False)
    fig.suptitle(f"{TITLES[name]}: minimum spanning tree of cgMLST allele profiles\n"
                 f"{len(D)} genomes, {k} types; distance = allele differences, pairwise deletion of missing loci; "
                 f"{'linear' if args.linear else 'logarithmic'} edge scale; dashed ring = public reference genome", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    for ext in ("png", "pdf"):
        fig.savefig(MST / f"mst_cgmlst_{name}.{ext}", dpi=200)
    plt.close(fig)
    pd.DataFrame({"node": range(k), "n": sizes, "species": [",".join(G.nodes[i]["species"]) for i in range(k)],
                  "samples": [",".join(g) for g in groups]}).to_csv(MST / f"mst_nodes_{name}.tsv", sep="\t", index=False)
    print(f"{name}: {len(D)} genomes, {k} types, wrote mst_cgmlst_{name}.png/pdf "
          f"(min shared loci {int(shared.to_numpy().min())})")
