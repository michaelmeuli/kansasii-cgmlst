"""Pairwise allele distances with pairwise deletion (as in Ridom SeqSphere): a locus counts for a pair only if both
genomes have an allele call; missing loci (LNF, NIPH, ASM, ... or 0 in chewBBACA's masked table) are ignored for that pair."""
from __future__ import annotations

import numpy as np
import numpy.typing as npt
import pandas as pd

IntArray = npt.NDArray[np.int64]


def encode(table: pd.DataFrame) -> IntArray:
    """Allele table (strings such as '12', 'INF-12', 'LNF') -> int matrix, 0 = missing (allele ids start at 1)."""
    s = table.astype(str).apply(lambda c: c.str.replace("INF-", "", regex=False).str.strip())
    out: IntArray = s.apply(pd.to_numeric, errors="coerce").fillna(0).clip(lower=0).to_numpy(dtype=np.int64)
    return out


def pairwise_deletion(a: IntArray) -> tuple[IntArray, IntArray]:
    """(distance, shared) matrices: differing alleles among loci called in both genomes, and the number of such loci."""
    n = a.shape[0]
    dist = np.zeros((n, n), dtype=np.int64)
    shared = np.zeros((n, n), dtype=np.int64)
    called = a > 0
    for i in range(n):
        both = called[i] & called
        dist[i] = (both & (a != a[i])).sum(axis=1)
        shared[i] = both.sum(axis=1)
    return dist, shared
