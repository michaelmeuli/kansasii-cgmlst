"""Hand-computed checks of the pairwise-deletion distance."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from allele_distance import encode, pairwise_deletion  # noqa: E402


def test_pairwise_deletion() -> None:
    t = pd.DataFrame([["1", "1", "1", "1"],
                      ["1", "2", "LNF", "1"],     # differs at locus 2; locus 3 missing -> ignored
                      ["INF-1", "NIPH", "2", "0"]],  # only loci 1 and 3 comparable with row 0, loci 1 with row 1
                     index=["a", "b", "c"])
    a = encode(t)
    assert a.tolist() == [[1, 1, 1, 1], [1, 2, 0, 1], [1, 0, 2, 0]]
    d, s = pairwise_deletion(a)
    assert d[0, 1] == 1 and s[0, 1] == 3
    assert d[0, 2] == 1 and s[0, 2] == 2  # loci 1 (same) and 3 (differs)
    assert d[1, 2] == 0 and s[1, 2] == 1  # only locus 1 is called in both
    assert (d == d.T).all() and (np.diag(d) == 0).all()
