"""Pairwise-deletion allele distance matrices of the final cgMLST table, for the whole table and per SUBSET.
Writes allele_distance_<subset>.csv and shared_loci_<subset>.csv (loci called in both genomes of a pair) to output/cgmlst/mst/."""
from __future__ import annotations

import argparse
import pandas as pd

from allele_distance import encode, pairwise_deletion
from common import CGMLST, OUT, SUBSETS

ap = argparse.ArgumentParser()
ap.add_argument("--threshold", default="0.95", help="cgMLST threshold directory to use (t<threshold>)")
args = ap.parse_args()

tsv = CGMLST / f"t{args.threshold}" / f"cgMLST{round(float(args.threshold) * 100)}.tsv"
table = pd.read_csv(tsv, sep="\t", index_col=0, dtype=str)
table.index = pd.Index([str(i).removesuffix(".fasta") for i in table.index], name="sample")  # file names -> sample ids
samples = pd.read_csv(OUT / "samples.tsv", sep="\t").set_index("sample")
missing = set(samples.index) - set(table.index)
if missing:  # genomes dropped by ExtractCgMLST (too many missing loci) -- report, do not hide
    print(f"WARNING {len(missing)} genomes not in the cgMLST table: {sorted(missing)}")
samples = samples.loc[[s for s in samples.index if s in table.index]]
table = table.loc[samples.index]
mst_dir = OUT / "mst"
mst_dir.mkdir(exist_ok=True)
for name, species in SUBSETS.items():
    keep = samples.index if species is None else samples.index[samples.species.isin(species)]
    A = encode(table.loc[keep])
    dist, shared = pairwise_deletion(A)
    pd.DataFrame(dist, index=keep, columns=keep).to_csv(mst_dir / f"allele_distance_{name}.csv")
    pd.DataFrame(shared, index=keep, columns=keep).to_csv(mst_dir / f"shared_loci_{name}.csv")
    print(f"{name}: {len(keep)} genomes x {A.shape[1]} loci; min shared loci per pair {shared.min()}, "
          f"median distance {pd.Series(dist[dist > 0]).median():.0f}")
