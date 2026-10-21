"""Pairwise-deletion allele distance matrices of the final cgMLST table, for the whole table and per SUBSET.
Writes allele_distance_<subset>.csv and shared_loci_<subset>.csv (loci called in both genomes of a pair) to output/cgmlst/mst/."""
from __future__ import annotations

import argparse
import pandas as pd

from allele_distance import encode, pairwise_deletion
from common import CGMLST, OUT, REFERENCES, SUBSETS, load_samples

ap = argparse.ArgumentParser()
ap.add_argument("--threshold", default="0.95", help="cgMLST threshold directory to use (t<threshold>)")
args = ap.parse_args()

tsv = CGMLST / f"t{args.threshold}" / f"cgMLST{round(float(args.threshold) * 100)}.tsv"
table = pd.read_csv(tsv, sep="\t", index_col=0, dtype=str)
table.index = pd.Index([str(i).removesuffix(".fasta") for i in table.index], name="sample")  # file names -> sample ids
samples = pd.read_csv(OUT / "samples.tsv", sep="\t").set_index("sample")
excluded = set(load_samples().query("excluded")["sample"])  # screening_map_results.csv: contaminated or < 90 % loci called
samples = samples[~samples.index.isin(excluded)]
print(f"excluded from distances (`excluded` column): {sorted(excluded & set(table.index))}")
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

# sanity check: each ATCC/DSM control strain vs the public genome of the same strain/type strain (expected ~0 alleles)
ref_of = {sp: acc for acc, sp in REFERENCES.items()}
dist_all = pd.read_csv(mst_dir / "allele_distance_complex.csv", index_col=0)
rows = [{"control": str(s), "label": r.label, "public_reference": ref_of[str(r.species)],
         "allele_differences": int(dist_all.to_numpy()[dist_all.index.get_loc(str(s)), dist_all.columns.get_loc(ref_of[str(r.species)])])}
        for s, r in samples[samples.is_control].iterrows() if r.species in ref_of]
chk = pd.DataFrame(rows)
chk.to_csv(mst_dir / "control_check.tsv", sep="\t", index=False)
print(chk.to_string(index=False))
