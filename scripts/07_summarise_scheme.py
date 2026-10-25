"""Percentage of good targets from the SchemaEvaluator report (data embedded in schema_report.html) and the cgMLST size per
presence threshold. A locus is a good target when none of its alleles is flagged (no invalid alleles: incomplete ORF,
ambiguous bases, missing start/stop, in-frame stop; none shorter than 201 bp or outside the size threshold)."""
from __future__ import annotations

import json
import re

import pandas as pd

from allele_distance import encode
from common import CGMLST, EVAL, OUTPUT

html = (EVAL / "schema_report.html").read_text()
m = re.search(r"preComputedData\s*=\s*(\{.*?\})\s*;?\s*</script>", html, flags=re.S)
assert m, "preComputedData not found in schema_report.html"
data = json.loads(m.group(1))
cols = data["analysis"][0]["columns"]
loci = pd.DataFrame(data["analysis"][1]["rows"], columns=cols)
flag_cols = ["Invalid Alleles", "Alleles < 201bp", "Alleles below threshold", "Alleles above threshold"]
flag_cols = [c for c in cols if c in flag_cols or c.startswith("Alleles <") or c.startswith("Alleles above")]
flagged = loci[flag_cols].apply(pd.to_numeric, errors="coerce").fillna(0).sum(axis=1) > 0
loci["good"] = ~flagged
loci.to_csv(OUTPUT / "schema_evaluator_loci.tsv", sep="\t", index=False)
n = len(loci)
rows = [("loci evaluated", n), ("good targets (no flagged allele)", int((~flagged).sum())),
        ("good targets (%)", round(100 * (~flagged).mean(), 2))]
valid = pd.to_numeric(loci["Proportion of Validated Alleles"], errors="coerce")
rows.append(("loci with all alleles valid (%)", round(100 * float((valid == 1).mean()), 2)))
for c in flag_cols:
    rows.append((f"loci with >=1 allele flagged: {c}", int((pd.to_numeric(loci[c], errors="coerce") > 0).sum())))
for t in ("1.0", "0.99", "0.95"):
    p = round(float(t) * 100)
    f = CGMLST / f"t{t}" / f"cgMLSTschema{p}.txt"
    if f.exists():
        rows.append((f"cgMLST loci at >= {p}% presence", sum(1 for _ in f.open(encoding="utf-8"))))
out = pd.DataFrame(rows, columns=["metric", "value"])
out.to_csv(OUTPUT / "scheme_summary.tsv", sep="\t", index=False)
print(out.to_string(index=False))

# per-genome estimate for screening_map_results.csv (immensekansasii 06): % of the cgMLST loci with an allele call
tsv = CGMLST / "t0.95" / "cgMLST95.tsv"
if tsv.exists():
    table = pd.read_csv(tsv, sep="\t", index_col=0, dtype=str)
    called = (encode(table) > 0).mean(axis=1) * 100
    qc = pd.DataFrame({"sample": [str(i).removesuffix(".fasta") for i in table.index],
                       "cgmlst_called_pct": called.round(2)}).sort_values("cgmlst_called_pct")
    qc.to_csv(OUTPUT / "sample_qc.tsv", sep="\t", index=False)
    print(f"\nsample_qc.tsv: {len(qc)} genomes x {table.shape[1]} loci; lowest:\n{qc.head(5).to_string(index=False)}")
