#!/bin/bash
# ExtractCgMLST at 100 %, 99 % and CGMLST_THRESHOLD (default 95 %); the last one is the scheme used downstream.
set -euo pipefail
source "$(dirname "$0")/env.sh"
IN="$OUT/no_paralogs/results_alleles_noparalogs.tsv"
test -s "$IN"
rm -rf "$OUT/cgmlst"; mkdir -p "$OUT/cgmlst"
for t in 1.0 0.99 "$CGMLST_THRESHOLD"; do
  chewBBACA.py ExtractCgMLST -i "$IN" -o "$OUT/cgmlst/t$t" --t "$t"
done
for t in 1.0 0.99 "$CGMLST_THRESHOLD"; do p=$(awk "BEGIN{printf \"%d\", $t*100}"); echo "t=$t loci=$(wc -l < "$OUT/cgmlst/t$t/cgMLSTschema$p.txt")"; done
