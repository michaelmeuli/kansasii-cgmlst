#!/bin/bash
# ExtractCgMLST at 100 %, 99 % and CGMLST_THRESHOLD (default 95 %); the last one is the scheme used downstream.
set -euo pipefail
source "$(dirname "$0")/env.sh"
IN="$WORK/no_paralogs/results_alleles_noparalogs.tsv"
test -s "$IN"
rm -rf "$WORK/cgmlst"; mkdir -p "$WORK/cgmlst"
for t in 1.0 0.99 "$CGMLST_THRESHOLD"; do
  chewBBACA.py ExtractCgMLST -i "$IN" -o "$WORK/cgmlst/t$t" --t "$t"
  p=$(awk "BEGIN{printf \"%d\", $t*100}")
  publish "cgmlst/t$t" "cgMLST$p.tsv" "cgMLSTschema$p.txt" missing_loci_stats.tsv presence_absence.tsv cgMLST.html
done
for t in 1.0 0.99 "$CGMLST_THRESHOLD"; do p=$(awk "BEGIN{printf \"%d\", $t*100}"); echo "t=$t loci=$(wc -l < "$WORK/cgmlst/t$t/cgMLSTschema$p.txt")"; done
