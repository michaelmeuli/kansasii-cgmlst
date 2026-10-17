#!/bin/bash
# SchemaEvaluator on the loci of the final cgMLST scheme: copy those loci's FASTAs from the called schema, then evaluate.
set -euo pipefail
source "$(dirname "$0")/env.sh"
T="$OUT/cgmlst/t$CGMLST_THRESHOLD"
test -s "$T/cgMLSTschema$PCT.txt"
rm -rf "$OUT/cgmlst_schema" "$OUT/schema_evaluator"; mkdir -p "$OUT/cgmlst_schema"
while read -r locus; do
  [ -z "$locus" ] && continue
  f="${locus%.fasta}.fasta"
  cp "$OUT/schema_called/$f" "$OUT/cgmlst_schema/$f"
done < "$T/cgMLSTschema$PCT.txt"
chewBBACA.py SchemaEvaluator -g "$OUT/cgmlst_schema" -o "$OUT/schema_evaluator" --cpu "$CPU" --ta 11 --loci-reports
