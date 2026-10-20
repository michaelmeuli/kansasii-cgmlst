#!/bin/bash
# SchemaEvaluator on the loci of the final cgMLST scheme: copy those loci's FASTAs from the called schema, then evaluate.
set -euo pipefail
source "$(dirname "$0")/env.sh"
T="$WORK/cgmlst/t$CGMLST_THRESHOLD"
test -s "$T/cgMLSTschema$PCT.txt"
rm -rf "$WORK/cgmlst_schema" "$WORK/schema_evaluator"; mkdir -p "$WORK/cgmlst_schema"
while read -r locus; do
  [ -z "$locus" ] && continue
  f="${locus%.fasta}.fasta"
  cp "$WORK/schema_called/$f" "$WORK/cgmlst_schema/$f"
done < "$T/cgMLSTschema$PCT.txt"
chewBBACA.py SchemaEvaluator -g "$WORK/cgmlst_schema" -o "$WORK/schema_evaluator" --cpu "$CPU" --ta 11 --loci-reports
publish schema_evaluator schema_report.html
