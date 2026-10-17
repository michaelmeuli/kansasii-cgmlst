#!/bin/bash
# RemoveGenes: drop the paralogous loci flagged by AlleleCall (paralogous_counts.tsv has them in column 1; header skipped).
set -euo pipefail
source "$(dirname "$0")/env.sh"
RES="$OUT/allele_call"
test -s "$RES/results_alleles.tsv"
ls "$RES"
tail -n +2 "$RES/paralogous_counts.tsv" | cut -f1 | sort -u > "$OUT/paralogous_loci.txt"
echo "paralogous loci: $(wc -l < "$OUT/paralogous_loci.txt")"
rm -rf "$OUT/no_paralogs"; mkdir -p "$OUT/no_paralogs"
if [ -s "$OUT/paralogous_loci.txt" ]; then
  chewBBACA.py RemoveGenes -i "$RES/results_alleles.tsv" -g "$OUT/paralogous_loci.txt" -o "$OUT/no_paralogs/results_alleles_noparalogs.tsv"
else
  cp "$RES/results_alleles.tsv" "$OUT/no_paralogs/results_alleles_noparalogs.tsv"
fi
