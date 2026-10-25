#!/bin/bash
# RemoveGenes: drop the paralogous loci flagged by AlleleCall (paralogous_counts.tsv has them in column 1; header skipped).
set -euo pipefail
source "$(dirname "$0")/env.sh"
RES="$WORK/allele_call"
test -s "$RES/results_alleles.tsv"
ls "$RES"
tail -n +2 "$RES/paralogous_counts.tsv" | cut -f1 | sort -u > "$WORK/paralogous_loci.txt"
echo "paralogous loci: $(wc -l < "$WORK/paralogous_loci.txt")"
rm -rf "$WORK/no_paralogs"; mkdir -p "$WORK/no_paralogs"
if [ -s "$WORK/paralogous_loci.txt" ]; then
  chewBBACA.py RemoveGenes -i "$RES/results_alleles.tsv" -g "$WORK/paralogous_loci.txt" -o "$WORK/no_paralogs/results_alleles_noparalogs.tsv"
else
  cp "$RES/results_alleles.tsv" "$WORK/no_paralogs/results_alleles_noparalogs.tsv"
fi
cp "$WORK/paralogous_loci.txt" "$OUTPUT/paralogous_loci.txt"
publish no_paralogs results_alleles_noparalogs.tsv
