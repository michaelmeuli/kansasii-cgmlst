#!/bin/bash
# AlleleCall of all genomes against the schema (adds new alleles to the schema copy, so the seed schema stays intact).
set -euo pipefail
source "$(dirname "$0")/env.sh"
rm -rf "$WORK/allele_call" "$WORK/schema_called"
cp -r "$WORK/schema/schema_seed" "$WORK/schema_called"
chewBBACA.py AlleleCall -i "$WORK/genomes" -g "$WORK/schema_called" -o "$WORK/allele_call" --ptf "$WORK/ATCC12478.trn" --cpu "$CPU" --output-masked
# result tables only; cds_coordinates.tsv and results_contigsInfo.tsv (50 MB) stay in $WORK
publish allele_call results_alleles.tsv results_alleles_masked.tsv results_statistics.tsv loci_summary_stats.tsv paralogous_counts.tsv paralogous_loci.tsv
