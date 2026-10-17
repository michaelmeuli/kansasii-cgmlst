#!/bin/bash
# AlleleCall of all genomes against the schema (adds new alleles to the schema copy, so the seed schema stays intact).
set -euo pipefail
source "$(dirname "$0")/env.sh"
rm -rf "$OUT/allele_call"
cp -r "$OUT/schema/schema_seed" "$OUT/schema_called"
chewBBACA.py AlleleCall -i "$OUT/genomes" -g "$OUT/schema_called" -o "$OUT/allele_call" --ptf "$OUT/ATCC12478.trn" --cpu "$CPU" --output-masked
