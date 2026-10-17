#!/bin/bash
# CreateSchema from the seed genomes (ATCC 12478 + one public representative per species) with the ATCC 12478 PTF.
set -euo pipefail
source "$(dirname "$0")/env.sh"
rm -rf "$OUT/schema"
chewBBACA.py CreateSchema -i "$OUT/seed_genomes" -o "$OUT/schema" --ptf "$OUT/ATCC12478.trn" --cpu "$CPU"
ls "$OUT/schema/schema_seed" | wc -l
