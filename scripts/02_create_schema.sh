#!/bin/bash
# CreateSchema from the seed genomes (ATCC 12478 + one public representative per species) with the ATCC 12478 PTF.
set -euo pipefail
source "$(dirname "$0")/env.sh"
rm -rf "$WORK/schema"
chewBBACA.py CreateSchema -i "$WORK/seed_genomes" -o "$WORK/schema" --ptf "$WORK/ATCC12478.trn" --cpu "$CPU"
ls "$WORK/schema/schema_seed" | wc -l
