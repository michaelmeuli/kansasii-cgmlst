#!/bin/bash
# Prodigal training file from the M. kansasii type strain ATCC 12478 (chromosome only: step 00 removes the plasmid pMK12478), single-genome mode.
set -euo pipefail
source "$(dirname "$0")/env.sh"
rm -f "$WORK/ATCC12478.trn"  # prodigal -t reads an existing file instead of writing it
prodigal -i "$WORK/reference_genomes/GCF_000157895.3.fasta" -t "$WORK/ATCC12478.trn" -p single -o /dev/null -q
ls -l "$WORK/ATCC12478.trn"
chewBBACA.py --version 2>&1 | head -1 || true
