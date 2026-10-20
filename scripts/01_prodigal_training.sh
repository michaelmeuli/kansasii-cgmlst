#!/bin/bash
# Prodigal training file from the M. kansasii type strain ATCC 12478 (chromosome + plasmid), single-genome mode.
set -euo pipefail
source "$(dirname "$0")/env.sh"
prodigal -i "$ATCC" -t "$WORK/ATCC12478.trn" -p single -o /dev/null -q
ls -l "$WORK/ATCC12478.trn"
chewBBACA.py --version 2>&1 | head -1 || true
