#!/usr/bin/env bash
# Whole pipeline (steps 00-10). STEPS="01 02 03" runs only those steps (default: all). Order after a new allele call:
# STEPS="00 01 02 03 04 05 06 07", then immensekansasii 06 (results table), then STEPS="08 09 10". Separate the steps with
# ":" when passing them through sbatch (--export splits at commas): sbatch --export=ALL,STEPS=08:09:10 submit_all.sbatch
#SBATCH --job-name=kansasii-cgmlst
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/runs/cgmlst/logs/%x-%j.log

set -euo pipefail
STEPS=${STEPS//:/ }
source "${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}/repos/kansasii-docs/templates/kansasii_paths.sh"
REPO=$REPOS/kansasii-cgmlst
mkdir -p "$RUNS/cgmlst/logs"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_cgmlst
cd "$REPO/scripts"
for f in [0-9][0-9]_*; do
  n=${f%%_*}
  [[ " ${STEPS:-$n} " == *" $n "* ]] || continue
  echo "=== $f $(date +%T)"
  case "$f" in *.sh) bash "$f" ;; *.py) python "$f" ;; esac
done

# provenance (see kansasii-docs/docs/PROVENANCE.md); the Ridom tarball is a loose file in output/
bash "$MANIFEST_SCRIPT" "$OUT_CGMLST/results" "$REPO" submit_all.sbatch "$SCREENING_MAP_RESULTS" "$ASSEMBLIES" "$GTDB_REPRESENTATIVES" "STEPS=${STEPS:-all}"
if [[ " ${STEPS:-10} " == *" 10 "* ]]; then
  MANIFEST_NAME=ridom_assemblies.MANIFEST.txt bash "$MANIFEST_SCRIPT" "$OUT_CGMLST" "$REPO" submit_all.sbatch "$SCREENING_MAP_RESULTS" "$ASSEMBLIES" "STEPS=${STEPS:-all}"
fi
