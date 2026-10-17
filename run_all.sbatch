#!/bin/bash
#SBATCH --job-name=kansasii-cgmlst
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.log
# Whole pipeline. STEPS="01 02 03" runs only those steps (default: all).
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_cgmlst
cd /shares/sander.imm.uzh/MM/kansasii/repos/kansasii-cgmlst/scripts
for f in [0-9][0-9]_*; do
  n=${f%%_*}
  [[ " ${STEPS:-$n} " == *" $n "* ]] || continue
  echo "=== $f $(date +%T)"
  case "$f" in *.sh) bash "$f" ;; *.py) python "$f" ;; esac
done
