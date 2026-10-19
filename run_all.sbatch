#!/bin/bash
#SBATCH --job-name=kansasii-cgmlst
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/cgmlst/logs/%x-%j.log
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

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash /shares/sander.imm.uzh/MM/kansasii/repos/kansasii-docs/bin/write_manifest.sh \
  "/shares/sander.imm.uzh/MM/kansasii/output/cgmlst" \
  /shares/sander.imm.uzh/MM/kansasii/repos/kansasii-cgmlst run_all.sbatch /shares/sander.imm.uzh/MM/kansasii/data/imm/screening_map_results.csv
