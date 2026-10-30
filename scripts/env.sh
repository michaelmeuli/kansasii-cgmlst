# Sourced by the shell steps: shared project paths, env_cgmlst on PATH, output dirs, thread count.
source "${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}/repos/kansasii-docs/templates/kansasii_paths.sh"
export KANSASII_ROOT=$ROOT
export PATH="$HOME/data/conda/envs/env_cgmlst/bin:$PATH"
OUTPUT=$OUT_CGMLST/results   # final tables and figures only
WORK=$RUNS/cgmlst            # chewBBACA working state
mkdir -p "$OUTPUT" "$WORK"
# publish DIR FILE...: copy the result tables FILE (relative to $WORK/DIR) into $OUTPUT/DIR
publish() { local d="$1"; shift; mkdir -p "$OUTPUT/$d"; for f in "$@"; do cp "$WORK/$d/$f" "$OUTPUT/$d/$f"; done; }
CPU="${SLURM_CPUS_PER_TASK:-8}"
CGMLST_THRESHOLD="${CGMLST_THRESHOLD:-0.95}"   # ExtractCgMLST: loci present in >= this fraction of genomes
PCT=$(awk "BEGIN{printf \"%d\", $CGMLST_THRESHOLD*100}")   # file-name suffix used by ExtractCgMLST (0.95 -> 95)
