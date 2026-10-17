# Sourced by the shell steps: env_cgmlst on PATH, output dirs, thread count.
export KANSASII_ROOT="${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}"
export PATH="$HOME/data/conda/envs/env_cgmlst/bin:$PATH"
OUT="$KANSASII_ROOT/output/cgmlst"
CPU="${SLURM_CPUS_PER_TASK:-8}"
ATCC="$KANSASII_ROOT/data/gtdb_genomes/Mycobacteriaceae/ncbi_dataset/data/GCF_000157895.3/GCF_000157895.3_ASM15789v2_genomic.fna"
CGMLST_THRESHOLD="${CGMLST_THRESHOLD:-0.95}"   # ExtractCgMLST: loci present in >= this fraction of genomes
PCT=$(awk "BEGIN{printf \"%d\", $CGMLST_THRESHOLD*100}")   # file-name suffix used by ExtractCgMLST (0.95 -> 95)
