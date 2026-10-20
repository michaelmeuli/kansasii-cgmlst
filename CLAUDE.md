# kansasii-cgmlst

chewBBACA cgMLST scheme + Ridom-style MSTs of the M. kansasii complex. Only what is easy to get wrong:

- chewBBACA 3.5.4 uses Pyrodigal for prediction; the `.trn` training file comes from the Prodigal binary (in `env_cgmlst`).
- Two trees: `runs/cgmlst/` (`WORK`: everything chewBBACA writes, steps read from here) and `output/cgmlst/` (`OUT`: result tables/figures only, filled by `publish`). Do not point a step's chewBBACA output at `OUT`.
- Output names differ from older docs: `cgMLST<pct>.tsv` / `cgMLSTschema<pct>.txt` (95, 99, 100), flat `allele_call/` dir, `SchemaEvaluator -g`.
- AlleleCall adds alleles to the schema, so step 03 runs on a copy (`schema_called`); step 06 reads loci from that copy.
- Genome ids with dots (GCF_000157895.3): never `Path.stem`; strip `.fasta` only.
- Distances use pairwise deletion (`allele_distance.py`); 0 in the cgMLST table = missing, allele ids start at 1.
- Good-target % is parsed from the JSON embedded in `schema_report.html` (07), not a separate table.
- `run_all.sbatch` uses `python`, not `python -I` (scripts import `common`). Run heavy steps via sbatch, not the login node.
- Sample QC is copied from `kansasii-mst/scripts/common.py` but differs on purpose: control strains (`LABEL`) of complex species (NR 128/132/133) are kept and shown by their LABEL; other rules must stay in sync.
- `layout_fix.py` guarantees no node on a non-incident edge (enlarges the layout if nudging gets stuck); `09` prints overlaps before -> after and warns if any remain.
- Type check: `~/data/conda/envs/env_typecheck/bin/mypy`. Commit directly on `master`; push only when asked.
