# kansasii-cgmlst

cgMLST scheme for the *Mycobacterium kansasii* complex built with chewBBACA (3.5.4), and Ridom SeqSphere-style
minimum spanning trees (MSTs) of the allele profiles. Complements `kansasii-mst` (core-SNP MSTs).

Final tables and figures go to `output/cgmlst/` (`OUTPUT`; `KANSASII_ROOT` moves the whole tree), figures in `output/cgmlst/mst/`.
The chewBBACA working state (genome links, schemas, `schema_called`, per-locus reports, `cds_coordinates.tsv`, logs; ~350 MB) lives in `runs/cgmlst/`;
the scripts copy the result tables from there to `output/cgmlst/` (`publish` in `scripts/env.sh`).
Run: `sbatch submit_all.sbatch` (env `env_cgmlst`, see `environment.yml`); `sbatch --export=ALL,STEPS=08:09:10 submit_all.sbatch` runs a subset (steps separated by `:`).
Cohort assemblies are read from the pipeline results (`ASSEMBLIES/<sample>/1_unicycler/`, flag `--assemblies`), the 7 public references from
`GTDB_REPRESENTATIVES` (flag `--gtdb-representatives`), the sample table from `SCREENING_MAP_RESULTS`; the former hand-made `output/ridom/` symlink
directory is no longer read (the mlst sequence type comes from `ASSEMBLIES/<sample>/5_typing/mlst/`).
## Layout

```
scripts/          numbered steps (00-10), env.sh, common.py, allele_distance.py, label_fix.py, layout_fix.py; kansasii_paths.py = shared constants (vendored)
submit_all.sbatch the whole pipeline (log in runs/cgmlst/logs/)
tests/            pytest
```

Paths are constants of `scripts/kansasii_paths.py` (same names in every repo, `ASSEMBLIES`, `OUTPUT_ROOT`, ...) re-used in
`scripts/common.py`, and flags of the scripts that read them (`--assemblies`, `--gtdb-representatives`, `--screening-map-results`, `--tarball`). Edit `kansasii-docs/templates/kansasii_paths.py`,
not the vendored copy (`bash repos/check_sync.sh --fix`).


## Method

| Step | Script | What |
|---|---|---|
| 00 | `00_prepare_inputs.py` | 175 QC-passing cohort assemblies (QC as in `kansasii-mst`: no contaminated (CheckM > 5 %), non-complex; unlike `kansasii-mst` the ATCC/DSM control strains of complex species, NR 128/132/133, are kept and shown by their `LABEL`) + 7 public GTDB representatives, plasmid contigs removed (only ATCC 12478 pMK12478, NC_022654.1, is labelled as a plasmid; `strip_plasmids` in `common.py` matches "plasmid" in the FASTA header) (ATCC 12478 = *kansasii*, plus persicum, pseudokansasii, attenuatum, innocens, ostraviense, gastri) = 182 genomes |
| 01 | `01_prodigal_training.sh` | `prodigal -p single` on ATCC 12478 (GCF_000157895.3) -> `ATCC12478.trn` |
| 02 | `02_create_schema.sh` | `CreateSchema` on the 7 public representatives (seed genomes) with `--ptf ATCC12478.trn` -> 9,072 loci |
| 03 | `03_allele_call.sh` | `AlleleCall` of all 182 genomes (`--ptf`, `--output-masked`); runs on a copy of the schema because new alleles are added |
| 04 | `04_remove_paralogs.sh` | `RemoveGenes` for the 69 loci in `paralogous_counts.tsv` |
| 05 | `05_extract_cgmlst.sh` | `ExtractCgMLST --t 1.0 / 0.99 / 0.95`; 0.95 is the scheme used downstream (`CGMLST_THRESHOLD` in `env.sh`) |
| 06 | `06_schema_evaluator.sh` | `SchemaEvaluator --ta 11` on the final cgMLST loci |
| 07 | `07_summarise_scheme.py` | good-target percentage -> `scheme_summary.tsv`, `schema_evaluator_loci.tsv` |
| 08 | `08_distance_matrix.py` | pairwise-deletion allele distances (`allele_distance.py`) per subset |
| 09 | `09_mst_figures.py` | MST figures |
| 10 | `10_export_ridom.py` | `output/ridom_assemblies.tar.gz` (manifest `output/ridom_assemblies.MANIFEST.txt`) for Ridom SeqSphere+: the same QC-passing cohort + 7 references as real FASTA files (`ridom/assemblies/`); isolates named `M001` (as in the MST labels), control strains by `LABEL` (`M_kansasii_ATCC_12478`, ...), references by accession; `ridom/metadata.tsv` (name, sample, species as `M. kansasii` ..., ST, CheckM, flags) for the isolate-data import. Written to `.part` and renamed when complete (~20 min, run via sbatch) |

## Sample exclusion (`excluded`)

Genomes that are allele-called (the *pool*: complex species, `contamination_flag` empty) define the scheme. Step 07 writes
`output/cgmlst/sample_qc.tsv` (`cgmlst_called_pct`: % of the cgMLST95 loci with an allele call); immensekansasii
`06_screening_map_results.py` copies it into `data/imm/screening_map_results.csv` and sets `excluded` = contaminated,
`cgmlst_called_pct` < 90 (Ridom's good-target threshold) or a non-complex clinical sample (NR 184). Steps 08-10 (distances, MSTs, Ridom export) skip `excluded`
genomes (`load_samples().keep`); the scheme itself still counts them in the 95 % presence rule (negligible: >= 173 of 182
genomes). Order: `STEPS=00:...:07` -> 06 -> `STEPS=08:09:10` (`run_downstream.sh` does this). M055 (81.5 %) is the only isolate excluded by coverage.

## Scheme

| metric | value |
|---|---|
| loci in schema (seed) | 9,072 |
| paralogous loci removed | 69 (9,003 remain) |
| cgMLST loci, present in 100 % / 99 % / 95 % of genomes | 1,494 / 2,569 / **3,647** |
| **good targets** (SchemaEvaluator: no flagged allele) | **3,596 / 3,647 = 98.6 %** (all alleles valid; 6 loci with an allele < 201 bp, 46 with an allele above the size threshold) |

98.6 % is above the 95-98 % usually quoted for published schemes. Schema and cgMLST were defined on the same 182 genomes, so the
95 % presence rule is satisfied by construction; validate on new isolates before treating the scheme as fixed.

## Distances: pairwise deletion

As in Ridom, a locus counts for a pair of genomes only if both have an allele call. Everything that is not an allele id
(LNF, NIPH, NIPHEM, ASM, ALM, PLOT3, PLOT5, LOTSC, 0 in the masked table) is missing. Distance = number of called-in-both loci with
different alleles (no normalisation). `shared_loci_<subset>.csv` gives the number of loci compared per pair
(minimum 3,136 of 3,647 in the complex). `tests/test_allele_distance.py` checks this against hand-computed values.

## Figures (`output/cgmlst/mst/`)

- `mst_cgmlst_complex.{png,pdf}`: whole complex incl. the public references (dashed ring).
- `mst_cgmlst_kansasii.{png,pdf}`, `mst_cgmlst_persicum.{png,pdf}`, `mst_cgmlst_pseudokansasii.{png,pdf}`: one MST per species
  (each with its own public reference genome, if any). Distances are computed within the subset, so shared loci differ per figure.
- `mst_cgmlst_persicum_attenuatum_ostraviense_innocens_gastri.{png,pdf}`: all non-kansasii, non-pseudokansasii species together (gastri = public reference only).
- Node labels are sample IDs as `M001` (Mkan329-001).
- `mst_nodes_<subset>.tsv`: node membership.

Layout: Kamada-Kawai on log edge lengths, then `layout_fix.py` nudges any node off non-incident edges and off other nodes (0 overlaps left in all figures).

Ridom-like features: node area ~ genome count, identical profiles share a node, pies for mixed-species nodes, edge labels = allele
differences, logarithmic edge length (`--linear` to switch), MST ties broken towards larger nodes, optional cluster shading
(`--cluster-threshold N`). Not identical to SeqSphere: its scheme and layout algorithm are proprietary, so distances and layout will differ.

## Control strains

The ATCC/DSM control strains NR 128 (*M. kansasii* ATCC 12478), 132 (*M. persicum* DSM 104278) and 133 (*M. pseudokansasii* DSM 107152) are
drawn with their `LABEL` text in a bold, thick-framed box and a star ring around the node (legend: control strain). NR 126, 127, 129-131
(Mtb, avium, abscessus, xenopi, ulcerans) are not complex members and are not included. Sanity check (`mst/control_check.tsv`), allele
differences to the public genome of the same species: 128 vs GCF_000157895.3 (ATCC 12478): 3; 133 vs GCF_900566075.1: 8; 132 vs
GCF_002086675.1: 117 (that public persicum genome is strain AFPC-000227, not DSM 104278, so a larger distance is expected).

## Caveats

- Only pMK12478 (labelled plasmid in the FASTA) is removed from the references. The 117 kb contig `NZ_UPHU01000002.1` of the pseudokansasii reference MK142 (conjugation genes, ESX cluster, no rRNA/tRNA; in 8 of 20 pseudokansasii and partly in 6 other isolates) is not labelled by NCBI and is kept on purpose.
- Nearly every isolate has a unique profile, so nodes are single genomes and MST edge labels are large
  (tens to hundreds of alleles); this reflects strong recombination in the complex, and allele distances are not recombination-filtered.
- ATCC 12478 carries unusual alleles; it is only the Prodigal training source and a seed genome.
- Hub regions have crowded edge labels; use the PDF to zoom. Node labels sit in a box in their species colour and are joined to their node by a dotted leader line (labels are pushed apart automatically, so they do not always touch the node).
- The cgMLST is defined at the complex level, so between-species distances (thousands of alleles) rest on loci conserved across species.

## Provenance

The sbatch files write `MANIFEST.txt` (commit, clean/dirty flag, script, conda env hash, inputs) into the output
directory when they finish; see `kansasii-docs/docs/PROVENANCE.md`. Project index: `kansasii-docs/README.md`.
