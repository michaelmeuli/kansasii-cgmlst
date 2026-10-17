# kansasii-cgmlst

cgMLST scheme for the *Mycobacterium kansasii* complex built with chewBBACA (3.5.4), and Ridom SeqSphere-style
minimum spanning trees (MSTs) of the allele profiles. Complements `kansasii-mst` (core-SNP MSTs).

Outputs go to `$KANSASII_ROOT/output/cgmlst/` (default `/shares/sander.imm.uzh/MM/kansasii`), figures in `output/cgmlst/mst/`.
Run: `sbatch run_all.sbatch` (env `env_cgmlst`, see `environment.yml`); `STEPS="07 08 09" sbatch --export=ALL,STEPS="07 08 09" run_all.sbatch` runs a subset.

## Method

| Step | Script | What |
|---|---|---|
| 00 | `00_prepare_inputs.py` | 172 QC-passing cohort assemblies (QC as in `kansasii-mst`: no contaminated, controls, non-complex, NR 10/21) + 7 public GTDB representatives (ATCC 12478 = *kansasii*, plus persicum, pseudokansasii, attenuatum, innocens, ostraviense, gastri) = 179 genomes |
| 01 | `01_prodigal_training.sh` | `prodigal -p single` on ATCC 12478 (GCF_000157895.3) -> `ATCC12478.trn` |
| 02 | `02_create_schema.sh` | `CreateSchema` on the 7 public representatives (seed genomes) with `--ptf ATCC12478.trn` -> 9,216 loci |
| 03 | `03_allele_call.sh` | `AlleleCall` of all 179 genomes (`--ptf`, `--output-masked`); runs on a copy of the schema because new alleles are added |
| 04 | `04_remove_paralogs.sh` | `RemoveGenes` for the 69 loci in `paralogous_counts.tsv` |
| 05 | `05_extract_cgmlst.sh` | `ExtractCgMLST --t 1.0 / 0.99 / 0.95`; 0.95 is the scheme used downstream (`CGMLST_THRESHOLD` in `env.sh`) |
| 06 | `06_schema_evaluator.sh` | `SchemaEvaluator --ta 11` on the final cgMLST loci |
| 07 | `07_summarise_scheme.py` | good-target percentage -> `scheme_summary.tsv`, `schema_evaluator_loci.tsv` |
| 08 | `08_distance_matrix.py` | pairwise-deletion allele distances (`allele_distance.py`) per subset |
| 09 | `09_mst_figures.py` | MST figures |

## Scheme

| metric | value |
|---|---|
| loci in schema (seed) | 9,216 |
| paralogous loci removed | 69 (9,147 remain) |
| cgMLST loci, present in 100 % / 99 % / 95 % of genomes | 1,498 / 2,578 / **3,610** |
| **good targets** (SchemaEvaluator: no flagged allele) | **3,558 / 3,610 = 98.6 %** (all alleles valid; 6 loci with an allele < 201 bp, 47 with an allele above the size threshold) |

98.6 % is above the 95-98 % usually quoted for published schemes. Schema and cgMLST were defined on the same 179 genomes, so the
95 % presence rule is satisfied by construction; validate on new isolates before treating the scheme as fixed.

## Distances: pairwise deletion

As in Ridom, a locus counts for a pair of genomes only if both have an allele call. Everything that is not an allele id
(LNF, NIPH, NIPHEM, ASM, ALM, PLOT3, PLOT5, LOTSC, 0 in the masked table) is missing. Distance = number of called-in-both loci with
different alleles (no normalisation). `shared_loci_<subset>.csv` gives the number of loci compared per pair
(minimum 2,668 of 3,610 in the complex). `tests/test_allele_distance.py` checks this against hand-computed values.

## Figures (`output/cgmlst/mst/`)

- `mst_cgmlst_complex.{png,pdf}`: whole complex incl. the public references (dashed ring).
- `mst_cgmlst_kansasii.{png,pdf}`, `mst_cgmlst_persicum.{png,pdf}`, `mst_cgmlst_pseudokansasii.{png,pdf}`: one MST per species
  (each with its own public reference genome, if any). Distances are computed within the subset, so shared loci differ per figure.
- `mst_nodes_<subset>.tsv`: node membership.

Layout: Kamada-Kawai on log edge lengths, then `layout_fix.py` nudges any node off non-incident edges and off other nodes (0 overlaps left in all figures).

Ridom-like features: node area ~ genome count, identical profiles share a node, pies for mixed-species nodes, edge labels = allele
differences, logarithmic edge length (`--linear` to switch), MST ties broken towards larger nodes, optional cluster shading
(`--cluster-threshold N`). Not identical to SeqSphere: its scheme and layout algorithm are proprietary, so distances and layout will differ.

## Caveats

- Nearly every isolate has a unique profile, so nodes are single genomes and MST edge labels are large
  (tens to hundreds of alleles); this reflects strong recombination in the complex, and allele distances are not recombination-filtered.
- ATCC 12478 carries unusual alleles; it is only the Prodigal training source and a seed genome.
- Hub regions have crowded edge labels; use the PDF to zoom. Node labels sit in a box in their species colour and are joined to their node by a dotted leader line (labels are pushed apart automatically, so they do not always touch the node).
- The cgMLST is defined at the complex level, so between-species distances (thousands of alleles) rest on loci conserved across species.
