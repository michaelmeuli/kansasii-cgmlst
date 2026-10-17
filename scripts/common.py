"""Shared paths, QC exclusion and subset definitions (sample QC mirrors kansasii-mst/scripts/common.py)."""
from __future__ import annotations

import os
from pathlib import Path

import pandas as pd

ROOT = Path(os.environ.get("KANSASII_ROOT", "/shares/sander.imm.uzh/MM/kansasii"))
OUT = ROOT / "output/cgmlst"
GENOMES = OUT / "genomes"  # symlinks <sample>.fasta, input of AlleleCall
SEEDS = OUT / "seed_genomes"  # symlinks, input of CreateSchema
SCHEMA = OUT / "schema"
ALLELE_CALL = OUT / "allele_call"
CGMLST = OUT / "cgmlst"
EVAL = OUT / "schema_evaluator"
ASSEMBLIES = ROOT / "output/ridom/assemblies"  # <sample>.fasta, the cohort assemblies collected for SeqSphere
GTDB_REPS = ROOT / "data/gtdb_genomes/Mycobacteriaceae/kansasii_complex_gtdb_representatives"
ATCC_FASTA = ROOT / "data/gtdb_genomes/Mycobacteriaceae/ncbi_dataset/data/GCF_000157895.3/GCF_000157895.3_ASM15789v2_genomic.fna"
PTF = OUT / "ATCC12478.trn"  # Prodigal training file

COMPLEX = ["kansasii", "persicum", "pseudokansasii", "attenuatum", "innocens", "ostraviense"]
# public GTDB representatives (file stem, species); the first is the ATCC 12478 type strain of M. kansasii
REFERENCES = {
    "GCF_000157895.3": "kansasii", "GCF_002086675.1": "persicum", "GCF_900566075.1": "pseudokansasii",
    "GCF_900566085.1": "attenuatum", "GCF_900566055.1": "innocens", "GCF_002705925.1": "ostraviense",
    "GCF_002102175.1": "gastri",
}
# one MST per subset: species to keep; None = every genome
SUBSETS: dict[str, list[str] | None] = {
    "complex": None,
    "kansasii": ["kansasii"],
    "persicum": ["persicum"],
    "pseudokansasii": ["pseudokansasii"],
}
# same colours as kansasii-mst / mlsa-kansasii
SPECIES_COLORS = {
    "kansasii": "#1b9e77", "persicum": "#d95f02", "pseudokansasii": "#7570b3", "innocens": "#e7298a",
    "attenuatum": "#66a61e", "ostraviense": "#e6ab02", "gastri": "#a6761d",
}
EXTRA_EXCL = {10, 21}  # low-quality genomes missed by the contamination flag
# chewBBACA allele-call classifications that are not an allele id; with pairwise deletion they count as missing
MISSING_CODES = ("LNF", "NIPH", "NIPHEM", "ASM", "ALM", "PLOT3", "PLOT5", "LOTSC")


def load_samples() -> pd.DataFrame:
    """Cohort genomes with columns sample, NR, species, label (LABEL of the ATCC/DSM control strains, else ''), keep (bool) and
    reason (why excluded, else ''). Control strains of complex species are kept (unlike kansasii-mst); non-complex ones are not."""
    m = pd.read_csv(ROOT / "data/imm/screening_map_results.csv",
                    usecols=["NR", "PROBENNUMMER", "LABEL", "species", "contamination_flag"])
    reason = pd.Series("", index=m.index)
    reason[~m.species.isin(COMPLEX)] = "not in complex"
    reason[(m.contamination_flag == "contaminated") & (reason == "")] = "contaminated"
    reason[m.NR.isin(EXTRA_EXCL) & (reason == "")] = "low quality (manual)"
    return pd.DataFrame({"sample": m.PROBENNUMMER, "NR": m.NR, "species": m.species,
                         "label": m.LABEL.fillna(""), "keep": reason == "", "reason": reason}).sort_values("NR").reset_index(drop=True)
