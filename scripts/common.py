"""Shared paths, QC exclusion and subset definitions (sample QC mirrors kansasii-mst/scripts/common.py)."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from kansasii_paths import ASSEMBLIES, COMPLEX, OUTPUT_ROOT, RUNS, SCREENING_MAP_RESULTS

OUTPUT = OUTPUT_ROOT / "cgmlst"  # final tables and figures only
WORK = RUNS / "cgmlst"  # chewBBACA working state (schemas, genome links, full reports)
GENOMES = WORK / "genomes"  # symlinks <sample>.fasta, input of AlleleCall
SEEDS = WORK / "seed_genomes"  # symlinks, input of CreateSchema
SCHEMA = WORK / "schema"
ALLELE_CALL = WORK / "allele_call"
CGMLST = WORK / "cgmlst"
EVAL = WORK / "schema_evaluator"
REF_FASTA = WORK / "reference_genomes"  # <accession>.fasta: the GTDB representatives without plasmid contigs (step 00)
RIDOM_TARBALL = OUTPUT_ROOT / "ridom_assemblies.tar.gz"
RIDOM_ZIP = OUTPUT_ROOT / "ridom_assemblies.zip"
PTF = WORK / "ATCC12478.trn"  # Prodigal training file

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
    "persicum_attenuatum_ostraviense_innocens_gastri": ["persicum", "attenuatum", "ostraviense", "innocens", "gastri"],
}
# chewBBACA allele-call classifications that are not an allele id; with pairwise deletion they count as missing
MISSING_CODES = ("LNF", "NIPH", "NIPHEM", "ASM", "ALM", "PLOT3", "PLOT5", "LOTSC")


def assembly_fasta(sample: str, assemblies: Path = ASSEMBLIES) -> Path:
    """Unicycler assembly of a cohort sample in the pipeline results."""
    return assemblies / sample / "1_unicycler" / f"{sample}.fasta"


def strip_plasmids(src: Path, dst: Path) -> list[str]:
    """Write src to dst without the contigs whose header says "plasmid"; return the headers of the removed contigs."""
    kept: list[str] = []
    removed: list[str] = []
    for rec in src.read_text().split("\n>"):
        rec = rec.lstrip(">")
        if not rec.strip():
            continue
        header = rec.split("\n", 1)[0]
        (removed if "plasmid" in header.lower() else kept).append(rec.rstrip("\n"))
    assert kept, f"{src}: only plasmid contigs"
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text("".join(f">{r}\n" for r in kept))
    return [r.split("\n", 1)[0] for r in removed]


def load_samples(results: Path = SCREENING_MAP_RESULTS) -> pd.DataFrame:
    """Cohort genomes with columns sample, NR, species, label (LABEL of the ATCC/DSM control strains, else ''), in_pool (complex
    species, not contaminated: the genomes that are allele-called and define the scheme), excluded (bool: `excluded` column of
    screening_map_results.csv = contaminated or < 90 % of the cgMLST loci called), keep (in_pool and not excluded: the genomes
    of the distance matrices, MSTs and Ridom export) and reason (why not kept, else ''). Control strains of complex species are
    kept (unlike kansasii-mst); non-complex ones are not. Without the `excluded` column (first run) nothing is excluded by it."""
    path = results
    cols = ["NR", "PROBENNUMMER", "LABEL", "species", "contamination_flag"]
    extra = [c for c in ("excluded", "cgmlst_called_pct") if c in pd.read_csv(path, nrows=0).columns]
    m = pd.read_csv(path, usecols=cols + extra, dtype={"excluded": str, "cgmlst_called_pct": str})
    pool_reason = pd.Series("", index=m.index)
    pool_reason[~m.species.isin(COMPLEX)] = "not in complex"
    pool_reason[(m.contamination_flag == "contaminated") & (pool_reason == "")] = "contaminated"
    excluded = m["excluded"].fillna("").str.lower().eq("true") if "excluded" in m else pd.Series(False, index=m.index)
    pct = m["cgmlst_called_pct"].fillna("") if "cgmlst_called_pct" in m else pd.Series("", index=m.index)
    reason = pool_reason.copy()
    low = excluded & (reason == "")
    reason[low] = "low cgMLST coverage (" + pct[low] + " % of loci called)"
    return pd.DataFrame({"sample": m.PROBENNUMMER, "NR": m.NR, "species": m.species, "label": m.LABEL.fillna(""),
                         "in_pool": pool_reason == "", "excluded": excluded, "keep": reason == "",
                         "reason": reason}).sort_values("NR").reset_index(drop=True)
