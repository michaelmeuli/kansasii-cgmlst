"""Link the QC-passing cohort assemblies and the public reference genomes into runs/cgmlst/genomes/, the species
representatives (plasmid contigs removed, runs/cgmlst/reference_genomes/) into runs/cgmlst/seed_genomes/ (CreateSchema input), and write samples.tsv (sample, NR, species, source, is_reference)."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from common import ASSEMBLIES, GENOMES, GTDB_REPS, OUT, REF_FASTA, REFERENCES, SEEDS, load_samples, strip_plasmids


def link(src: Path, dst: Path) -> None:
    if dst.is_symlink() or dst.exists():
        dst.unlink()
    dst.symlink_to(src)


GENOMES.mkdir(parents=True, exist_ok=True)
SEEDS.mkdir(parents=True, exist_ok=True)
for old in GENOMES.glob("*.fasta"):  # drop links of genomes that no longer pass QC (AlleleCall reads the whole directory)
    if old.is_symlink():
        old.unlink()
rows: list[dict[str, object]] = []
s = load_samples()
for r in s[s.keep].itertuples():
    src = ASSEMBLIES / f"{r.sample}.fasta"
    assert src.exists(), src
    link(src, GENOMES / f"{r.sample}.fasta")
    rows.append({"sample": r.sample, "NR": r.NR, "species": r.species, "label": r.label,
                 "is_control": bool(r.label), "source": str(src), "is_reference": False})
for acc, sp in REFERENCES.items():
    raw = GTDB_REPS / f"{acc}_query.fasta"
    assert raw.exists(), raw
    src = REF_FASTA / f"{acc}.fasta"
    for h in strip_plasmids(raw, src):  # plasmid genes are not core genome: keep them out of seeds, training and allele calls
        print(f"{acc}: plasmid contig removed: {h}")
    link(src, GENOMES / f"{acc}.fasta")
    link(src, SEEDS / f"{acc}.fasta")
    rows.append({"sample": acc, "NR": None, "species": sp, "label": "", "is_control": False,
                 "source": str(src), "is_reference": True})
df = pd.DataFrame(rows)
OUT.mkdir(parents=True, exist_ok=True)
df.to_csv(OUT / "samples.tsv", sep="\t", index=False)
print(f"{len(df)} genomes ({int(df.is_reference.sum())} references); species:\n{df.species.value_counts().to_string()}")
