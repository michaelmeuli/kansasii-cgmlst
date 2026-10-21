"""Build output/ridom_assemblies.tar.gz for Ridom SeqSphere+: the QC-passing cohort assemblies (as real files, named M001; the
ATCC/DSM control strains by their LABEL) + the 7 public reference genomes (named by accession), and metadata.tsv for the
isolate-data import. Names match the node labels of the cgMLST MSTs (09_mst_figures.py)."""
from __future__ import annotations

import re
import shutil
import tarfile
import tempfile
from pathlib import Path

import pandas as pd

from common import ASSEMBLIES, REF_FASTA, REFERENCES, ROOT, load_samples

TARBALL = ROOT / "output/ridom_assemblies.tar.gz"
MANIFEST = ROOT / "output/ridom/manifest.tsv"


def species_label(sp: str) -> str:
    return f"M. {sp}"


def safe_name(label: str) -> str:
    """'M. kansasii (ATCC 12478)' -> 'M_kansasii_ATCC_12478' (no spaces, dots or brackets in isolate/file names)."""
    return re.sub(r"_+", "_", re.sub(r"[^A-Za-z0-9]+", "_", label)).strip("_")


def copy_fasta(src: Path, dst: Path) -> None:
    shutil.copyfile(src, dst)  # follows symlinks: the tarball holds real files
    assert dst.stat().st_size > 0, src


s = load_samples()
s = s[s.keep]
man = pd.read_csv(MANIFEST, sep="\t", dtype=str).set_index("sample")
rows: list[dict[str, object]] = []
for sample, species, label in zip(s["sample"].astype(str), s["species"].astype(str), s["label"].astype(str)):
    name = safe_name(label) if label else sample.replace("Mkan329-", "M")
    rows.append({"name": name, "sample": sample, "species": species_label(species), "mlst_ST": man.at[sample, "mlst_ST"],
                 "checkm_contamination": man.at[sample, "checkm_contamination"], "is_control": bool(label),
                 "is_reference": False, "src": ASSEMBLIES / f"{sample}.fasta"})
for acc, sp in REFERENCES.items():
    rows.append({"name": acc, "sample": acc, "species": species_label(sp), "mlst_ST": "", "checkm_contamination": "",
                 "is_control": False, "is_reference": True, "src": REF_FASTA / f"{acc}.fasta"})
df = pd.DataFrame(rows)
assert df.name.is_unique, df.name[df.name.duplicated()].tolist()

with tempfile.TemporaryDirectory(prefix="ridom_export_") as tmp:
    top = Path(tmp) / "ridom"
    (top / "assemblies").mkdir(parents=True)
    for name, src in zip(df["name"].astype(str), df["src"]):
        copy_fasta(Path(str(src)), top / "assemblies" / f"{name}.fasta")
    df.drop(columns="src").to_csv(top / "metadata.tsv", sep="\t", index=False)
    TARBALL.parent.mkdir(parents=True, exist_ok=True)
    tmp_tar = TARBALL.with_name(TARBALL.name + ".part")  # renamed only when complete
    with tarfile.open(tmp_tar, "w:gz") as tf:
        tf.add(top, arcname="ridom", recursive=False)
        tf.add(top / "assemblies", arcname="ridom/assemblies", recursive=False)
        for f in sorted((top / "assemblies").iterdir()):
            tf.add(f, arcname=f"ridom/assemblies/{f.name}")
        tf.add(top / "metadata.tsv", arcname="ridom/metadata.tsv")
    tmp_tar.replace(TARBALL)

with tarfile.open(TARBALL) as tf:  # verify
    fasta = [m for m in tf.getmembers() if m.name.endswith(".fasta")]
    assert all(m.isfile() and m.size > 0 for m in fasta), "link or empty file in tarball"
    assert len(fasta) == len(df), (len(fasta), len(df))
print(f"{TARBALL}: {len(df)} genomes ({int(df.is_reference.sum())} references, {int(df.is_control.sum())} controls), "
      f"{TARBALL.stat().st_size / 1e6:.1f} MB; species:\n{df.species.value_counts().to_string()}")
