"""Build output/ridom_assemblies.tar.gz and .zip (same content) for Ridom SeqSphere+: the QC-passing cohort assemblies (as real files, named M001; the
ATCC/DSM control strains by their LABEL) + the 7 public reference genomes (named by accession), and metadata.tsv for the
isolate-data import. Names match the node labels of the cgMLST MSTs (09_mst_figures.py)."""
from __future__ import annotations

import argparse
import re
import shutil
import tarfile
import tempfile
import zipfile
from pathlib import Path

import pandas as pd

from common import REF_FASTA, REFERENCES, RIDOM_TARBALL, RIDOM_ZIP, assembly_fasta, load_samples
from kansasii_paths import ASSEMBLIES, SCREENING_MAP_RESULTS

ap = argparse.ArgumentParser()
ap.add_argument("--assemblies", type=Path, default=ASSEMBLIES, help="<sample>/1_unicycler/ and 5_typing/mlst/ per sample")
ap.add_argument("--screening-map-results", type=Path, default=SCREENING_MAP_RESULTS, help="keep set and CheckM contamination")
ap.add_argument("--tarball", type=Path, default=RIDOM_TARBALL)
ap.add_argument("--zip", type=Path, default=RIDOM_ZIP)
args = ap.parse_args()
TARBALL: Path = args.tarball
ZIPFILE: Path = args.zip


def mlst_st(sample: str) -> str:
    """Sequence type of the pipeline's mlst run (scheme mycobacteria_2); '-' when mlst assigns none."""
    f = args.assemblies / sample / "5_typing" / "mlst" / f"mlst_output_{sample}.tsv"
    return str(f.read_text(encoding="utf-8").splitlines()[1].split("\t")[2])


def species_label(sp: str) -> str:
    return f"M. {sp}"


def safe_name(label: str) -> str:
    """'M. kansasii (ATCC 12478)' -> 'M_kansasii_ATCC_12478' (no spaces, dots or brackets in isolate/file names)."""
    return re.sub(r"_+", "_", re.sub(r"[^A-Za-z0-9]+", "_", label)).strip("_")


def copy_fasta(src: Path, dst: Path) -> None:
    shutil.copyfile(src, dst)  # follows symlinks: the tarball holds real files
    assert dst.stat().st_size > 0, src


s = load_samples(args.screening_map_results)
s = s[s.keep]
chk = pd.read_csv(args.screening_map_results, dtype=str).fillna("").set_index("PROBENNUMMER")["checkm_contamination"]
rows: list[dict[str, object]] = []
for sample, species, label in zip(s["sample"].astype(str), s["species"].astype(str), s["label"].astype(str)):
    name = safe_name(label) if label else sample.replace("Mkan329-", "M")
    rows.append({"name": name, "sample": sample, "species": species_label(species), "mlst_ST": mlst_st(sample),
                 "checkm_contamination": chk.at[sample], "is_control": bool(label),
                 "is_reference": False, "src": assembly_fasta(sample, args.assemblies)})
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
    tmp_zip = ZIPFILE.with_name(ZIPFILE.name + ".part")
    with zipfile.ZipFile(tmp_zip, "w", zipfile.ZIP_DEFLATED) as zf:
        for f in sorted((top / "assemblies").iterdir()):
            zf.write(f, arcname=f"ridom/assemblies/{f.name}")
        zf.write(top / "metadata.tsv", arcname="ridom/metadata.tsv")
    tmp_zip.replace(ZIPFILE)

with tarfile.open(TARBALL) as tf:  # verify
    fasta = [m for m in tf.getmembers() if m.name.endswith(".fasta")]
    assert all(m.isfile() and m.size > 0 for m in fasta), "link or empty file in tarball"
    assert len(fasta) == len(df), (len(fasta), len(df))
with zipfile.ZipFile(ZIPFILE) as zf:
    zfasta = [i for i in zf.infolist() if i.filename.endswith(".fasta")]
    assert all(i.file_size > 0 for i in zfasta), "empty file in zip"
    assert sorted(i.filename for i in zfasta) == sorted(m.name for m in fasta), "zip and tarball differ"
print(f"{TARBALL} + {ZIPFILE.name}: {len(df)} genomes ({int(df.is_reference.sum())} references, {int(df.is_control.sum())} controls), "
      f"{TARBALL.stat().st_size / 1e6:.1f} MB; species:\n{df.species.value_counts().to_string()}")
