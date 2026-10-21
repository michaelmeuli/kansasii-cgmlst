"""strip_plasmids drops contigs whose header says plasmid and keeps the rest unchanged."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from common import strip_plasmids  # noqa: E402


def test_strip_plasmids(tmp_path: Path) -> None:
    src = tmp_path / "in.fasta"
    src.write_text(">c1 Mycobacterium kansasii ATCC 12478, complete sequence\nACGT\nAC\n"
                   ">p1 Mycobacterium kansasii ATCC 12478 plasmid pMK12478, complete sequence\nGGGG\n>c2 draft\nTT\n")
    dst = tmp_path / "out" / "out.fasta"
    removed = strip_plasmids(src, dst)
    assert removed == ["p1 Mycobacterium kansasii ATCC 12478 plasmid pMK12478, complete sequence"]
    assert dst.read_text() == ">c1 Mycobacterium kansasii ATCC 12478, complete sequence\nACGT\nAC\n>c2 draft\nTT\n"
