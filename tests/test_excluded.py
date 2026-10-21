"""The `excluded` column of data/imm/screening_map_results.csv: contaminated or < 90 % of the cgMLST loci called.
Compared with the previous rule (contamination only), M055 must be the only additionally excluded complex isolate."""
from __future__ import annotations

import os
from pathlib import Path

import pandas as pd
import pytest

CSV = Path(os.environ.get("KANSASII_ROOT", "/shares/sander.imm.uzh/MM/kansasii")) / "data/imm/screening_map_results.csv"
COMPLEX = ["kansasii", "persicum", "pseudokansasii", "attenuatum", "innocens", "ostraviense"]

pytestmark = pytest.mark.skipif(not CSV.exists(), reason="screening_map_results.csv not available")


def load() -> pd.DataFrame:
    m = pd.read_csv(CSV, dtype=str).fillna("")
    assert "excluded" in m.columns and "cgmlst_called_pct" in m.columns, "rerun immensekansasii 06_screening_map_results.py"
    assert set(m.excluded) <= {"true", "false"}
    return m[m.species.isin(COMPLEX)]


def test_only_m055_additionally_excluded() -> None:
    m = load()
    excluded = set(m.PROBENNUMMER[m.excluded == "true"])
    contaminated = set(m.PROBENNUMMER[m.contamination_flag == "contaminated"])
    assert contaminated <= excluded
    assert excluded - contaminated == {"Mkan329-055"}


def test_threshold() -> None:
    m = load()
    pct = pd.to_numeric(m.cgmlst_called_pct, errors="coerce")
    assert (m.excluded[pct < 90] == "true").all()
    assert (m.excluded[(pct >= 90) & (m.contamination_flag == "")] == "false").all()
