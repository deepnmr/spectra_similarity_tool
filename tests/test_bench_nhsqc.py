"""Guards the expanded two-decoy dense benchmark: LCC must win and keep both decoys below every
same-protein point. Skips when the Nhsqc Bruker data are not present (they are not redistributed)."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import bench_nhsqc  # noqa: E402

DATA = bench_nhsqc.PRL3 / "2" / "pdata" / "1" / "2rr"


def test_dense_registry_exposes_local_contrast():
    assert "local_contrast" in bench_nhsqc._methods()


@pytest.mark.skipif(not DATA.exists(), reason="Nhsqc dense data not present")
def test_local_contrast_wins_and_separates_two_decoys():
    decoys = {"OAA": Path(bench_nhsqc.DEF_OAA), "EphB3": Path(bench_nhsqc.DEF_EPHB3)}
    rows = bench_nhsqc.run(bench_nhsqc.PRL3, decoys)["rows"]
    candidate = rows["local_contrast"]
    assert candidate["separation"] == max(r["separation"] for r in rows.values())
    # both decoy proteins fall below every same-protein titration point (positive margin)
    assert candidate["margin"] > 0
    assert candidate["max_diff"] < candidate["min_same"]
    assert candidate["separation"] > rows["lcc_new"]["separation"]
