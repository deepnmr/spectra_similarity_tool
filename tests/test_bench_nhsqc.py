"""Guards the expanded two-decoy dense benchmark: LCC must win and keep both decoys below every
same-protein point. Skips when the Nhsqc Bruker data are not present (they are not redistributed)."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import bench_nhsqc  # noqa: E402

DATA = bench_nhsqc.PRL3 / "2" / "pdata" / "1" / "2rr"


@pytest.mark.skipif(not DATA.exists(), reason="Nhsqc dense data not present")
def test_lcc_wins_and_separates_two_decoys():
    decoys = {"OAA": Path(bench_nhsqc.DEF_OAA), "EphB3": Path(bench_nhsqc.DEF_EPHB3)}
    rows = bench_nhsqc.run(bench_nhsqc.PRL3, decoys)["rows"]
    lcc = rows["lcc_new"]
    # LCC has the top separation of all methods
    assert lcc["separation"] == max(r["separation"] for r in rows.values())
    # both decoy proteins fall below every same-protein titration point (positive margin)
    assert lcc["margin"] > 0
    assert lcc["max_diff"] < lcc["min_same"]
    # mean-centring still earns the gap over the un-centred cosine ablation
    assert lcc["separation"] > rows["cosine_uncentred"]["separation"]
