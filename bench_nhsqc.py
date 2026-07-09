#!/usr/bin/env python3
"""Expanded dense 1H-15N benchmark with MULTIPLE decoy proteins (addresses the n=1 negative).

A PRL3 apo->olsalazine titration (23 same-protein points, should score HIGH) scored against TWO
distinct decoy folds -- OAA and EphB3 (should score LOW). Unlike the single-negative benchmark the
negative side is now a distribution over two folds, so separation/margin gain a real negative spread.

Data paths default to the ./Nhsqc tree but, like bench.py, are overridable via flags / environment
(the raw Bruker data are not redistributed); the harness errors cleanly when they are absent.
Reuses bench.py's method set and window. Run with python3.11.
"""
from __future__ import annotations

import argparse
import json
import os
import statistics
from pathlib import Path

from bench import _methods, WINDOW  # same 6 methods, same 6.5-10 x 105-130 ppm window
from hsqc_similarity import read_bruker_2d

HERE = Path(__file__).parent
NH = HERE / "Nhsqc"
REF = 2
SAME = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48]

# defaults point at the in-repo Nhsqc tree; override with flags or $NHSQC_* env vars
DEF_TITRATION = os.environ.get("NHSQC_TITRATION", str(NH / "PRL3_olsalazine_titration_hsqcs"))
DEF_OAA = os.environ.get("NHSQC_OAA", str(NH / "OAA_CEST_277K_04may15" / "103"))
DEF_EPHB3 = os.environ.get("NHSQC_EPHB3", str(NH / "jhchoi_230105_ephb3_check_900" / "11"))

# module-level so the test suite can probe availability without parsing args
PRL3 = Path(DEF_TITRATION)


def run(titration: Path, decoys: dict[str, Path]) -> dict:
    ref = read_bruker_2d(titration / str(REF))
    same = {e: read_bruker_2d(titration / str(e)) for e in SAME}
    decoy = {name: read_bruker_2d(p) for name, p in decoys.items()}

    rows: dict[str, dict] = {}
    for name, fn in _methods().items():
        self_s = fn(ref, ref)
        assert abs(self_s - 1.0) < 1e-6, f"{name} self-score {self_s} != 1"
        same_scores = {e: fn(ref, same[e]) for e in SAME}
        decoy_scores = {d: fn(ref, decoy[d]) for d in decoys}
        ss = list(same_scores.values())
        ds = list(decoy_scores.values())
        mean_same, min_same = statistics.mean(ss), min(ss)
        mean_diff, max_diff = statistics.mean(ds), max(ds)  # worst (highest) decoy = hardest negative
        rows[name] = {
            "self": self_s,
            "same": {f"2v{e}": round(v, 4) for e, v in same_scores.items()},
            "decoy": {d: round(v, 4) for d, v in decoy_scores.items()},
            "mean_same": round(mean_same, 4), "min_same": round(min_same, 4),
            "mean_diff": round(mean_diff, 4), "max_diff": round(max_diff, 4),
            "separation": round(mean_same - mean_diff, 4),
            "margin": round(min_same - max_diff, 4),
        }
    return {"reference": f"PRL3 apo exp {REF}", "n_same": len(SAME),
            "decoys": list(decoys), "window": WINDOW, "rows": rows}


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--titration", type=Path, default=DEF_TITRATION,
                   help="PRL3 apo+olsalazine titration dir (or $NHSQC_TITRATION)")
    p.add_argument("--oaa", type=Path, default=DEF_OAA, help="OAA decoy experiment dir (or $NHSQC_OAA)")
    p.add_argument("--ephb3", type=Path, default=DEF_EPHB3, help="EphB3 decoy experiment dir (or $NHSQC_EPHB3)")
    p.add_argument("--json", action="store_true")
    args = p.parse_args()

    # clean failure on a public checkout that lacks the (non-redistributed) Bruker data
    missing = [str(d) for d in (args.titration / str(REF), args.oaa, args.ephb3) if not Path(d).exists()]
    if missing:
        p.error("expanded dense benchmark needs the Nhsqc Bruker data, not found at: "
                + "; ".join(missing) + ". Pass --titration/--oaa/--ephb3 (or set $NHSQC_TITRATION/"
                "$NHSQC_OAA/$NHSQC_EPHB3). The data are available from the author on request.")

    result = run(args.titration, {"OAA": args.oaa, "EphB3": args.ephb3})
    (HERE / "results" / "nhsqc_dense.json").write_text(json.dumps(result, indent=2))
    if args.json:
        print(json.dumps(result, indent=2))
        return 0
    rows = result["rows"]
    cols = ["mean_same", "min_same", "OAA", "EphB3", "mean_diff", "max_diff", "separation", "margin"]
    print(f"{'method':18} " + " ".join(f"{c:>9}" for c in cols))
    for name in sorted(rows, key=lambda n: -rows[n]["separation"]):
        r = rows[name]
        vals = [r["mean_same"], r["min_same"], r["decoy"]["OAA"], r["decoy"]["EphB3"],
                r["mean_diff"], r["max_diff"], r["separation"], r["margin"]]
        print(f"{name:18} " + " ".join(f"{v:9.4f}" for v in vals))
    print("\nwrote results/nhsqc_dense.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
