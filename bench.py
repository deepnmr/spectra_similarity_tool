#!/usr/bin/env python3
"""Benchmark harness: score every method on the PRL3-vs-OAA discrimination task and
report the SEPARATION (mean same-protein score minus different-protein score).

Ground truth: reference PRL3 exp 2 (1H-15N amide HSQC). The same protein plus ligand
(exps 4..14, a titration series) should score HIGH; a different protein (OAA exp 103)
should score LOW. A better method maximises the separation while self-scoring exactly 1.

Data lives outside the repo; paths are the defaults below, override with --prl3 / --oaa.
Run with python3.11 (the interpreter that has numpy/matplotlib here).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from hsqc_similarity import hsqc_similarity, read_bruker_2d
from hsqc_methods import nn_peak_similarity, quadtree_similarity

WINDOW = dict(range_f2=(6.5, 10.0), range_f1=(105.0, 130.0))
SAME = [4, 6, 8, 10, 12, 14]  # PRL3 + ligand titration
REF = 2


def _methods():
    """name -> callable(x, y) -> similarity float. New methods get appended here."""
    m = {
        "bin_Bodis09": lambda x, y: hsqc_similarity(x, y, **WINDOW)["similarity"],
        "bin_rot45": lambda x, y: hsqc_similarity(x, y, rotate_deg=45.0, **WINDOW)["similarity"],
        "tree_Castillo13": lambda x, y: quadtree_similarity(x, y, **WINDOW)["similarity"],
        "nn_Pierens12": lambda x, y: nn_peak_similarity(x, y, **WINDOW)["similarity"],
    }
    try:
        from hsqc_lcc import lcc_similarity  # the new synthesized method
        m["lcc_new"] = lambda x, y: lcc_similarity(x, y, **WINDOW)["similarity"]
    except Exception:
        pass
    return m


def run(prl3: Path, oaa: Path) -> dict:
    ref = read_bruker_2d(prl3 / str(REF))
    same = {e: read_bruker_2d(prl3 / str(e)) for e in SAME}
    diff = read_bruker_2d(oaa / "103")

    methods = _methods()
    rows: dict[str, dict[str, float]] = {}
    sep: dict[str, float] = {}
    for name, fn in methods.items():
        self_s = fn(ref, ref)
        same_scores = {f"2v{e}": fn(ref, same[e]) for e in SAME}
        diff_s = fn(ref, diff)
        mean_same = sum(same_scores.values()) / len(same_scores)
        rows[name] = {"self": self_s, **same_scores, "2vOAA": diff_s,
                      "mean_same": mean_same, "separation": mean_same - diff_s}
        sep[name] = mean_same - diff_s
    return {"rows": rows, "separation": sep}


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--prl3", type=Path, default=Path("/Users/donghanlee/PRL3_mark_5nov24"))
    p.add_argument("--oaa", type=Path, default=Path("/Users/donghanlee/Downloads/OAA_CEST_277K_04may15"))
    p.add_argument("--json", action="store_true")
    args = p.parse_args()

    result = run(args.prl3, args.oaa)
    if args.json:
        print(json.dumps(result, indent=2))
        return 0

    rows = result["rows"]
    cols = ["self", "2v4", "2v6", "2v8", "2v10", "2v12", "2v14", "2vOAA", "separation"]
    print(f"{'method':18} " + " ".join(f"{c:>8}" for c in cols))
    # rank by separation, best first
    for name in sorted(rows, key=lambda n: -rows[n]["separation"]):
        r = rows[name]
        print(f"{name:18} " + " ".join(f"{r[c]:8.4f}" for c in cols))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
