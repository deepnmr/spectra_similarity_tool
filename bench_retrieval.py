#!/usr/bin/env python3
"""Dereplication / library-retrieval experiment on the sparse 1H-13C set + bootstrap CIs.

Answers the reviewer question the separation score does not: as a LIBRARY SEARCH, does each
measure retrieve the correct same-compound partner? For every one of the 10 spectra we rank the
other 9 by similarity and ask whether the true same-compound partner comes out on top.

Metrics per method:
  - top1  : fraction of queries whose rank-1 hit is the same-compound partner (retrieval accuracy)
  - mrr   : mean reciprocal rank of the true partner (1.0 = always rank 1)
  - separation with a bootstrap 95% CI (resampling the 45 unordered pairs)

Reuses bench_13c's loader/methods/pairs (olivetol duplicate already excluded there). Slow because
tree/NN are slow on the large grid; run once, writes results/retrieval_13c.json.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np

from bench_13c import METHODS, PAIRS, load_peaklist, run as _agg_run  # noqa: F401


def _matrix(specs, names, fn):
    n = len(names)
    S = np.eye(n)
    for i, j in itertools.combinations(range(n), 2):
        v = fn(specs[names[i]], specs[names[j]])
        S[i, j] = S[j, i] = v
    return S


def retrieval(S, names, comp):
    """top-1 accuracy and MRR of the same-compound partner over all queries."""
    n = len(names)
    top1 = 0
    rr = []
    for i in range(n):
        order = sorted((k for k in range(n) if k != i), key=lambda k: -S[i, k])
        # the (single) true partner: the other spectrum of the same compound
        partner = [k for k in range(n) if k != i and comp[names[k]] == comp[names[i]]]
        assert len(partner) == 1
        rank = order.index(partner[0]) + 1
        rr.append(1.0 / rank)
        top1 += int(rank == 1)
    return {"top1": top1 / n, "mrr": float(np.mean(rr))}


def bootstrap_sep(same_vals, diff_vals, iters=10000, seed=0):
    rng = np.random.default_rng(seed)
    s, d = np.asarray(same_vals), np.asarray(diff_vals)
    seps = np.empty(iters)
    for b in range(iters):
        seps[b] = s[rng.integers(0, s.size, s.size)].mean() - d[rng.integers(0, d.size, d.size)].mean()
    lo, hi = np.percentile(seps, [2.5, 97.5])
    return {"sep": float(s.mean() - d.mean()), "ci95": [float(lo), float(hi)]}


def main() -> int:
    data_dir = Path("data_13c")
    names = [v for _, a, b in PAIRS for v in (a[0], b[0])]
    comp = {v: c for c, a, b in PAIRS for v in (a[0], b[0])}
    specs = {n: load_peaklist(data_dir / f"{n}.json") for n in names}

    same_idx = [(names.index(a[0]), names.index(b[0])) for _, a, b in PAIRS]
    same_set = {frozenset(p) for p in same_idx}

    out = {}
    for m, fn in METHODS.items():
        S = _matrix(specs, names, fn)
        same_vals = [S[i, j] for i, j in same_idx]
        diff_vals = [S[i, j] for i, j in itertools.combinations(range(len(names)), 2)
                     if frozenset((i, j)) not in same_set]
        r = retrieval(S, names, comp)
        b = bootstrap_sep(same_vals, diff_vals)
        out[m] = {**r, **b}
        print(f"{m:18} top1={r['top1']:.2f} mrr={r['mrr']:.3f} "
              f"sep={b['sep']:.3f} CI95=[{b['ci95'][0]:.3f},{b['ci95'][1]:.3f}]")

    Path("results").mkdir(exist_ok=True)
    json.dump({"rows": out, "n_query": len(names), "n_same": len(same_idx)},
              open("results/retrieval_13c.json", "w"), indent=2)
    print("\nwrote results/retrieval_13c.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
