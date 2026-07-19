#!/usr/bin/env python3
"""Dereplication / library-retrieval experiment on the sparse 1H-13C set + bootstrap CIs.

Answers the reviewer question the separation score does not: as a LIBRARY SEARCH, does each
measure retrieve the correct same-compound partner? For every one of the 10 spectra we rank the
other 9 by similarity and ask whether the true same-compound partner comes out on top.

Metrics per method:
  - top1  : fraction of queries whose rank-1 hit is the same-compound partner (retrieval accuracy)
  - mrr   : mean reciprocal rank of the true partner (1.0 = always rank 1)
  - AUROC / AUPRC and leave-one-compound-out threshold classification
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


def ranking_metrics(same_vals, diff_vals):
    """Pair-level AUROC and tie-grouped average precision (AUPRC)."""
    s, d = map(np.asarray, (same_vals, diff_vals))
    auroc = ((s[:, None] > d).mean() + 0.5 * (s[:, None] == d).mean())
    scores = np.concatenate((s, d))
    positive = np.arange(scores.size) < s.size
    auprc = previous_recall = 0.0
    for threshold in np.unique(scores)[::-1]:
        selected = scores >= threshold
        tp = np.count_nonzero(selected & positive)
        recall = tp / s.size
        auprc += (recall - previous_recall) * tp / np.count_nonzero(selected)
        previous_recall = recall
    return {"auroc": float(auroc), "auprc": float(auprc)}


def _best_threshold(scores, positive):
    values = np.unique(scores)
    thresholds = np.r_[np.nextafter(values[0], -np.inf),
                       (values[:-1] + values[1:]) / 2,
                       np.nextafter(values[-1], np.inf)]

    def cost(threshold):
        predicted = scores >= threshold
        fpr = np.count_nonzero(predicted & ~positive) / np.count_nonzero(~positive)
        fnr = np.count_nonzero(~predicted & positive) / np.count_nonzero(positive)
        return ((fpr + fnr) / 2, fpr, -threshold)

    return float(min(thresholds, key=cost))


def loo_threshold_classification(S, names, comp):
    """Learn a threshold without each compound, then classify every pair containing it.

    Training minimizes balanced error; ties prefer lower false-positive rate. Cross-compound
    pairs are evaluated once for each held-out endpoint compound, each with its own threshold.
    """
    pairs = [(i, j, comp[names[i]], comp[names[j]])
             for i, j in itertools.combinations(range(len(names)), 2)]
    scores = np.asarray([S[i, j] for i, j, _, _ in pairs])
    positive = np.asarray([a == b for _, _, a, b in pairs])
    totals = {"fp": 0, "fn": 0, "positive": 0, "negative": 0}
    thresholds = {}
    for held_out in dict.fromkeys(comp[name] for name in names):
        train = np.asarray([held_out not in (a, b) for _, _, a, b in pairs])
        test = np.asarray([held_out in (a, b) for _, _, a, b in pairs])
        if not np.any(train & positive) or not np.any(train & ~positive):
            raise ValueError("leave-one-compound-out training needs both same and different pairs")
        threshold = _best_threshold(scores[train], positive[train])
        predicted = scores[test] >= threshold
        actual = positive[test]
        totals["fp"] += np.count_nonzero(predicted & ~actual)
        totals["fn"] += np.count_nonzero(~predicted & actual)
        totals["positive"] += np.count_nonzero(actual)
        totals["negative"] += np.count_nonzero(~actual)
        thresholds[str(held_out)] = threshold
    fpr = float(totals["fp"] / totals["negative"])
    fnr = float(totals["fn"] / totals["positive"])
    return {
        "loo_error": (totals["fp"] + totals["fn"]) / (totals["positive"] + totals["negative"]),
        "loo_balanced_error": (fpr + fnr) / 2,
        "rejection_fpr": fpr,
        "false_negative_rate": fnr,
        "loo_thresholds": thresholds,
    }


def cluster_sep(same_vals, cp_mean, comps, iters=10000, seed=1):
    """Compound-level cluster bootstrap: resample the 5 compounds (not the 45 pairs), so the
    dependence induced by each spectrum recurring in 8 different-compound pairs is respected.
    same_vals is indexed by compound; cp_mean maps a compound-pair to its mean cross-spectrum score."""
    rng = np.random.default_rng(seed)
    sv = np.asarray(same_vals)
    K = len(comps)
    seps = np.empty(iters)
    for t in range(iters):
        idx = rng.integers(0, K, K)
        cs = [comps[i] for i in idx]
        dv = [cp_mean[frozenset((cs[a], cs[b]))]
              for a in range(K) for b in range(a + 1, K) if cs[a] != cs[b]]
        seps[t] = sv[idx].mean() - (np.mean(dv) if dv else sv[idx].mean())
    lo, hi = np.percentile(seps, [2.5, 97.5])
    return [float(x) for x in (lo, hi)]


def paired_diff(sa, da, sb, db, iters=10000, seed=0):
    """Paired-difference bootstrap: resample same/diff pairs with IDENTICAL indices for both
    methods, so the separation difference is the right quantity to test between two methods on the
    SAME pairs. Reports Delta = sep_a - sep_b, its 95% CI and P(Delta <= 0)."""
    rng = np.random.default_rng(seed)
    sa, da, sb, db = map(np.asarray, (sa, da, sb, db))
    d = np.empty(iters)
    for t in range(iters):
        si = rng.integers(0, sa.size, sa.size)
        di = rng.integers(0, da.size, da.size)
        d[t] = (sa[si].mean() - da[di].mean()) - (sb[si].mean() - db[di].mean())
    lo, hi = np.percentile(d, [2.5, 97.5])
    return {"delta": float(d.mean()), "ci95": [float(lo), float(hi)], "p_le_0": float((d <= 0).mean())}


def cluster_paired_diff(sa, cpa, sb, cpb, comps, iters=10000, seed=2):
    """Compound-cluster bootstrap of separation(A) - separation(B)."""
    rng = np.random.default_rng(seed)
    sa, sb = map(np.asarray, (sa, sb))
    K = len(comps)
    base_a = [cpa[frozenset((comps[i], comps[j]))]
              for i, j in itertools.combinations(range(K), 2)]
    base_b = [cpb[frozenset((comps[i], comps[j]))]
              for i, j in itertools.combinations(range(K), 2)]
    point = (sa.mean() - np.mean(base_a)) - (sb.mean() - np.mean(base_b))
    delta = np.empty(iters)
    for t in range(iters):
        idx = rng.integers(0, K, K)
        cs = [comps[i] for i in idx]
        da = [cpa[frozenset((cs[i], cs[j]))]
              for i, j in itertools.combinations(range(K), 2) if cs[i] != cs[j]]
        db = [cpb[frozenset((cs[i], cs[j]))]
              for i, j in itertools.combinations(range(K), 2) if cs[i] != cs[j]]
        ma, mb = sa[idx].mean(), sb[idx].mean()
        delta[t] = (ma - (np.mean(da) if da else ma)) - (mb - (np.mean(db) if db else mb))
    lo, hi = np.percentile(delta, [2.5, 97.5])
    return {"delta": float(point), "ci95": [float(lo), float(hi)],
            "p_le_0": float((delta <= 0).mean())}


def main() -> int:
    data_dir = Path("data_13c")
    names = [v for _, a, b in PAIRS for v in (a[0], b[0])]
    comp = {v: c for c, a, b in PAIRS for v in (a[0], b[0])}
    specs = {n: load_peaklist(data_dir / f"{n}.json") for n in names}

    same_idx = [(names.index(a[0]), names.index(b[0])) for _, a, b in PAIRS]
    same_set = {frozenset(p) for p in same_idx}
    comps = [c for c, _, _ in PAIRS]                       # compound label per same-pair
    diff_idx = [(i, j) for i, j in itertools.combinations(range(len(names)), 2)
                if frozenset((i, j)) not in same_set]
    def cpair(i, j):
        return frozenset((comp[names[i]], comp[names[j]]))

    out, sv, dv, cpv = {}, {}, {}, {}
    for m, fn in METHODS.items():
        S = _matrix(specs, names, fn)
        same_vals = [S[i, j] for i, j in same_idx]
        diff_vals = [S[i, j] for i, j in diff_idx]
        sv[m], dv[m] = same_vals, diff_vals
        cp_mean = {}
        for (i, j), val in zip(diff_idx, diff_vals):
            cp_mean.setdefault(cpair(i, j), []).append(val)
        cp_mean = {k: float(np.mean(v)) for k, v in cp_mean.items()}
        cpv[m] = cp_mean
        r = retrieval(S, names, comp)
        rank = ranking_metrics(same_vals, diff_vals)
        loo = loo_threshold_classification(S, names, comp)
        b = bootstrap_sep(same_vals, diff_vals)                       # pair-level (legacy)
        cl = cluster_sep(same_vals, cp_mean, comps)                   # compound-cluster CI
        out[m] = {**r, **rank, **loo, **b, "cluster_ci95": cl}
        print(f"{m:18} top1={r['top1']:.2f} mrr={r['mrr']:.3f} "
              f"AUROC={rank['auroc']:.3f} AUPRC={rank['auprc']:.3f} "
              f"LOOerr={loo['loo_error']:.3f} rejectFPR={loo['rejection_fpr']:.3f} "
              f"sep={b['sep']:.3f} "
              f"pairCI=[{b['ci95'][0]:.3f},{b['ci95'][1]:.3f}] clusterCI=[{cl[0]:.3f},{cl[1]:.3f}]")

    # paired-difference of STCC vs each baseline (identical resample indices per draw)
    paired = {b: paired_diff(sv["lcc_new"], dv["lcc_new"], sv[b], dv[b])
              for b in METHODS if b != "lcc_new"}
    print("\npaired STCC - baseline (Delta, CI, P<=0):")
    for b, v in paired.items():
        print(f"  {b:18} {v['delta']:+.3f} [{v['ci95'][0]:.3f},{v['ci95'][1]:.3f}] P={v['p_le_0']:.3f}")

    cluster_paired = {
        baseline: cluster_paired_diff(
            sv["local_contrast"], cpv["local_contrast"],
            sv[baseline], cpv[baseline], comps,
        )
        for baseline in METHODS if baseline != "local_contrast"
    }
    print("compound-cluster local_contrast - baseline (Delta, CI, P<=0):")
    for baseline, value in cluster_paired.items():
        print(f"  {baseline:18} {value['delta']:+.3f} "
              f"[{value['ci95'][0]:.3f},{value['ci95'][1]:.3f}] P={value['p_le_0']:.3f}")

    Path("results").mkdir(exist_ok=True)
    json.dump({"rows": out, "paired_vs_STCC": paired, "cluster_paired": cluster_paired,
               "n_query": len(names), "n_same": len(same_idx)},
              open("results/retrieval_13c.json", "w"), indent=2)
    print("\nwrote results/retrieval_13c.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
