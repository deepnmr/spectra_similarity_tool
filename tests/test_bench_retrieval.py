import numpy as np

from bench_retrieval import bootstrap_sep, retrieval


def test_retrieval_ranks_true_partner_first():
    # 4 spectra, 2 compounds (a,a',b,b'); same-compound partner is the closest.
    names = ["a", "a2", "b", "b2"]
    comp = {"a": "A", "a2": "A", "b": "B", "b2": "B"}
    S = np.array([
        [1.0, 0.9, 0.2, 0.1],
        [0.9, 1.0, 0.1, 0.2],
        [0.2, 0.1, 1.0, 0.8],
        [0.1, 0.2, 0.8, 1.0],
    ])
    r = retrieval(S, names, comp)
    assert r["top1"] == 1.0 and r["mrr"] == 1.0


def test_retrieval_penalises_wrong_partner():
    # a's nearest is b (wrong compound); true partner a2 ranks 2nd -> top1=0.5, mrr=0.75.
    names = ["a", "a2", "b", "b2"]
    comp = {"a": "A", "a2": "A", "b": "B", "b2": "B"}
    S = np.array([
        [1.0, 0.5, 0.9, 0.1],   # a closest to b, not a2
        [0.5, 1.0, 0.1, 0.2],
        [0.9, 0.1, 1.0, 0.3],
        [0.1, 0.2, 0.3, 1.0],
    ])
    r = retrieval(S, names, comp)
    assert 0.0 < r["top1"] < 1.0
    assert 0.0 < r["mrr"] < 1.0


def test_bootstrap_ci_brackets_point_estimate():
    same = [0.9, 0.8, 1.0, 0.95, 0.85]
    diff = [0.1, 0.05, 0.2, 0.0, 0.15, 0.1]
    b = bootstrap_sep(same, diff, iters=2000, seed=1)
    lo, hi = b["ci95"]
    assert lo <= b["sep"] <= hi
    assert b["sep"] == float(np.mean(same) - np.mean(diff))
