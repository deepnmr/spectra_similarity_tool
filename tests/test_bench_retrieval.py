import numpy as np

from bench_retrieval import (
    bootstrap_sep,
    cluster_paired_diff,
    loo_threshold_classification,
    ranking_metrics,
    retrieval,
)


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


def test_three_compound_classification_metrics_are_perfect():
    names = ["a1", "a2", "b1", "b2", "c1", "c2"]
    comp = {name: name[0] for name in names}
    S = np.full((6, 6), 0.1)
    np.fill_diagonal(S, 1.0)
    for i in range(0, 6, 2):
        S[i, i + 1] = S[i + 1, i] = 0.9

    rank = ranking_metrics([0.9] * 3, [0.1] * 12)
    loo = loo_threshold_classification(S, names, comp)

    assert rank == {"auroc": 1.0, "auprc": 1.0}
    assert loo["loo_error"] == loo["loo_balanced_error"] == 0.0
    assert loo["rejection_fpr"] == loo["false_negative_rate"] == 0.0
    assert set(loo["loo_thresholds"]) == {"a", "b", "c"}


def test_identical_score_maps_have_zero_cluster_paired_delta_and_ci():
    comps = ["A", "B", "C"]
    same = [0.9, 0.8, 0.7]
    cross = {
        frozenset(("A", "B")): 0.2,
        frozenset(("A", "C")): 0.3,
        frozenset(("B", "C")): 0.1,
    }

    result = cluster_paired_diff(same, cross, same, cross, comps, iters=100, seed=3)

    assert result["delta"] == 0.0
    assert result["ci95"] == [0.0, 0.0]
