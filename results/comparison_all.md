# Master comparison: both regimes at a glance

Separation is mean same − mean different; margin is worst same − best
different. Higher is better; every method self-scores exactly 1.00.

## Primary method comparison

The dense column is the original PRL3 single-decoy benchmark. The sparse column
uses the current stick-input protocol with 5 same-compound and 40
different-compound pairs.

| method | dense `1H-15N` sep | margin | sparse `1H-13C` sep | margin |
| --- | ---: | ---: | ---: | ---: |
| Bin (Bodis 2009) | 0.292 | 0.237 | 0.675 | 0.212 |
| Bin + 45° | 0.250 | 0.197 | 0.672 | 0.194 |
| Quad-tree (Castillo 2013) | 0.028 | −0.010 | 0.413 | 0.069 |
| Nearest neighbour (Pierens 2012) | 0.036 | 0.030 | 0.099 | 0.041 |
| Cosine, un-centred (STCC ablation) | 0.590 | — | 0.744 | 0.373 |
| ECHOS (Amezcua 2013; box bin, no lineshape) | 0.769 | 0.681 | 0.666 | 0.060 |
| **STCC (default)** | **0.755** | **0.710** | **0.745** | **0.374** |

Among the primary methods, STCC gives the strongest separation on the sparse
regime and ties the ECHOS baseline on the dense one. The two ablations are
complementary: mean-centring matters in the crowded dense fingerprint
(cosine 0.590→0.755) but changes little on sparse sticks (0.744→0.745), whereas
the Gaussian lineshape matters on sparse sticks (ECHOS 0.666/0.060→0.745/0.374)
but not on the dense fingerprint (0.769/0.681 vs 0.755/0.710). ECHOS is the
Pearson correlation of Amezcua & Szabo (J Pharm Sci 2013) on a 5-point box bin —
STCC without the blur.

![primary separation across both regimes](comparison_all.png)

## Experimental Local-Contrast extension

Local Contrast is evaluated on the expanded dense benchmark (23 titration points,
two decoys) and the same sparse stick benchmark:

| benchmark | STCC sep / margin | Local Contrast sep / margin |
| --- | ---: | ---: |
| Dense: 23 same + 2 decoys | 0.5666 / 0.4051 | **0.6719 / 0.5230** |
| Sparse: 5 same + 40 different | 0.7447 / 0.3736 | **0.8203 / 0.4336** |

The candidate improves descriptive separation and margin in both examples. Its
sparse compound-cluster paired gain over STCC is 0.0757 (95% CI
0.0014–0.2057), but AUROC/AUPRC/top-1/MRR remain 1.00 and both methods retain
zero leave-one-compound-out threshold error and rejection FPR. It therefore
remains experimental; `--method lcc` remains the default.

Details:

- Dense original and expanded results: [`README.md`](README.md)
- Sparse aggregate, per-pair and held-out results:
  [`comparison_13c.md`](comparison_13c.md)
- Machine-readable sparse values: [`comparison_13c.json`](comparison_13c.json)
