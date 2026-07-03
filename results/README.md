# Method comparison results

Reference spectrum: **PRL3 experiment 2**, a `1H-15N` HSQC. Compared against the
same protein plus ligand (experiments 4–14, a titration series that should score
high) and against a **different protein** (OAA experiment 103, which should score
low).

- Window: F2 6.5–10 ppm, F1 105–130 ppm
- Bin widths: `min_bin_width_f2 = 0.1`, `min_bin_width_f1 = 1.0`

Raw data: [`method_comparison.csv`](method_comparison.csv),
[`method_comparison.json`](method_comparison.json).

| comparison | bin (Bodis 2009) | bin + 45° | tree (Castillo 2013) | NN (Pierens 2012) |
| --- | --- | --- | --- | --- |
| 2 vs 2 (self) | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| 2 vs 4 (+B12) | 0.8106 | 0.8385 | 0.9401 | 0.9925 |
| 2 vs 6 (+NB91) | 0.7739 | 0.7990 | 0.8619 | 0.9981 |
| 2 vs 8 (+benz) | 0.8274 | 0.8570 | 0.9443 | 0.9877 |
| 2 vs 10 (+GB001) | 0.8164 | 0.8555 | 0.9179 | 0.9848 |
| 2 vs 12 (+PF1) | 0.7637 | 0.7912 | 0.8760 | 0.9978 |
| 2 vs 14 (+PF2) | 0.7320 | 0.7654 | 0.8624 | 0.9851 |
| **2 vs OAA (different protein)** | **0.4949** | 0.5681 | 0.8721 | 0.9550 |
| separation (same − different) | **0.32** | 0.29 | 0.05 | 0.03 |

## Conclusion

All methods self-score exactly 1. On this **dense protein `1H-15N` amide HSQC**
the **bin method separates same-protein from different-protein spectra best**
(≈0.32). The quad-tree and nearest-neighbour methods saturate near 1 (separation
≈0.03–0.05): they were designed for **sparse small-molecule `1H-13C` HSQC** and
for shift insensitivity, so with a crowded amide fingerprint (150+ peaks in one
region) every spectrum has a near neighbour for every peak and a similar
mass-centre structure.

Pick the method by regime: **bins for dense protein fingerprints**,
**tree / nearest-neighbour for sparse small-molecule spectra** where shift
tolerance matters.
