# Master comparison: both regimes at a glance

All five similarity methods, on both benchmarks. **Separation** = mean(same) − mean(different);
**margin** = worst same − best different (what a single classification threshold actually sees).
Higher is better; all methods self-score exactly 1.00.

| method | `1H-15N` sep | `1H-15N` margin | `1H-13C` sep | `1H-13C` margin |
| --- | --- | --- | --- | --- |
| Bin (Bodis 2009) | 0.29 | 0.24 | 0.69 | 0.20 |
| Bin + 45° | 0.25 | 0.20 | 0.69 | 0.18 |
| Quad-tree (Castillo 2013) | 0.03 | −0.01 | 0.39 | 0.00 |
| Nearest neighbour (Pierens 2012) | 0.04 | 0.03 | 0.10 | 0.04 |
| **LCC (this work)** | **0.75** | **0.71** | **0.81** | **0.37** |

![separation across both regimes](comparison_all.png)

- **Dense protein `1H-15N`** (reference PRL3 vs same-protein titration vs a different protein):
  full detail in [`README.md`](README.md).
- **Sparse small-molecule `1H-13C`** (six compounds each recorded twice, same-compound pairs vs
  different-compound pairs): full detail in [`comparison_13c.md`](comparison_13c.md).

## Reading it

**LCC wins in both regimes** — the only method that keeps same-class similarity high while
pushing different-class similarity down, whether the fingerprint is a dense protein amide region
or a handful of scattered small-molecule crosspeaks.

The bin method is a solid second and, notably, does *much* better on sparse `1H-13C` (0.69) than
on dense `1H-15N` (0.29): with fewer, well-separated peaks its hard bins rarely straddle a peak.
The quad-tree and nearest-neighbour methods **saturate in both regimes** — their margins stay
near zero — because a crowded amide region and a wide small-molecule window both give every peak
a coincidental near neighbour. Their built-in shift tolerance is real (it keeps solvent-shifted
same-compound pairs high, see [`comparison_13c.md`](comparison_13c.md)) but it tolerates
everything, so it does not translate into discrimination.

Numbers: [`comparison_all.json`](comparison_all.json).
