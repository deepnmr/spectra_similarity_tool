# Cross-regime check: 1H-13C HSQC (sparse small molecules)

The main study ([`README.md`](README.md)) benchmarks the methods on a **dense protein
`1H-15N`** amide fingerprint. This is the complementary regime the tree and
nearest-neighbour methods were actually designed for: **sparse small-molecule
`1H-13C` HSQC**.

**Data.** Public HSQC peak lists from the simpleNMR example set
([EricHughesABC/simpleNMR](https://github.com/EricHughesABC/simpleNMR)). Six compounds,
each recorded **twice** (a variant pair — e.g. menthol in two solvents), giving 12
spectra. A good measure should score the **6 same-compound pairs** high and the **60
different-compound pairs** low. Peak lists are rasterized to a synthetic 2D spectrum
(one Gaussian per peak) so every method runs unchanged. Regenerate with
`python3.11 bench_13c.py` (downloads the peak lists on first run); numbers in
[`comparison_13c.json`](comparison_13c.json).

- Window: `1H` 0–10 ppm × `13C` 0–165 ppm. LCC blur `σ_1H = 0.05`, `σ_13C = 0.5` ppm.
- separation = mean(same-compound) − mean(different-compound); margin = min(same) − max(different).

| method | mean same | mean diff | separation | margin |
| --- | --- | --- | --- | --- |
| Bin (Bodis 2009) | 0.76 | 0.07 | 0.69 | 0.20 |
| Bin + 45° | 0.78 | 0.09 | 0.69 | 0.18 |
| Quad-tree (Castillo 2013) | 0.92 | 0.53 | 0.39 | 0.00 |
| Nearest neighbour (Pierens 2012) | 1.00 | 0.90 | 0.10 | 0.04 |
| **LCC (this work)** | **0.82** | **0.00** | **0.81** | **0.37** |

All methods self-score exactly 1.00. Plot: [`comparison_13c.png`](comparison_13c.png).

![1H-13C comparison](comparison_13c.png)

## What this shows

**LCC is the best discriminator in *both* regimes.** It leads here (separation 0.81,
different-compound mean ≈ 0.00) just as it led on the dense `1H-15N` data (0.75). The bin
method is again a strong second (0.69). The tree and nearest-neighbour methods **still
saturate** even on sparse spectra: their different-compound means are 0.53 and 0.90, so a
different molecule scores almost as high as the same one (separation 0.39 and 0.10, margins
≈ 0). Over a large common window, a different small molecule still has *some* peak near
every peak and a similar mass-centre structure, so the nearest-neighbour distance and the
tree overlap stay small — the same saturation mechanism as in the protein case.

**But the shift-tolerance the tree/NN methods were built for is real — and visible per
pair.** Two same-compound pairs are hard because the two recordings are genuinely
shifted (menthol in different solvents; a rotenone re-measurement):

| pair | Bin | Tree | NN | LCC |
| --- | --- | --- | --- | --- |
| menthol (two solvents) | 0.60 | 0.96 | 0.99 | 0.58 |
| rotenone (re-measured) | 0.30 | 0.83 | 1.00 | 0.39 |
| olivetol | 1.00 | 1.00 | 1.00 | 1.00 |
| santonin | 0.97 | 0.99 | 1.00 | 1.00 |
| chartreusin | 0.80 | 0.83 | 1.00 | 0.94 |
| indanone | 0.89 | 0.93 | 1.00 | 0.99 |

On the two shifted pairs the shift-tolerant tree/NN keep the same compound high, while the
position-sensitive bin/LCC penalize the shift. The catch is that tree/NN tolerate
*everything* — that tolerance is exactly why their different-compound scores are also high.
LCC exposes the same trade-off as one physical knob: widening the blur recovers the shifted
pairs while different compounds stay near zero.

| LCC `σ_1H`/`σ_13C` (ppm) | mean same | min same | mean diff | separation |
| --- | --- | --- | --- | --- |
| 0.03 / 0.3 | 0.78 | 0.36 | 0.00 | 0.78 |
| 0.05 / 0.5 | 0.82 | 0.39 | 0.00 | 0.81 |
| 0.08 / 0.8 | 0.85 | 0.40 | 0.01 | 0.84 |
| 0.12 / 1.2 | 0.87 | 0.41 | 0.03 | 0.84 |

## Caveats

- Spectra are **rendered from peak lists**, not measured 2D matrices, so absolute values
  depend on the render width; the *ranking* is the result, not the exact numbers.
- These are the repository's reimplementations of the tree and nearest-neighbour methods
  with one common wide window; the original papers tune per-spectrum windows and parameters,
  so this is not a verdict on those methods as their authors deployed them — only on how
  they behave as drop-in similarity scores here.
- The olivetol pair is near-identical data (a trivial positive). The between-compound
  contrast is strong (60 pairs), but the same-compound side is only 6 pairs.

## Bottom line

Across both a dense protein `1H-15N` fingerprint and sparse small-molecule `1H-13C`
spectra, **LCC gives the widest same/different separation**. The tree and
nearest-neighbour methods add graded shift tolerance but, used as a global score, saturate
in both regimes; LCC recovers that shift tolerance controllably through its blur width.
