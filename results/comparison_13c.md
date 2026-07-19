# Cross-regime check: 1H-13C HSQC (sparse small molecules)

This benchmark complements the dense protein `1H-15N` example with sparse
small-molecule `1H-13C` HSQC peak lists.

## Data and protocol

Public peak lists come from the
[simpleNMR example set](https://github.com/EricHughesABC/simpleNMR). Five
compounds have two genuinely distinct recordings, giving 5 same-compound pairs
and 40 different-compound pairs. A nominal olivetol pair is excluded because its
two files are byte-identical; `bench_13c.py` rejects any such duplicate.

Peak lists are rasterized as a synthetic 2D **stick spectrum** using absolute
intensity, so every method applies its own smoothing or shift-tolerance
processing. This removes the former method-neutral Gaussian pre-blur and avoids
double blurring STCC-family methods. Every comparison uses the same fixed window:
`1H` 0–10 ppm × `13C` 0–165 ppm. STCC uses
`sigma_1H = 0.05`, `sigma_13C = 0.5` ppm.

Regenerate with:

```bash
python3.11 bench_13c.py
python3.11 bench_retrieval.py
```

## Aggregate results

Separation is mean same − mean different; margin is minimum same − maximum
different. All methods self-score exactly 1.00.

| method | mean same | mean different | separation | margin |
| --- | ---: | ---: | ---: | ---: |
| **Local Contrast (experimental)** | **0.8245** | 0.0041 | **0.8203** | **0.4336** |
| STCC (default) | 0.7462 | **0.0016** | 0.7447 | 0.3736 |
| Cosine, un-centred | 0.7466 | 0.0021 | 0.7445 | 0.3727 |
| Bin (Bodis 2009) | 0.7480 | 0.0733 | 0.6747 | 0.2123 |
| Bin + 45° | 0.7675 | 0.0955 | 0.6720 | 0.1943 |
| Quad-tree (Castillo 2013) | 0.9145 | 0.5018 | 0.4127 | 0.0687 |
| Nearest neighbour (Pierens 2012) | 0.9972 | 0.8983 | 0.0990 | 0.0406 |

The local-contrast candidate reuses the STCC-rendered image $G$, computes
$H=\sqrt{\max(G,0)}$ and
$F=H-\operatorname{GaussianBlur}(H,3\sigma)$, then takes the clipped zero-lag
cosine of the two feature images. The factor 3 and square root are fixed; the
method adds no alignment, peak picking, dependency, model, or tuning option.

## Same-compound examples

| pair | Bin | Bin +45° | Tree | NN | Cosine | STCC | Local Contrast |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| menthol, two solvents | 0.734 | 0.807 | 0.967 | 0.988 | 0.421 | 0.421 | **0.467** |
| rotenone, re-measured | 0.323 | 0.344 | 0.839 | 1.000 | 0.389 | 0.388 | **0.727** |
| santonin | 0.973 | 0.976 | 0.993 | 1.000 | 1.000 | 1.000 | 0.999 |
| chartreusin | 0.813 | 0.812 | 0.852 | 0.999 | 0.930 | 0.930 | **0.932** |
| indanone | 0.897 | 0.900 | 0.923 | 1.000 | 0.992 | 0.992 | **0.997** |

Tree and nearest-neighbour scoring keeps shifted positive pairs high, but their
different-compound means are also high. Local Contrast notably lifts the shifted
rotenone pair while retaining a near-zero different-compound mean.

## Held-out evaluation

| method | cluster separation 95% CI | AUROC | AUPRC | LOO error | rejection FPR | top-1 | MRR |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| STCC | 0.508–0.982 | 1.00 | 1.00 | 0.00 | 0.00 | 1.00 | 1.00 |
| Local Contrast | 0.623–0.985 | 1.00 | 1.00 | 0.00 | 0.00 | 1.00 | 1.00 |

The compound-cluster paired Local Contrast − STCC separation difference is
0.0757 (95% CI 0.0014–0.2057, $P(\Delta\leq0)=0.0067$). This supports a
descriptive separation improvement on this small set, but the held-out threshold
error and rejection false-positive rate do not improve: both are already zero.
STCC therefore remains the default and Local Contrast remains experimental.

## Caveats

- The spectra are peak-list sticks, not measured 2D matrices; exact values remain
  benchmark-specific.
- All methods use one fixed wide ppm window. The original tree and nearest-neighbour
  papers used their own parameterization, so these results characterize the repository
  implementations as drop-in global scores.
- Only five compounds contribute positive pairs. The cluster bootstrap respects
  compound dependence, but independent libraries are needed for a promotion decision.

Raw aggregates and pair scores are in
[`comparison_13c.json`](comparison_13c.json); retrieval and held-out metrics are
in [`retrieval_13c.json`](retrieval_13c.json).
