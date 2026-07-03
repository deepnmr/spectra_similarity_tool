# Spectrum similarity for Bruker NMR

This tool implements the bin-based spectrum similarity measure described in:

Lorant Bodis, Alfred Ross, Erno Pretsch, "A novel spectra similarity measure",
Chemometrics and Intelligent Laboratory Systems 85 (2007) 1-8.

`spectrum_similarity.py` handles 1D spectra; `hsqc_similarity.py` extends the
same method to 2D (HSQC and other processed 2D experiments).

## Method

Two spectra $x$ and $y$ are compared inside a common chemical-shift window. Each
processed point $p$ has a chemical shift $\delta(p)$ and an intensity $I_a(p)$
($a \in \{x, y\}$). With the default `clip` baseline, negative intensities are
zeroed first, $I_a(p) \leftarrow \max(I_a(p), 0)$.

Every point carries a weight $w(p)$ equal to its local integration element, so
spectra of different digital resolution integrate consistently:

$$
w_{\text{1D}}(p) = \lvert \Delta\delta(p) \rvert,
\qquad
w_{\text{2D}}(p) = \lvert \Delta\delta_{F1}(p)\rvert \cdot \lvert \Delta\delta_{F2}(p)\rvert .
$$

**Binning.** At resolution $n$ the window is split into equal bins. In 1D, a
window $[\delta_\text{lo}, \delta_\text{hi}]$ of width $W = \delta_\text{hi} -
\delta_\text{lo}$ gives $n$ bins with edges $e_k = \delta_\text{lo} + kW/n$
($k = 0,\dots,n$). In 2D the same is done independently along $F2$ and $F1$,
producing an $n \times n$ grid. The raw integral of bin $k$ and its
normalization to a target integral $N_a$ (`--norm-x` / `--norm-y`, default $1$):

$$
A^{a}_{k}(n) = \sum_{p \,:\, \delta(p) \in \text{bin } k} I_a(p)\, w(p),
\qquad
b^{a}_{k}(n) = N_a \, \frac{A^{a}_{k}(n)}{\sum_{k'} A^{a}_{k'}(n)},
\qquad
\sum_{k} b^{a}_{k}(n) = N_a .
$$

**Similarity index.** Each resolution yields a weighted Jaccard (Ruzicka) index,
clamped to $[0, 1]$:

$$
\mathrm{SI}_n
= \frac{\sum_{k} \min\!\big(b^{x}_{k}, b^{y}_{k}\big)}
       {N_x + N_y - \sum_{k} \min\!\big(b^{x}_{k}, b^{y}_{k}\big)}
= \frac{\sum_{k} \min\!\big(b^{x}_{k}, b^{y}_{k}\big)}
       {\sum_{k} \max\!\big(b^{x}_{k}, b^{y}_{k}\big)},
$$

using $\sum_k \min + \sum_k \max = N_x + N_y$. In 2D the sums run over grid bins
$(i, j)$ with $b^{a}_{ij}(n)$ replacing $b^{a}_{k}(n)$.

The resolution runs $n = 1, \dots, N_\text{bins}$, where the finest bins are
bounded by the minimum bin width $w_\text{min}$:

$$
N_\text{bins}^{\text{1D}} = \left\lfloor \frac{W}{w_\text{min}} \right\rfloor,
\qquad
N_\text{bins}^{\text{2D}} = \min\!\left(
  \left\lfloor \frac{W_{F2}}{w_{\min,F2}} \right\rfloor,\;
  \left\lfloor \frac{W_{F1}}{w_{\min,F1}} \right\rfloor
\right).
$$

**Upper envelope.** At $n = 1$ every pair of normalized spectra collapses to a
single bin, so $\mathrm{SI}_1 = 1$ is uninformative; the envelope $\mathrm{SI}^{*}_n$
smooths the descending $\mathrm{SI}_n$ curve. Anchors are chosen greedily: set
$\mathrm{SI}^{*}_1 = 1$; given an anchor at position $m$, the next anchor is the
largest index $j > m$ that attains $\max_{k > m} \mathrm{SI}_k$. Between
consecutive anchors $\mathrm{SI}^{*}_n$ is linear in $n$, and finally
$\mathrm{SI}^{*}_n = \max(\text{interpolation}, \mathrm{SI}_n)$ so that
$\mathrm{SI}^{*}_n \ge \mathrm{SI}_n$ everywhere.

**Score.** The reported similarity is the mean of the envelope over all
resolutions:

$$
S = \frac{1}{N_\text{bins}} \sum_{n=1}^{N_\text{bins}} \mathrm{SI}^{*}_n
\;\in\; [0, 1].
$$

Identical spectra give $\mathrm{SI}_n = 1$ for every $n$, hence $S = 1$. Because
$S$ depends on $N_\text{bins}$, only compare scores computed with the same
minimum bin width.

## 1D

The input is a processed Bruker 1D spectrum. Pass either an experiment directory,
a `pdata` directory, or a processed `1r` file. The program reads `procs` and `1r`.

## Usage

```bash
python3 spectrum_similarity.py /path/to/exp1 /path/to/exp2
```

Useful options:

```bash
python3 spectrum_similarity.py exp1 exp2 --procno 1 --ppm-min 0 --ppm-max 10
python3 spectrum_similarity.py exp1 exp2 --min-bin-width 0.4 --json
python3 spectrum_similarity.py exp1 exp2 --norm-x 19 --norm-y 19
python3 spectrum_similarity.py exp1 exp2 --plot result.png
python3 spectrum_similarity.py exp1 exp2 --plot --show
```

Defaults:

- `--min-bin-width 0.4`: the value recommended in the paper for 1H NMR.
- `--baseline clip`: negative intensities are set to zero before integration.
- If no ppm range is supplied, the common ppm overlap of the two spectra is used.
- Each spectrum is normalized to total integral 1 unless `--norm-x` and
  `--norm-y` are supplied. For predicted or assigned proton spectra, set these
  to the proton counts when desired.
- `--plot` saves a graph with the two normalized spectra and the `SI_n` /
  `SI*_n` similarity curves. Without a path it writes `spectrum_similarity.png`.
- `--show` displays the same graph interactively when a GUI backend is available.

## 2D HSQC

`hsqc_similarity.py` applies the same bin-method concept to processed Bruker 2D
spectra (HSQC, HMBC, COSY, …). Each spectrum is subdivided into an `n` by `n`
grid of bins whose widths shrink as `n` grows, and the per-grid `SI_n` index and
its `SI*_n` envelope are averaged into one score in `[0, 1]`.

```bash
python3 hsqc_similarity.py /path/to/exp1 /path/to/exp2
python3 hsqc_similarity.py exp1 exp2 --f2-min 0 --f2-max 10 --f1-min 0 --f1-max 160
python3 hsqc_similarity.py exp1 exp2 --min-bin-width-f2 0.1 --min-bin-width-f1 1.0
python3 hsqc_similarity.py exp1 exp2 --rotate 45                 # 1H-13C HSQC
python3 hsqc_similarity.py exp1 exp2 --smooth-f2 0.05 --smooth-f1 0.5
python3 hsqc_similarity.py exp1 exp2 --plot result.png
python3 hsqc_similarity.py exp1 exp2 --json
```

The reader targets processed 2D Bruker data and reconstructs the submatrix
(tile) layout:

```text
experiment/
  pdata/
    1/
      2rr
      procs    # direct dimension F2 (e.g. 1H)
      proc2s   # indirect dimension F1 (e.g. 13C or 15N)
```

It reads `SI`, `OFFSET`, `SW_p`, `SF`, and `XDIM` from `procs`/`proc2s`, plus
`BYTORDP`, `DTYPP`, and `NC_proc` from `procs`, to decode `2rr` and build both
ppm axes.

Defaults:

- `--min-bin-width-f2 0.1`: near the 1H peak linewidth, fine enough to resolve
  small chemical-shift changes without dropping below the digital resolution.
- `--min-bin-width-f1 1.0`: the matching value for the 13C/15N axis, whose ppm
  range is roughly ten times larger. Reduce both together for more discrimination
  between clearly different spectra; the absolute score drops as bins shrink, so
  only compare scores taken at the same bin widths.
- If no ppm range is supplied for a dimension, the common overlap of the two
  spectra is used in that dimension.

### Optional refinements (from the literature)

The core `SI_n` binning above is the exact 2D method of Bodis et al. (2009). Two
paper-grounded options are available and default to **off**, because on protein
`1H-15N` amide HSQC (peaks scattered, no `1H-13C` correlation) they do not help;
they are intended for small-molecule `1H-13C` HSQC.

- `--rotate 45` rotates the binning grid. Bodis et al. (2009) rotate `1H-13C`
  HSQC spectra by 45° so the square bins cross the CH<sub>n</sub> correlation
  diagonal obliquely, which makes a small chemical-shift change less likely to
  move a peak into a neighbouring bin. For `1H-13C` HSQC this improved
  discrimination in the paper; on `1H-15N` data it mainly coarsens the grid.
- `--smooth-f2` / `--smooth-f1` apply a Gaussian shift tolerance (sigma in ppm)
  before binning, a lightweight version of the shift-insensitivity that Castillo
  et al. (2013) build into their tree representation. It raises the score of
  near-matches but also blurs genuine differences.

On the test data (base `1H-15N` HSQC vs the same-protein ligand titration and a
different protein), the plain binning gave the widest separation between
same-protein and different-protein spectra; rotation and smoothing narrowed it.
Both options preserve a self-similarity of exactly 1.

## References

1. L. Bodis, A. Ross, E. Pretsch, *A novel spectra similarity measure*,
   Chemometrics and Intelligent Laboratory Systems 85 (2007) 1-8. — the 1D
   bin method.
2. L. Bodis, A. Ross, J. Bodis, E. Pretsch, *Automatic compatibility tests of
   HSQC NMR spectra with proposed structures of chemical compounds*, Talanta 79
   (2009) 1379-1386. — the 2D (`n × n` bin) extension implemented here,
   including the 45° rotation.
3. A.M. Castillo, L. Uribe, L. Patiny, J. Wist, *Fast and shift-insensitive
   similarity comparisons of NMR using a tree-representation of spectra*,
   Chemometrics and Intelligent Laboratory Systems 127 (2013) 1-6. — the
   shift-insensitive node similarity that motivates the smoothing option.
4. G.K. Pierens, S. Brossi, Z. Yang, D.C. Reutens, V. Vegh, *HSQC spectral based
   similarity matching of compounds using nearest neighbours and a fast discrete
   genetic algorithm*, Journal of Cheminformatics 4 (2012) 25. — a complementary
   peak-list matching approach (not bin-based).

## Bruker assumptions

The reader targets processed 1D Bruker data:

```text
experiment/
  pdata/
    1/
      1r
      procs
```

It uses `SI`, `OFFSET`, `SW_p`, `SF`, `BYTORDP`, `DTYPP`, and `NC_proc` from
`procs` to decode the data and build the ppm axis.

## Test

```bash
python3 -m pytest
```
