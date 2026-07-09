---
title: "Supporting Information — A Lineshape Correlation Coefficient for Robust Similarity of Two-Dimensional HSQC NMR Spectra"
author: "Donghan Lee"
date: " "
geometry: margin=1in
colorlinks: true
linkcolor: RoyalBlue
urlcolor: RoyalBlue
fontsize: 11pt
header-includes:
  - |
    \usepackage{url}
    \makeatletter
    \g@addto@macro\UrlBreaks{\do\_\do\/\do\-\do\.}
    \makeatother
    \renewcommand{\topfraction}{0.9}
    \renewcommand{\bottomfraction}{0.9}
    \renewcommand{\textfraction}{0.08}
    \renewcommand{\floatpagefraction}{0.7}
    \setcounter{topnumber}{3}
    \setcounter{bottomnumber}{2}
    \setcounter{totalnumber}{5}
---

Korea Basic Science Institute (KBSI), Ochang, Republic of Korea. Correspondence:
kbsi.bionmr@gmail.com

Software (LCC and the three reference methods), both benchmark harnesses, all tabulated scores
and the test suite are freely available under the MIT license at
https://github.com/deepnmr/spectra_similarity_tool. LCC is implemented in `hsqc_lcc.py`.

## Contents

- **S1** Setup and preprocessing
- **S2** Rendering to a shared grid
- **S3** Lineshape blur
- **S4** Mean-centring and the similarity score
- **S5** Properties (with proofs)
- **S6** Why LCC beats the existing methods
- **S7** Parameters
- **S8** Benchmark 1 — dense protein $^1$H–$^{15}$N
- **S9** Benchmark 2 — sparse small-molecule $^1$H–$^{13}$C
- **S10** Robustness and ablation
- **S11** Relationship to the source methods
- **S12** Software, reproduction and test suite
- **S13** Statistical strength: bootstrap CIs and a retrieval test
- **S14** Supporting references

---

## S1 Setup and preprocessing

Two spectra $a \in \{x, y\}$ are compared inside a common window
$\Omega = [\delta^{F2}_\text{lo}, \delta^{F2}_\text{hi}] \times
[\delta^{F1}_\text{lo}, \delta^{F1}_\text{hi}]$
(direct dimension $F2$, e.g. $^1\mathrm{H}$; indirect dimension $F1$, e.g. $^{15}\mathrm{N}$ or
$^{13}\mathrm{C}$). Each processed point $p$ has shifts $\delta_{F2}(p), \delta_{F1}(p)$ and
intensity $I_a(p)$. With the default `clip` baseline, negative intensities are zeroed and every
point is weighted by its local integration area so that spectra of different digital resolution
integrate consistently:

$$
w_a(p) \;=\; \max\!\big(I_a(p),\,0\big)\;\cdot\;
\lvert \Delta\delta_{F1}(p)\rvert \;\cdot\; \lvert \Delta\delta_{F2}(p)\rvert . \tag{S1}
$$

**Negative peaks (edited HSQC).** `clip` is correct for a non-edited spectrum, where negatives are
baseline noise. In a **multiplicity-edited** HSQC, however, the sign encodes multiplicity (CH/CH$_3$
positive, CH$_2$ negative), not absence [2] — `clip` would silently delete every CH$_2$ crosspeak. For
such data use `baseline="abs"`, which weights by $|I_a(p)|$ (replacing $\max(I_a,0)$ in Eq. S1)
and preserves the CH$_2$ peaks; all methods and the CLI expose this option. The sparse
$^1$H–$^{13}$C peak-list benchmark already renders with $|$intensity$|$ for exactly this reason
(§S9).

This is exactly the area weighting of the bin method [1]. The mass is normalized to
$\sum_p w_a(p) = 1$; the final score is invariant to this scale, so normalization only guards
against numerical overflow.

![**Figure S1.** Baseline modes on a multiplicity-edited HSQC cross-section. **Left:** the raw
signal, with a positive CH/CH$_3$ peak and a **negative** CH$_2$ peak (the sign encodes
multiplicity, not absence). **Middle (`clip`, default):** $\max(I,0)$ zeroes negatives — correct
when negatives are baseline noise, but it silently **deletes the CH$_2$ crosspeak**. **Right
(`abs`):** $|I|$ preserves the CH$_2$ peak. Use `clip` for a non-edited spectrum and `abs` for
edited data; the sparse $^1$H–$^{13}$C benchmark renders with $|$intensity$|$ for exactly this
reason (§S9).](../results/si_baseline.png)

## S2 Rendering to a shared grid

The window is discretized once into a fixed grid shared by **both** spectra: $K$ bins along $F1$
of width $\Delta_1$ (default `step_f1` $= 0.10$ ppm) and $L$ bins along $F2$ of width $\Delta_2$
(default `step_f2` $= 0.01$ ppm). Each spectrum is rendered by summing its weighted points into
cells $(k,l)$:

$$
R_a[k,l] \;=\!\!\sum_{p \,:\, \delta(p)\,\in\,\text{cell}(k,l)}\!\! w_a(p),
\qquad k = 1,\dots,K,\;\; l = 1,\dots,L. \tag{S2}
$$

Rendering both spectra onto the *same* edges is what makes the two images directly comparable
regardless of their native resolutions [3]. **Optional intensity compression:** to stop a single
dominant amide from dominating the correlation, the rendered intensities may be raised to a power
$\rho$ (`--p`, default $\rho = 1$, i.e. off): $R_a[k,l] \leftarrow R_a[k,l]^{\rho}$. Both the
unit-mass normalization and this optional compression are the standard scaling/normalization
operations of quantitative NMR profiling [4,5].

![**Figure S2.** Rendering onto a shared grid. The **same** four crosspeaks are recorded by two
instruments at different native digital resolutions (**left**, fine $\Delta$; **middle**, coarse
$\Delta$). Histogramming both onto one **shared** grid with identical cell edges (**right**) places
corresponding peaks in the **same cell**, so the two images become directly comparable regardless of
their native resolutions (Eq. S2).](../results/si_sharedgrid.png)

## S3 Lineshape blur

The discrete image is convolved with a separable Gaussian whose width is the physical NMR
linewidth *plus* the expected chemical-shift drift, one $\sigma$ per axis:

$$
\sigma_{F2} = \sqrt{\ell_{F2}^{2} + d_{F2}^{2}},
\qquad
\sigma_{F1} = \sqrt{\ell_{F1}^{2} + d_{F1}^{2}}, \tag{S3}
$$

with linewidth $\ell$ and expected drift $d$ (defaults $\sigma_{F2} = 0.03$,
$\sigma_{F1} = 0.30$ ppm for $^1$H–$^{15}$N). The 1-D kernel, with $s = \sigma/\Delta$ the width
in pixels and radius $r = \lceil 3s \rceil$, is

$$
g_\sigma(m) \;=\; \frac{\exp\!\big(-m^{2}/2s^{2}\big)}
{\sum_{m'=-r}^{r}\exp\!\big(-m'^{2}/2s^{2}\big)},
\qquad m = -r,\dots,r, \tag{S4}
$$

and the blurred ("lineshape-rendered") image is the separable convolution

$$
G_a \;=\; g_{\sigma_{F1}} \;*_{F1}\; \big(\, g_{\sigma_{F2}} \;*_{F2}\; R_a \,\big). \tag{S5}
$$

The continuous Gaussian is precisely what replaces the bin method's *hard* bin edges with a
*soft* response (see §S6).

![**Figure S3.** The lineshape blur width. **(a)** Each axis blur $\sigma$ combines the physical NMR
linewidth $\ell$ and the expected chemical-shift drift $d$ in **quadrature**,
$\sigma=\sqrt{\ell^2+d^2}$ (Eq. S3). **(b)** Convolving a rendered peak (a spike) with this Gaussian
sets the shift-tolerance scale: a larger $\sigma$ (linewidth *plus* drift) tolerates a larger
displacement than the linewidth alone would.](../results/si_blur.png)

## S4 Mean-centring and the similarity score

Each image is centred by subtracting its global mean

$$
\bar G_a = \frac{1}{KL}\sum_{k,l} G_a[k,l],
\qquad
\widehat G_a[k,l] = G_a[k,l] - \bar G_a, \tag{S6}
$$

and the similarity is the mean-centred normalized cross-correlation at **zero lag** (the Pearson
correlation coefficient of the two flattened images), clamped to $[0,1]$:

$$
\boxed{\;
S(x,y) \;=\; \max\!\left(0,\;
\frac{\displaystyle\sum_{k,l} \widehat G_x[k,l]\,\widehat G_y[k,l]}
{\sqrt{\displaystyle\sum_{k,l}\widehat G_x[k,l]^{2}}\;
 \sqrt{\displaystyle\sum_{k,l}\widehat G_y[k,l]^{2}}}
\right) \;} \tag{S7}
$$

No shift search / alignment is performed — this is deliberate (§S6c).

![**Figure S4.** The four LCC steps end to end, on a toy pair. A reference $X$ and a same-compound
copy $Y$ with a small titration drift are each (1) rendered as area-weighted points on a single
**shared grid**, (2) blurred by the physical linewidth (a separable Gaussian), and (3) mean-centred
(red $=$ above the image mean, blue $=$ below). (4) The similarity is the sum of the elementwise
product $\widehat{G}_X\!\cdot\widehat{G}_Y$ (normalized): co-located intensity gives a **positive**
(red) product and drives the score up, so this same-compound pair scores $S = 0.76$. No shift search
is performed — the images are correlated at zero lag.](../results/si_pipeline.png)

## S5 Properties

**Self-similarity is exactly 1.** For $y = x$ the numerator and each factor of the denominator
are the same sum $\lVert\widehat G_x\rVert^{2}$, so

$$
S(x,x) = \frac{\lVert\widehat G_x\rVert^{2}}
{\sqrt{\lVert\widehat G_x\rVert^{2}}\,\sqrt{\lVert\widehat G_x\rVert^{2}}} = 1, \tag{S8}
$$

and the $\min(1,\cdot)$ clamp in the code absorbs the last-ULP floating-point wobble so an
identical spectrum returns `1.0` bit-exactly.

**Symmetry.** $S(x,y) = S(y,x)$ exactly, since the score is a symmetric inner product of the two
centred images.

**Range.** As a Pearson coefficient $S \in [-1,1]$ before clamping; the clamp maps
anticorrelation to $0$, giving $S \in [0,1]$.

**Monotonicity in drift.** For an isolated peak on a grid large enough that mean-centring is
negligible, $S$ is a strictly decreasing function of the displacement magnitude (Eq. S9); this is
checked numerically in the test suite.

**Degenerate guard.** A flat image ($\widehat G_a \equiv 0$) has zero denominator and returns
$0$; the windowing step already excludes empty windows upstream.

## S6 Why LCC beats the existing methods

The three source methods fail in *opposite* ways: the bin method [1] is shift-**brittle** (a peak
straddling a bin edge splits its mass, so a small titration drift drops the score in a
discontinuous jump), while the tree [6,8] and nearest-neighbour (NN) [7] methods are shift-**blind** (in a
dense fingerprint every peak has a near neighbour, so both saturate near 1). LCC sits between them
via three mechanisms.

![**Figure S5.** The two opposite failure modes, on the dense $^1$H–$^{15}$N benchmark (Table S2).
For each method the shaded **band** spans the same-protein similarities (the titration series) and
the red marker is the **different protein**. The bin method is shift-**brittle** — a physical
titration drift erodes its same-protein scores down to 0.73–0.83, the lowest band shown. The tree
and NN methods are shift-**blind** — the different protein sits *inside* (tree) or right at (NN) the
same-protein band, so it is not separated. Only LCC pushes the different protein (0.18) below
*every* same-protein score, opening a wide clean gap.](../results/si_failuremodes.png)

**(a) Graded shift tolerance.** Model a peak as a point mass; after the blur it is a Gaussian of
width $\sigma$. Two copies of one peak displaced by $(\Delta_{F2}, \Delta_{F1})$ contribute to the
centred, normalized correlation a factor

$$
\exp\!\left(-\frac{\Delta_{F2}^{2}}{4\sigma_{F2}^{2}}
              -\frac{\Delta_{F1}^{2}}{4\sigma_{F1}^{2}}\right), \tag{S9}
$$

which is **smooth and monotone** in the drift — no bin-edge discontinuity. A physical titration
drift ($\Delta_{^1\mathrm{H}} \sim 0.01$–$0.05$, $\Delta_{^{15}\mathrm{N}} \sim 0.2$–$0.8$ ppm) is
a fraction of $\sigma$ and costs little; a random relocation ($\Delta \gg \sigma$) drives the
factor toward $0$. The single knob $\sigma$ sets the tolerance scale (main-text Figure 1b).

![**Figure S6.** Graded vs brittle response to a shift. **(a)** The bin method's edges are hard: as
a peak drifts, a small displacement splits its mass across a bin boundary. **(b)** Similarity
contribution vs peak displacement $\Delta$. The Gaussian blur makes LCC decay **smoothly and
monotonically** ($\exp(-\Delta^2/4\sigma^2)$, blue), so a physical drift costs little; the bin
method (orange) instead **jumps discontinuously** each time the peak crosses an edge — a small drift
can cause a sudden score drop. This is why the bin method is shift-brittle and LCC is
not.](../results/si_shifttol.png)

![**Figure S7.** The single tuning knob. The same-peak contribution $\exp(-\Delta^2/4\sigma^2)$
(Eq. S9) versus peak displacement $\Delta$, for three blur widths $\sigma$. One number sets the
entire tolerance scale: a small $\sigma$ is strict (a modest drift already costs similarity), a
large $\sigma$ is forgiving. The default is chosen so a typical $^{15}$N titration drift (shaded)
costs little while a random relocation drives the contribution toward
zero.](../results/si_sigma.png)

**(b) Coverage penalty from mean-centring.** Writing the score as a covariance,

$$
S \propto \operatorname{cov}(G_x, G_y)
= \frac{1}{KL}\sum_{k,l}\big(G_x[k,l]-\bar G_x\big)\big(G_y[k,l]-\bar G_y\big), \tag{S10}
$$

a cell where one spectrum has a peak (value above its mean) and the other is empty (value below
its mean) contributes a **negative** product and *lowers* $S$. Intensity that is not co-located is
therefore actively penalized — a differently-scattered protein cannot hide behind coincidental
proximity the way it does under NN's one-way nearest-neighbour distance. This mean-centring step is
worth $\approx +0.16$ of separation over an un-centred overlap on the protein benchmark (§S10).

![**Figure S8.** Why mean-centring is the discriminating step (1-D cross-section). Two images
$\widehat{G}_X$ (blue) and $\widehat{G}_Y$ (orange) and their elementwise product (red $=$ positive,
blue $=$ negative). **Left (LCC, mean-centred):** where both have a peak the product is **positive**
(reward); where one has a peak and the other is empty, the empty image sits *below its mean* so the
product is **negative** — non-co-located intensity is actively **penalized**. **Right (un-centred
cosine):** every product is $\ge 0$, so a mismatched peak contributes nothing and cannot be
penalized. Mean-centring is exactly what converts "overlap" into
"discrimination".](../results/si_meancentre.png)

**(c) Zero lag on purpose.** Adding a global shift search (FFT / block matching), as used in
spectral-alignment tools [9], would let a
different protein slide into registration and re-saturate — the same failure as tree/NN. Scoring
at zero lag keeps the position information that carries the discrimination.

![**Figure S9.** Why the score is read at **zero lag**. **(a)** Two genuinely different compounds
(unequal peak spacing, so no rigid shift aligns them all). **(b)** Their normalized correlation as a
function of a rigid lag $\tau$. At zero lag the score is correctly **low** (0.08, blue); a global
shift search would instead report the **maximum** over all lags (0.75, red), sliding the different
spectra into a partial false registration — the same over-tolerant saturation as tree/NN. LCC scores
only at $\tau = 0$, keeping the position information that carries the
discrimination.](../results/si_zerolag.png)

## S7 Parameters

**Table S1.** LCC parameters. The two blur widths are the only tuning levers, and both are set by
spectroscopy: the linewidth sets the floor, the expected shift perturbation sets the tolerance.

| symbol | flag | default ($^1$H–$^{15}$N) | meaning |
| --- | --- | :---: | --- |
| $\sigma_{F2}$ | `--sigma-f2` | 0.03 ppm | $^1$H linewidth + drift blur width (main $F2$ tolerance) |
| $\sigma_{F1}$ | `--sigma-f1` | 0.30 ppm | $^{15}$N/$^{13}$C linewidth + drift (main tuning lever) |
| $\Delta_2$ | `--step-f2` | 0.01 ppm | $F2$ render pixel width ($\geq 3$ px per $\sigma$; score grid-invariant) |
| $\Delta_1$ | `--step-f1` | 0.10 ppm | $F1$ render pixel width |
| $\rho$ | `--p` | 1.0 | intensity compression exponent (off) |

Fixed by design: centring **on**, shift search **off**. Negatives are handled by the baseline
mode (`clip` default; use `abs` for multiplicity-edited HSQC — §S1); `center=False` gives the
un-centred cosine ablation (§S10). For the $^1$H–$^{13}$C benchmark the blur widths are scaled to
that regime ($\sigma_{^1\mathrm{H}} = 0.05$, $\sigma_{^{13}\mathrm{C}} = 0.5$ ppm).

![**Figure S10.** The render grid is a discretization, not a tuning parameter. A fixed
same-compound toy pair is scored while the render pixel size is swept: the LCC score is **flat**
while the Gaussian is sampled by $\gtrsim 3$ pixels per $\sigma$ (red zone $=$ under-sampled).
The choice of grid does not affect the result.](../results/si_gridinvariance.png)

## S8 Benchmark 1 — dense protein $^1$H–$^{15}$N

**Design.** Reference spectrum: PRL3 experiment 2, a $^1$H–$^{15}$N HSQC. It is compared with the
**same protein plus ligand** (experiments 4–14, a titration series that should score high — the
chemical-shift-perturbation and fragment-based-screening setting [10–13]) and
with a **different protein** (OAA experiment 103, which should score low). This
same-protein-versus-different-protein discrimination is the core operation of higher-order-structure
(HOS) comparability assessment, for which 2D-HSQC fingerprinting is an established atomic-resolution
readout (§S13). Window $F2\,6.5$–$10 \times F1\,105$–$130$ ppm. Bin widths for the bin method
`min_bin_width_f2 = 0.1`, `min_bin_width_f1 = 1.0`; LCC blur $\sigma_{F2} = 0.03$,
$\sigma_{F1} = 0.30$ ppm. Separation $\Sigma = \overline{S}_\text{same} - S_\text{diff}$; margin
$= \min_e S_\text{same} - S_\text{diff}$.

**Table S2.** Per-comparison similarity, dense protein $^1$H–$^{15}$N benchmark. All methods
self-score exactly 1.

| comparison | Bin (Bodis 2009) | Bin + 45° | Tree (Castillo 2013) | NN (Pierens 2012) | **LCC** |
| --- | :---: | :---: | :---: | :---: | :---: |
| 2 vs 2 (self) | 1.0000 | 1.0000 | 1.0000 | 1.0000 | **1.0000** |
| 2 vs 4 (+B12) | 0.8106 | 0.8385 | 0.9401 | 0.9925 | **0.9662** |
| 2 vs 6 (+NB91) | 0.7739 | 0.7990 | 0.8619 | 0.9981 | **0.9594** |
| 2 vs 8 (+benz) | 0.8274 | 0.8570 | 0.9443 | 0.9877 | **0.9485** |
| 2 vs 10 (+GB001) | 0.8164 | 0.8555 | 0.9179 | 0.9848 | **0.8972** |
| 2 vs 12 (+PF1) | 0.7637 | 0.7912 | 0.8760 | 0.9978 | **0.9558** |
| 2 vs 14 (+PF2) | 0.7320 | 0.7654 | 0.8624 | 0.9851 | **0.8911** |
| **2 vs OAA (diff. protein)** | 0.4949 | 0.5681 | 0.8721 | 0.9550 | **0.1815** |
| mean same − different (sep.) | 0.29 | 0.25 | 0.03 | 0.04 | **0.75** |
| worst-case margin | 0.24 | 0.20 | −0.01 | 0.03 | **0.71** |

LCC separates same-protein from different-protein spectra $\approx 2.6\times$ better in separation
(0.75 vs 0.29) than the previous best, pushing the different protein (0.18) below *every*
same-protein score. The old bin
method left the different protein (0.49) close to the worst same-protein point (0.73); LCC opens a
wide gap.

**Expanded benchmark — a second decoy protein.** The benchmark above scores against a single decoy
(OAA). To give the negative side a distribution over distinct folds — and to make the dense test
self-contained — we added an independent PRL3 apo→olsalazine titration (500 µM PRL3(4–159), 23
points from exp 4 to exp 48) and a second decoy protein, the kinase domain **EphB3**, all read as
processed Bruker spectra in the same $^1$H 6.5–10 $\times$ $^{15}$N 105–130 window
(`bench_nhsqc.py`). As the ligand is titrated the chemical-shift perturbation grows and LCC decays
**smoothly** from 0.98 to a plateau near 0.82 (Table S9), the graded shift tolerance of §S6a acting
on real titration data; both decoy proteins fall far below the entire titration (OAA 0.19,
EphB3 0.41), so LCC's margin (0.41) is positive across two distinct folds and the negative side is a
distribution, not one spectrum. Mean-centring still earns the gap (LCC separation 0.57 vs un-centred
cosine 0.46), and the tree and NN measures again saturate (margins 0.12, 0.02). The absolute margin
is smaller than the six-point single-decoy benchmark (0.71) because this test is harder on both
axes — a longer titration reaching larger shifts (min same 0.82 vs 0.89) and a tougher decoy
(EphB3 0.41 vs OAA 0.19) — but the method ordering is unchanged.

**Table S8.** Expanded dense $^1$H–$^{15}$N benchmark: PRL3 apo–olsalazine titration (23
same-protein points) vs two decoy proteins (OAA, EphB3). Separation $= \overline{S}_\text{same} -
\overline{S}_\text{decoy}$; margin $= \min_e S_\text{same} - \max_d S_\text{decoy}$. Every method
self-scores 1.

| method | mean same | min same | OAA | EphB3 | separation | margin |
| --- | :---: | :---: | :---: | :---: | :---: | :---: |
| Bin (Bodis 2009) | 0.817 | 0.753 | 0.461 | 0.500 | 0.337 | 0.253 |
| Bin + 45° | 0.865 | 0.810 | 0.521 | 0.563 | 0.323 | 0.247 |
| Quad-tree (Castillo 2013) | 0.946 | 0.912 | 0.793 | 0.749 | 0.176 | 0.119 |
| Nearest neighbour (Pierens 2012) | 0.986 | 0.983 | 0.959 | 0.962 | 0.026 | 0.021 |
| Cosine, un-centred (LCC ablation) | 0.897 | 0.853 | 0.360 | 0.507 | 0.464 | 0.346 |
| **LCC (this work)** | **0.871** | **0.820** | **0.195** | **0.414** | **0.567** | **0.405** |

**Table S9.** LCC along the PRL3–olsalazine titration (reference = apo exp 2): a smooth, near-monotone
decay as ligand — and thus CSP — increases. Every titration point stays above both decoys
(OAA 0.19, EphB3 0.41), so the worst same-protein score (0.82) still clears the worst decoy by 0.41.

| exp | 4 | 8 | 12 | 16 | 20 | 24 | 28 | 32 | 36 | 40 | 44 | 48 |
| --- | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: |
| LCC | 0.98 | 0.95 | 0.93 | 0.90 | 0.87 | 0.85 | 0.85 | 0.84 | 0.83 | 0.83 | 0.83 | 0.82 |

## S9 Benchmark 2 — sparse small-molecule $^1$H–$^{13}$C

**Design.** Public HSQC peak lists from the simpleNMR example set [14] (EricHughesABC/simpleNMR). Six
nominal compounds, each recorded **twice** (a variant pair — e.g. menthol in two solvents). One
nominal pair, **olivetol**, was excluded: its two source files (`Olivetol.json` and
`Olivetol_A.json`) are **byte-identical** (same MD5, identical peak lists), so it was a
self-comparison scoring 1.00 for every method — an information-free positive that inflated
`mean_same`. The benchmark now uses the **5 genuinely distinct** same-compound pairs (10 spectra)
and asks each measure to score them above the **40 different-compound pairs**. `bench_13c.py`
raises at load time on any same-pair that renders to an identical image, so a duplicate cannot
silently re-enter. Peak lists are rasterized to a synthetic 2D spectrum (one Gaussian per peak,
using $|$intensity$|$ so edited-HSQC CH$_2$ peaks are kept — §S1) so every method runs unchanged.
Window $^1$H 0–10 $\times$ $^{13}$C 0–165 ppm; LCC blur $\sigma_{^1\mathrm{H}} = 0.05$,
$\sigma_{^{13}\mathrm{C}} = 0.5$ ppm.

**Table S3.** Aggregate scores, sparse small-molecule $^1$H–$^{13}$C benchmark (5 same-compound
pairs, 40 different). All methods self-score exactly 1. The un-centred cosine is LCC without
mean-centring (§S6b), included here as the ablation baseline (bug-fix: it is now a reproducible
benchmark row, not an asserted number).

| method | mean same | mean diff | separation | margin |
| --- | :---: | :---: | :---: | :---: |
| Bin (Bodis 2009) | 0.712 | 0.067 | 0.646 | 0.195 |
| Bin + 45° | 0.740 | 0.091 | 0.650 | 0.181 |
| Quad-tree (Castillo 2013) | 0.907 | 0.494 | 0.414 | 0.089 |
| Nearest neighbour (Pierens 2012) | 0.997 | 0.896 | 0.101 | 0.043 |
| Cosine, un-centred (LCC ablation) | 0.781 | 0.004 | 0.777 | 0.370 |
| **LCC (this work)** | **0.781** | **0.003** | **0.778** | **0.370** |

**Table S4.** Per-pair same-compound scores. Two pairs are genuinely shifted (menthol in two
solvents; a rotenone re-measurement); on these the shift-tolerant tree/NN keep the same compound
high while the position-sensitive bin/LCC penalize the shift. The catch is that tree/NN tolerate
*everything*, which is why their different-compound scores (Table S3) are also high. The un-centred
cosine column is identical to LCC to two decimals — in this sparse regime mean-centring changes
nothing (see the ablation, §S10).

| pair | Bin | Tree | NN | Cosine | LCC |
| --- | :---: | :---: | :---: | :---: | :---: |
| menthol (two solvents) | 0.60 | 0.96 | 0.99 | 0.58 | 0.58 |
| rotenone (re-measured) | 0.30 | 0.83 | 1.00 | 0.39 | 0.39 |
| santonin | 0.97 | 0.99 | 1.00 | 1.00 | 1.00 |
| chartreusin | 0.80 | 0.83 | 1.00 | 0.94 | 0.94 |
| indanone | 0.89 | 0.93 | 1.00 | 0.99 | 0.99 |

Over a large common window a different small molecule still supplies *some* peak near every peak
and a similar mass-centre structure, so the NN distance and the tree overlap stay small — the same
saturation mechanism as in the protein case.

## S10 Robustness and ablation

**Ablation — mean-centring is the discriminating step, in the dense regime.** Replacing the
mean-centred correlation (Eq. S7) with an un-centred cosine overlap on the *same* rendered/blurred
images drops the dense-protein separation from **0.75 to 0.59**, confirming mean-centring (§S6b) as
the operation that turns overlap into discrimination when the amide fingerprint is crowded. In the
**sparse** $^1$H–$^{13}$C regime, by contrast, the un-centred cosine and centred LCC coincide
(separation **0.777 vs 0.778**; Table S3, Table S4): where crosspeaks rarely overlap there is
almost no non-co-located intensity for mean-centring to penalize, so its benefit is specific to the
dense fingerprint that is the hard case this measure targets. The cosine baseline is a real,
reproducible benchmark row in both harnesses (`cosine_similarity` in `hsqc_lcc.py`), not an
asserted number.

**Robustness — the gain is not a tuned parameter.** On the protein benchmark LCC beats the bin
method across the whole physical blur range: separation is $\Sigma = 0.71$–$0.77$ at
$\sigma = 0.02$–$0.04\,/\,0.20$–$0.40$ ppm, and still $\Sigma = 0.36$ (vs the bin method's 0.29)
when the grid is coarsened all the way to the bin method's own resolution ($0.10\,/\,1.0$ ppm).

**Table S5.** LCC blur sweep, sparse $^1$H–$^{13}$C benchmark — widening the blur recovers the
shifted same-compound pairs while the different-compound mean stays near zero. This is the single
physical knob that recovers the tree/NN shift tolerance without surrendering discrimination.

| $\sigma_{^1\mathrm{H}}\,/\,\sigma_{^{13}\mathrm{C}}$ (ppm) | mean same | min same | mean diff | separation |
| --- | :---: | :---: | :---: | :---: |
| 0.03 / 0.3 | 0.736 | 0.363 | 0.002 | 0.734 |
| 0.05 / 0.5 | 0.781 | 0.392 | 0.003 | 0.778 |
| 0.08 / 0.8 | 0.822 | 0.402 | 0.008 | 0.814 |
| 0.12 / 1.2 | 0.846 | 0.411 | 0.024 | 0.821 |

**Caveats.** (i) In the sparse benchmark, spectra are **rendered from peak lists**, not measured
2D matrices, so absolute values depend on the render width; the *ranking* is the result, not the
exact numbers. (ii) The tree and NN methods here are the repository's representative
reimplementations run inside one common wide window; the original papers tune per-spectrum windows
and parameters, so this is not a verdict on those methods as their authors deployed them — only on
how they behave as drop-in global similarity scores. (iii) The between-compound contrast is strong
(40 pairs) but the same-compound side is only 5 pairs, of which two (menthol, rotenone) are the
genuinely-shifted, discriminating tests; a byte-identical sixth pair (olivetol) was excluded
(§S9).

## S11 Relationship to the source methods

**Table S6.** LCC keeps the physically meaningful ingredient of each predecessor and adds the
mean-centred zero-lag correlation as the new, discriminating step.

| ingredient | taken from |
| --- | --- |
| window + negative handling (`clip`/`abs`) + area weighting; normalization to unit mass | Bodis et al. 2009 (bins) |
| physical blur width as shift tolerance ($\sigma = \sqrt{\ell^2 + d^2}$) | Castillo et al. 2013 (tree shift-insensitivity) |
| a peak / lineshape picture of the spectrum | Pierens et al. 2012 (NN peaks) |
| mean-centred zero-lag correlation (the discriminating step) | this work |

The un-centred cosine of the same rendered/blurred images (`center=False`) is the natural ablation
of the last row — it isolates what mean-centring buys (§S10) and is the contrast-angle similarity
of mass-spectral library search [15,16].

## S12 Software, reproduction and test suite

**Implementation.** LCC is pure NumPy [17] (`hsqc_lcc.py`); it reads the processed Bruker data (the
`2rr` matrix and its `procs` parameters) directly through the repository's own reader — without a
general-purpose package such as nmrglue [19] — needs no peak picking and no external solver, and
runs in about 6 ms per spectrum pair on a $250 \times 350$ grid (Apple Silicon, single core). The
processed input can come from any standard NMR processing pipeline [20–23]. Figures use
Matplotlib [18].

**Reproduction.** Both benchmark harnesses are in the repository:

- `python3.11 bench.py` — dense protein $^1$H–$^{15}$N benchmark (Table S2), writing
  `results/method_comparison.{csv,json}` and `results/lcc_comparison.png`.
- `python3.11 bench_nhsqc.py` — expanded dense benchmark with two decoy proteins (Tables S8, S9;
  main-text Table 2), reading the `Nhsqc/` PRL3–olsalazine titration and OAA/EphB3 spectra and
  writing `results/nhsqc_dense.json`.
- `python3.11 bench_13c.py` — sparse small-molecule $^1$H–$^{13}$C benchmark (Tables S3, S4),
  downloading the simpleNMR peak lists on first run and writing `results/comparison_13c.json`. It
  excludes the byte-identical olivetol pair and raises on any same-pair that renders identically.
- `python3.11 bench_retrieval.py` — dereplication/retrieval metrics (top-1, MRR) and bootstrap 95%
  CIs on the sparse regime (Table S7), writing `results/retrieval_13c.json`.
- `python3.11 results/make_plots.py` — regenerates Figure 2 (and `comparison_13c.png`) from the
  JSON; `python3.11 results/make_fig3.py` — regenerates Figure 3 (bootstrap CIs + ablation) from
  `retrieval_13c.json`; `python3.11 results/make_schematic.py` — regenerates Figure 1;
  `python3.11 results/make_si_figs.py` — regenerates the ten SI method figures (S1–S10).

**Test suite.** A 27-test suite (all passing) checks the properties proved in §S5 (self-similarity
is exactly 1, related spectra score above unrelated ones, monotonicity in drift) plus the three
bug-fix regressions: the un-centred cosine scores an unrelated pair *higher* than mean-centred LCC
(so mean-centring is the discriminating step), the `abs` baseline preserves negative CH$_2$ peaks
that `clip` deletes, and a byte-identical same-pair is rejected. Run with `python3.11 -m pytest`
(files under `tests/`: `test_hsqc_lcc.py`, `test_hsqc_methods.py`, `test_hsqc_similarity.py`,
`test_bench_13c.py`, `test_bench_nhsqc.py`, `test_bench_retrieval.py`, `test_spectrum_similarity.py`).

## S13 Statistical strength: bootstrap CIs and a retrieval test

Both additions run on the public sparse $^1$H–$^{13}$C data and are reproducible with
`python3.11 bench_retrieval.py` (writes `results/retrieval_13c.json`).

**Bootstrap confidence intervals.** The sparse-regime separation is resampled $10^4$ times
(independently resampling the 5 same-compound and 40 different-compound pairs) to give a
non-parametric 95% CI [24] (Table S7). With only five same-compound pairs the CIs are wide: LCC's
interval overlaps the bin method's, so **on this five-compound set LCC's edge over the bin method
is not statistically resolved** — both are strong. LCC's interval is, however, cleanly above the
saturating tree and nearest-neighbour intervals, so its advantage over the two shift-tolerant
methods **is** significant. The honest reading is therefore: in the sparse regime LCC and the bin
method are comparable and both beat tree/NN; LCC's decisive, non-overlapping advantage is in the
**dense** protein regime (margin 0.71). That margin was first measured against a single decoy; the
expanded benchmark (§S8, Tables S8–S9) now scores an independent 23-point PRL3 titration against
**two** decoy proteins (OAA and EphB3) and LCC keeps both below every same-protein point
(margin 0.41), so the dense result no longer rests on one negative spectrum — while the sparse
benchmark and this retrieval test broaden the evidence further with a larger public negative set.

**Table S7.** Sparse $^1$H–$^{13}$C separation with a bootstrap 95% CI, and library-retrieval
metrics. top-1 = fraction of the 10 query spectra whose rank-1 neighbour is the true same-compound
partner; MRR = mean reciprocal rank [25] of that partner.

| method | separation | 95% CI | top-1 | MRR |
| --- | :---: | :---: | :---: | :---: |
| Bin (Bodis 2009) | 0.646 | 0.42–0.84 | 1.00 | 1.00 |
| Bin + 45° | 0.650 | 0.44–0.82 | 1.00 | 1.00 |
| Quad-tree (Castillo 2013) | 0.414 | 0.34–0.49 | 1.00 | 1.00 |
| Nearest neighbour (Pierens 2012) | 0.101 | 0.09–0.11 | 1.00 | 1.00 |
| Cosine, un-centred (LCC ablation) | 0.777 | 0.55–0.98 | 1.00 | 1.00 |
| **LCC (this work)** | **0.778** | **0.55–0.98** | **1.00** | **1.00** |

**Retrieval (dereplication).** Cast as a library search — for each of the 10 spectra, rank the
other 9 by similarity and ask whether the true same-compound partner is retrieved first — **every**
method achieves top-1 = 1.00 and MRR = 1.00. This is an honest null result: on a small, clean
library the correct partner is always the single most-similar spectrum, even for the saturating
tree/NN methods whose *different*-compound scores are also high. Retrieval accuracy therefore does
**not** discriminate the methods at this scale; the discriminating quantity is the thresholded
same/different decision (separation and margin), which is what would matter against a large decoy
library or with a rejection threshold. We report the null result rather than omit an experiment
that did not favour the proposed method. A prospective, at-scale dereplication study [26,27] against
a community library (BMRB [28], HMDB [29], NMRShiftDB [30]) is the natural way to turn this into a
discriminating chemical application; LCC provides the scoring primitive, complementary to the
multivariate workflows of NMR-based profiling [31]. A second, equally natural application is the
**higher-order-structure (HOS) comparability** of protein biotherapeutics — biosimilar-versus-
originator, or before/after a formulation or process change — where the fold is confirmed by
matching 2D-HSQC fingerprints and where 2D-NMR is an established, high-precision readout
[32–37]; the dense-protein benchmark (§S8) is precisely this same/different discrimination task.

## S14 Supporting references

1. L. Bodis, A. Ross, J. Bodis, E. Pretsch, *Automatic compatibility tests of HSQC NMR spectra
   with proposed structures of chemical compounds*, Talanta **2009**, *79*, 1379–1386.
   doi:10.1016/j.talanta.2009.06.017.
2. W. Willker, D. Leibfritz, R. Kerssebaum, W. Bermel, *Gradient selection in inverse
   heteronuclear correlation spectroscopy*, Magn. Reson. Chem. **1993**, *31*, 287–292.
   doi:10.1002/mrc.1260310315.
3. L. Bodis, A. Ross, E. Pretsch, *A novel spectra similarity measure*, Chemometr. Intell. Lab.
   Syst. **2007**, *85*, 1–8. doi:10.1016/j.chemolab.2005.10.002.
4. A. Craig, O. Cloarec, E. Holmes, J. K. Nicholson, J. C. Lindon, *Scaling and normalization
   effects in NMR spectroscopic metabonomic data sets*, Anal. Chem. **2006**, *78*, 2262–2267.
   doi:10.1021/ac0519312.
5. F. Dieterle, A. Ross, G. Schlotterbeck, H. Senn, *Probabilistic quotient normalization as
   robust method to account for dilution of complex biological mixtures*, Anal. Chem. **2006**,
   *78*, 4281–4290. doi:10.1021/ac051632c.
6. A. M. Castillo, L. Uribe, L. Patiny, J. Wist, *Fast and shift-insensitive similarity
   comparisons of NMR using a tree-representation of spectra*, Chemometr. Intell. Lab. Syst.
   **2013**, *127*, 1–6. doi:10.1016/j.chemolab.2013.05.009.
7. G. K. Pierens, S. Brossi, Z. Yang, D. C. Reutens, V. Vegh, *HSQC spectral based similarity
   matching of compounds using nearest neighbours and a fast discrete genetic algorithm*,
   J. Cheminform. **2012**, *4*, 25. doi:10.1186/1758-2946-4-25.
8. H. Samet, *The quadtree and related hierarchical data structures*, ACM Comput. Surv. **1984**,
   *16*, 187–260. doi:10.1145/356924.356930.
9. F. Savorani, G. Tomasi, S. B. Engelsen, *icoshift: a versatile tool for the rapid alignment of
   1D NMR spectra*, J. Magn. Reson. **2010**, *202*, 190–202. doi:10.1016/j.jmr.2009.11.012.
10. M. P. Williamson, *Using chemical shift perturbation to characterise ligand binding*, Prog.
    Nucl. Magn. Reson. Spectrosc. **2013**, *73*, 1–16. doi:10.1016/j.pnmrs.2013.02.001.
11. S. B. Shuker, P. J. Hajduk, R. P. Meadows, S. W. Fesik, *Discovering high-affinity ligands for
    proteins: SAR by NMR*, Science **1996**, *274*, 1531–1534. doi:10.1126/science.274.5292.1531.
12. M. Pellecchia, I. Bertini, D. Cowburn, C. Dalvit, E. Giralt, W. Jahnke, T. L. James,
    S. W. Homans, H. Kessler, C. Luchinat, B. Meyer, H. Oschkinat, J. Peng, H. Schwalbe,
    G. Siegal, *Perspectives on NMR in drug discovery: a technique comes of age*, Nat. Rev. Drug
    Discov. **2008**, *7*, 738–745. doi:10.1038/nrd2606.
13. D. A. Erlanson, S. W. Fesik, R. E. Hubbard, W. Jahnke, H. Jhoti, *Twenty years on: the impact
    of fragments on drug discovery*, Nat. Rev. Drug Discov. **2016**, *15*, 605–619.
    doi:10.1038/nrd.2016.109.
14. E. Hughes, A. M. Kenwright, *SimpleNMR: an interactive graph network approach to aid
    constitutional isomer verification using standard 1D and 2D NMR experiments*, Magn. Reson.
    Chem. **2024**, *62*, 556–565. doi:10.1002/mrc.5441.
15. S. E. Stein, D. R. Scott, *Optimization and testing of mass spectral library search algorithms
    for compound identification*, J. Am. Soc. Mass Spectrom. **1994**, *5*, 859–866.
    doi:10.1016/1044-0305(94)87009-8.
16. K. X. Wan, I. Vidavsky, M. L. Gross, *Comparing similar spectra: from similarity index to
    spectral contrast angle*, J. Am. Soc. Mass Spectrom. **2002**, *13*, 85–88.
    doi:10.1016/S1044-0305(01)00327-0.
17. C. R. Harris, K. J. Millman, S. J. van der Walt, R. Gommers, P. Virtanen, D. Cournapeau,
    E. Wieser, J. Taylor, S. Berg, N. J. Smith, R. Kern, M. Picus, S. Hoyer, M. H. van Kerkwijk,
    M. Brett, A. Haldane, J. F. del Río, M. Wiebe, P. Peterson, P. Gérard-Marchant, K. Sheppard,
    T. Reddy, W. Weckesser, H. Abbasi, C. Gohlke, T. E. Oliphant, *Array programming with NumPy*,
    Nature **2020**, *585*, 357–362. doi:10.1038/s41586-020-2649-2.
18. J. D. Hunter, *Matplotlib: a 2D graphics environment*, Comput. Sci. Eng. **2007**, *9*, 90–95.
    doi:10.1109/MCSE.2007.55.
19. J. J. Helmus, C. P. Jaroniec, *Nmrglue: an open source Python package for the analysis of
    multidimensional NMR data*, J. Biomol. NMR **2013**, *55*, 355–367.
    doi:10.1007/s10858-013-9718-x.
20. F. Delaglio, S. Grzesiek, G. W. Vuister, G. Zhu, J. Pfeifer, A. Bax, *NMRPipe: a multidimensional
    spectral processing system based on UNIX pipes*, J. Biomol. NMR **1995**, *6*, 277–293.
    doi:10.1007/BF00197809.
21. W. F. Vranken, W. Boucher, T. J. Stevens, R. H. Fogh, A. Pajon, M. Llinas, E. L. Ulrich,
    J. L. Markley, J. Ionides, E. D. Laue, *The CCPN data model for NMR spectroscopy: development of
    a software pipeline*, Proteins **2005**, *59*, 687–696. doi:10.1002/prot.20449.
22. S. P. Skinner, R. H. Fogh, W. Boucher, T. J. Ragan, L. G. Mureddu, G. W. Vuister, *CcpNmr
    AnalysisAssign: a flexible platform for integrated NMR analysis*, J. Biomol. NMR **2016**, *66*,
    111–124. doi:10.1007/s10858-016-0060-y.
23. C. Ludwig, U. L. Günther, *MetaboLab — advanced NMR data processing and analysis for
    metabolomics*, BMC Bioinformatics **2011**, *12*, 366. doi:10.1186/1471-2105-12-366.
24. B. Efron, *Bootstrap methods: another look at the jackknife*, Ann. Stat. **1979**, *7*, 1–26.
    doi:10.1214/aos/1176344552.
25. E. M. Voorhees, *The TREC-8 question answering track report*, in *Proceedings of the 8th Text
    REtrieval Conference (TREC-8)*, NIST Special Publication 500-246, Gaithersburg, MD, **1999**,
    pp. 77–82.
26. S. L. Robinette, R. Brüschweiler, F. C. Schroeder, A. S. Edison, *NMR in metabolomics and
    natural products research: two sides of the same coin*, Acc. Chem. Res. **2012**, *45*,
    288–297. doi:10.1021/ar2001606.
27. K. Bingol, R. Brüschweiler, *Multidimensional approaches to NMR-based metabolomics*, Anal.
    Chem. **2014**, *86*, 47–57. doi:10.1021/ac403520j.
28. E. L. Ulrich, H. Akutsu, J. F. Doreleijers, Y. Harano, Y. E. Ioannidis, J. Lin, M. Livny,
    S. Mading, D. Maziuk, Z. Miller, E. Nakatani, C. F. Schulte, D. E. Tolmie, R. K. Wenger,
    H. Yao, J. L. Markley, *BioMagResBank*, Nucleic Acids Res. **2008**, *36*, D402–D408.
    doi:10.1093/nar/gkm957.
29. D. S. Wishart, A. Guo, E. Oler, F. Wang, A. Anjum, H. Peters, R. Dizon, Z. Sayeeda, S. Tian,
    B. L. Lee, M. Berjanskii, R. Mah, M. Yamamoto, J. Jovel, C. Torres-Calzada, M. Hiebert-Giesbrecht,
    V. W. Lui, D. Varshavi, D. Varshavi, D. Allen, D. Arndt, N. Khetarpal, A. Sivakumaran, K. Harford,
    S. Sanford, K. Yee, X. Cao, Z. Budinski, J. Liigand, L. Zhang, J. Zheng, R. Mandal, N. Karu,
    M. Dambrova, H. B. Schiöth, R. Greiner, V. Gautam, *HMDB 5.0: the Human Metabolome Database for
    2022*, Nucleic Acids Res. **2022**, *50*, D622–D631. doi:10.1093/nar/gkab1062.
30. C. Steinbeck, S. Kuhn, *NMRShiftDB — compound identification and structure elucidation support
    through a free community-built web database*, Phytochemistry **2004**, *65*, 2711–2717.
    doi:10.1016/j.phytochem.2004.08.027.
31. B. Worley, R. Powers, *Multivariate analysis in metabolomics*, Curr. Metabolomics **2013**,
    *1*, 92–107. doi:10.2174/2213235X11301010092.
32. S. A. Berkowitz, J. R. Engen, J. R. Mazzeo, G. B. Jones, *Analytical tools for characterizing
    biopharmaceuticals and the implications for biosimilars*, Nat. Rev. Drug Discov. **2012**, *11*,
    527–540. doi:10.1038/nrd3746.
33. J. P. Gabrielson, W. F. Weiss IV, *Technical decision-making with higher order structure data:
    starting a new dialogue*, J. Pharm. Sci. **2015**, *104*, 1240–1245. doi:10.1002/jps.24393.
34. L. W. Arbogast, R. G. Brinson, J. P. Marino, *Mapping monoclonal antibody structure by 2D
    $^{13}$C NMR at natural abundance*, Anal. Chem. **2015**, *87*, 3556–3561.
    doi:10.1021/ac504804m.
35. H. Ghasriani, D. J. Hodgson, R. G. Brinson, I. McEwen, L. F. Buhse, S. Kozlowski, J. P. Marino,
    Y. Aubin, D. A. Keire, *Precision and robustness of 2D-NMR for structure assessment of
    filgrastim biosimilars*, Nat. Biotechnol. **2016**, *34*, 139–141. doi:10.1038/nbt.3474.
36. L. W. Arbogast, F. Delaglio, J. E. Schiel, J. P. Marino, *Multivariate analysis of
    two-dimensional $^1$H,$^{13}$C methyl NMR spectra of monoclonal antibody therapeutics to
    facilitate assessment of higher order structure*, Anal. Chem. **2017**, *89*, 11839–11845.
    doi:10.1021/acs.analchem.7b03571.
37. R. G. Brinson, J. P. Marino, F. Delaglio, L. W. Arbogast, et al., *Enabling adoption of 2D-NMR
    for the higher order structure assessment of monoclonal antibody therapeutics*, mAbs **2019**,
    *11*, 94–105. doi:10.1080/19420862.2018.1544454.
