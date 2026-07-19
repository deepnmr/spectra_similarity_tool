---
title: "Supporting Information — Robust Sameness Scoring Using a Shift-Tolerant Correlation Coefficient"
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
    \usepackage{caption}
    \captionsetup[figure]{labelformat=empty}
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

Software (STCC, experimental Local Contrast and the three reference methods), all benchmark
harnesses, tabulated scores and the test suite are freely available under the MIT license at
https://github.com/deepnmr/spectra_similarity_tool. Both STCC variants are implemented in
`hsqc_lcc.py`.

## Contents

- **S1** Setup and preprocessing
- **S2** Rendering to a shared grid
- **S3** Lineshape blur
- **S4** Mean-centring, the STCC score and experimental Local Contrast
- **S5** Properties (with proofs)
- **S6** Why STCC beats the existing methods
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
positive, CH$_2$ negative), not absence [1] — `clip` would silently delete every CH$_2$ crosspeak. For
such data use `baseline="abs"`, which weights by $|I_a(p)|$ (replacing $\max(I_a,0)$ in Eq. S1)
and preserves the CH$_2$ peaks; all methods and the CLI expose this option. The sparse
$^1$H–$^{13}$C peak-list benchmark already renders with $|$intensity$|$ for exactly this reason
(§S9).

This is exactly the area weighting of the bin method [2]. The mass is normalized to
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
regardless of their native resolutions [3]. **Optional STCC intensity compression:** to stop a single
dominant amide from dominating the correlation, the rendered intensities may be raised to a power
$\rho$ (`--p`, default $\rho = 1$, i.e. off): $R_a[k,l] \leftarrow R_a[k,l]^{\rho}$. Both the
unit-mass normalization and this optional compression are the standard scaling/normalization
operations of quantitative NMR profiling [4,5]. The experimental Local-Contrast method does not
use `--p`; its fixed square-root transform is applied after the lineshape blur (§S4).

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

## S4 Mean-centring, the STCC score and experimental Local Contrast

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

![**Figure S4.** The four STCC steps end to end, on a toy pair. A reference $X$ and a same-compound
copy $Y$ with a small titration drift are each (1) rendered as area-weighted points on a single
**shared grid**, (2) blurred by the physical linewidth (a separable Gaussian), and (3) mean-centred
(red $=$ above the image mean, blue $=$ below). (4) The similarity is the sum of the elementwise
product $\widehat{G}_X\!\cdot\widehat{G}_Y$ (normalized): co-located intensity gives a **positive**
(red) product and drives the score up, so this same-compound pair scores $S = 0.76$. No shift search
is performed — the images are correlated at zero lag.](../results/si_pipeline.png)

### Experimental Local-Contrast feature

The experimental candidate reuses the completed STCC render $G_a$ from Eq. S5. It compresses the
dynamic range, estimates a broad local background with a Gaussian three times wider than the STCC
lineshape blur on each axis, and subtracts that background:

$$
H_a=\sqrt{\max(G_a,0)}, \qquad
B_a=g_{3\sigma_{F1}}*_{F1}\big(g_{3\sigma_{F2}}*_{F2}H_a\big), \qquad
F_a=H_a-B_a. \tag{SLC1}
$$

The Local-Contrast similarity is the un-centred cosine of the two signed feature images at zero
lag, clipped to $[0,1]$:

$$
S_{\mathrm{LC}}(x,y)=\operatorname{clip}\!\left(
\frac{\langle F_x,F_y\rangle}{\lVert F_x\rVert\lVert F_y\rVert},0,1\right). \tag{SLC2}
$$

If either feature norm is zero, the score is 0. The background factor 3 and square-root transform
are fixed by design and are not new tuning options. The candidate preserves the same fixed ppm
window, physical STCC blur, grid and zero-lag comparison; it adds no peak picking, alignment,
dependency or fitted model. Square-root intensity weighting has precedent in spectral-library
search [38], but that precedent is a design rationale rather than direct evidence of NMR performance;
the evidence here is restricted to the two repository benchmarks.

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

**Local-Contrast contract.** For non-degenerate features, Eq. SLC2 gives exact self-similarity,
symmetry and a score in $[0,1]$. Its normalization makes the score invariant to a global intensity
multiplier. The tests additionally verify monotone decline under increasing chemical-shift drift,
related-above-unrelated ordering after broad baseline and noise are added, a score change below
$10^{-3}$ when each render step is halved, and explicit errors for non-finite data or non-positive
finite blur/step values. A zero feature norm returns 0 as stated above.

## S6 Why STCC beats the existing methods

The three source methods fail in *opposite* ways: the bin method [2] is shift-**brittle** (a peak
straddling a bin edge splits its mass, so a small titration drift drops the score in a
discontinuous jump), while the tree [6,7] and nearest-neighbour (NN) [8] methods are shift-**blind** (in a
dense fingerprint every peak has a near neighbour, so both saturate near 1). STCC sits between them
via three mechanisms.

![**Figure S5.** The two opposite failure modes, on the dense $^1$H–$^{15}$N benchmark (Table S2).
For each method the shaded **band** spans the same-protein similarities (the titration series) and
the red marker is the **different protein**. The bin method is shift-**brittle** — a physical
titration drift erodes its same-protein scores down to 0.73–0.83, the lowest band shown. The tree
and NN methods are shift-**blind** — the different protein sits *inside* (tree) or right at (NN) the
same-protein band, so it is not separated. Only STCC pushes the different protein (0.18) below
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
contribution vs peak displacement $\Delta$. The Gaussian blur makes STCC decay **smoothly and
monotonically** ($\exp(-\Delta^2/4\sigma^2)$, blue), so a physical drift costs little; the bin
method (orange) instead **jumps discontinuously** each time the peak crosses an edge — a small drift
can cause a sudden score drop. This is why the bin method is shift-brittle and STCC is
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
blue $=$ negative). **Left (STCC, mean-centred):** where both have a peak the product is **positive**
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
spectra into a partial false registration — the same over-tolerant saturation as tree/NN. STCC scores
only at $\tau = 0$, keeping the position information that carries the
discrimination.](../results/si_zerolag.png)

## S7 Parameters

**Table S1.** STCC parameters. The two blur widths are the only tuning levers, and both are set by
spectroscopy: the linewidth sets the floor, the expected shift perturbation sets the tolerance.

| symbol | flag | default ($^1$H–$^{15}$N) | meaning |
| --- | --- | :---: | --- |
| $\sigma_{F2}$ | `--sigma-f2` | 0.03 ppm | $^1$H linewidth + drift blur width (main $F2$ tolerance) |
| $\sigma_{F1}$ | `--sigma-f1` | 0.30 ppm | $^{15}$N/$^{13}$C linewidth + drift (main tuning lever) |
| $\Delta_2$ | `--step-f2` | 0.01 ppm | $F2$ render pixel width ($\geq 3$ px per $\sigma$; score grid-invariant) |
| $\Delta_1$ | `--step-f1` | 0.10 ppm | $F1$ render pixel width |
| $\rho$ | `--p` | 1.0 | intensity compression exponent (off) |

For default STCC, centring is **on** and shift search is **off**; `center=False` gives the
un-centred cosine ablation (§S10). Select the experimental candidate with
`--method local-contrast`; its square root and background factor 3 are fixed, and `--p` does not
apply. Negatives are handled by the baseline mode (`clip` default; use `abs` for
multiplicity-edited HSQC — §S1). The CLI default remains `--method lcc`. For the
$^1$H–$^{13}$C benchmark the blur widths are scaled to that regime
($\sigma_{^1\mathrm{H}} = 0.05$, $\sigma_{^{13}\mathrm{C}} = 0.5$ ppm).

![**Figure S10.** The render grid is a discretization, not a tuning parameter. A fixed
same-compound toy pair is scored while the render pixel size is swept: the STCC score is **flat**
while the Gaussian is sampled by $\gtrsim 3$ pixels per $\sigma$ (red zone $=$ under-sampled).
The choice of grid does not affect the result.](../results/si_gridinvariance.png)

## S8 Benchmark 1 — dense protein $^1$H–$^{15}$N

**Design.** Reference spectrum: PRL3 experiment 2, a $^1$H–$^{15}$N HSQC. It is compared with the
**same protein plus ligand** (experiments 4–14, a titration series that should score high — the
chemical-shift-perturbation and fragment-based-screening setting [10–13]) and
with a **different protein** (OAA experiment 103, which should score low). This
same-protein-versus-different-protein discrimination — high for the same fold, low for a different
one — is a component operation of higher-order-structure (HOS) comparability workflows, for which
2D-HSQC fingerprinting is an established atomic-resolution readout (§S13); it does not by itself
exercise the subtle originator-versus-biosimilar comparison, the harder intermediate-score task the
main text leaves to future work. Window $F2\,6.5$–$10 \times F1\,105$–$130$ ppm. Bin widths for the bin method
`min_bin_width_f2 = 0.1`, `min_bin_width_f1 = 1.0`; STCC blur $\sigma_{F2} = 0.03$,
$\sigma_{F1} = 0.30$ ppm. Separation $\Sigma = \overline{S}_\text{same} - S_\text{diff}$; margin
$= \min_e S_\text{same} - S_\text{diff}$.

**Table S2.** Per-comparison similarity, dense protein $^1$H–$^{15}$N benchmark. All methods
self-score exactly 1.

| comparison | Bin (Bodis 2009) | Bin + 45° | Tree (Castillo 2013) | NN (Pierens 2012) | **STCC** |
| --- | :---: | :---: | :---: | :---: | :---: |
| 2 vs 2 (self) | 1.0000 | 1.0000 | 1.0000 | 1.0000 | **1.0000** |
| 2 vs 4 (+B12) | 0.8106 | 0.8385 | 0.9401 | 0.9925 | **0.9662** |
| 2 vs 6 (+NB91) | 0.7739 | 0.7990 | 0.8619 | 0.9981 | **0.9594** |
| 2 vs 8 (+benz) | 0.8274 | 0.8570 | 0.9443 | 0.9877 | **0.9485** |
| 2 vs 10 (+GB001) | 0.8164 | 0.8555 | 0.9179 | 0.9848 | **0.8972** |
| 2 vs 12 (+PF1) | 0.7637 | 0.7912 | 0.8760 | 0.9978 | **0.9558** |
| 2 vs 14 (+PF2) | 0.7320 | 0.7654 | 0.8624 | 0.9851 | **0.8911** |
| **2 vs OAA** | 0.4949 | 0.5681 | 0.8721 | 0.9550 | **0.1815** |
| separation | 0.29 | 0.25 | 0.03 | 0.04 | **0.75** |
| margin | 0.24 | 0.20 | −0.01 | 0.03 | **0.71** |

STCC separates same-protein from different-protein spectra $\approx 2.6\times$ better in separation
(0.75 vs 0.29) than the previous best, pushing the different protein (0.18) below *every*
same-protein score. The old bin
method left the different protein (0.49) close to the worst same-protein point (0.73); STCC opens a
wide gap.

**Expanded benchmark — a second decoy protein.** The benchmark above scores against a single decoy
(OAA). To give the negative side a distribution over distinct folds — and to make the dense test
self-contained — we added an independent PRL3 apo→olsalazine titration (500 µM PRL3(4–159), 23
points from exp 4 to exp 48) and a second decoy protein, the kinase domain **EphB3**, all read as
processed Bruker spectra in the same $^1$H 6.5–10 $\times$ $^{15}$N 105–130 window
(`bench_nhsqc.py`). As the ligand is titrated the chemical-shift perturbation grows and STCC decays
**smoothly** from 0.98 to a plateau near 0.82 (Table S4), the graded shift tolerance of §S6a acting
on real titration data; both decoy proteins fall far below the entire titration (OAA 0.19,
EphB3 0.41), so STCC's margin (0.41) is positive across two distinct folds and the negative side is a
distribution, not one spectrum. Mean-centring still earns the gap (STCC separation 0.57 vs un-centred
cosine 0.46), and the tree and NN measures again saturate (margins 0.12, 0.02). Local Contrast
suppresses the decoys further (0.02 and 0.14), raising separation to 0.67 and margin to 0.52.
STCC's absolute margin is smaller than the six-point single-decoy benchmark (0.71) because this test is harder on both
axes — a longer titration reaching larger shifts (min same 0.82 vs 0.89) and a tougher decoy
(EphB3 0.41 vs OAA 0.19). Local Contrast is reported separately as a post hoc experimental
candidate rather than used to redefine the primary comparison.

**Table S3.** Expanded dense $^1$H–$^{15}$N benchmark: PRL3 apo–olsalazine titration (23
same-protein points) vs two decoy proteins (OAA, EphB3). Separation $= \overline{S}_\text{same} -
\overline{S}_\text{decoy}$; margin $= \min_e S_\text{same} - \max_d S_\text{decoy}$. Every method
self-scores 1.

| method | mean same | min same | OAA | EphB3 |
| --- | :---: | :---: | :---: | :---: |
| **Local Contrast (experimental)** | 0.752 | 0.661 | 0.022 | 0.138 |
| Bin (Bodis 2009) | 0.817 | 0.753 | 0.461 | 0.500 |
| Bin + 45° | 0.865 | 0.810 | 0.521 | 0.563 |
| Quad-tree (Castillo 2013) | 0.946 | 0.912 | 0.793 | 0.749 |
| Nearest neighbour (Pierens 2012) | 0.986 | 0.983 | 0.959 | 0.962 |
| Cosine, un-centred (STCC ablation) | 0.897 | 0.853 | 0.360 | 0.507 |
| **STCC (default)** | 0.871 | 0.820 | 0.195 | 0.414 |

**Table S3 (continued).** Discrimination summary for the expanded dense benchmark.

| method | separation | margin |
| --- | :---: | :---: |
| **Local Contrast (experimental)** | **0.672** | **0.523** |
| Bin (Bodis 2009) | 0.337 | 0.253 |
| Bin + 45° | 0.323 | 0.247 |
| Quad-tree (Castillo 2013) | 0.176 | 0.119 |
| Nearest neighbour (Pierens 2012) | 0.026 | 0.021 |
| Cosine, un-centred (STCC ablation) | 0.464 | 0.346 |
| **STCC (default)** | 0.567 | 0.405 |

**Table S4.** Selected points along the PRL3–olsalazine titration for STCC and experimental Local
Contrast (reference = apo exp 2). Across all 23 titration points, every score stays above both
decoys: STCC has a worst
same-protein score of 0.82 versus worst decoy 0.41; Local Contrast has 0.66 versus 0.14.

| exp | STCC | Local Contrast |
| --- | :---: | :---: |
| 4 | 0.98 | 0.95 |
| 8 | 0.95 | 0.89 |
| 12 | 0.93 | 0.85 |
| 16 | 0.90 | 0.80 |
| 20 | 0.87 | 0.76 |
| 24 | 0.85 | 0.73 |
| 28 | 0.85 | 0.71 |
| 32 | 0.84 | 0.69 |
| 36 | 0.83 | 0.68 |
| 40 | 0.83 | 0.67 |
| 44 | 0.83 | 0.67 |
| 48 | 0.82 | 0.66 |

\clearpage

## S9 Benchmark 2 — sparse small-molecule $^1$H–$^{13}$C

**Design.** Public HSQC peak lists from the simpleNMR example set [14] (EricHughesABC/simpleNMR). Six
nominal compounds, each recorded **twice** (a variant pair — e.g. menthol in two solvents). One
nominal pair, **olivetol**, was excluded: its two source files (`Olivetol.json` and
`Olivetol_A.json`) are **byte-identical** (same MD5, identical peak lists), so it was a
self-comparison scoring 1.00 for every method — an information-free positive that inflated
`mean_same`. The benchmark now uses the **5 genuinely distinct** same-compound pairs (10 spectra)
and asks each measure to score them above the **40 different-compound pairs**. `bench_13c.py`
raises at load time on any same-pair that renders to an identical image, so a duplicate cannot
silently re-enter. Peak lists are rasterized to a synthetic 2D **stick** spectrum, using
$|$intensity$|$ so edited-HSQC CH$_2$ peaks are kept (§S1); each method then applies its own
shift-tolerance processing. This avoids the former method-neutral pre-blur and the resulting double
blur for STCC-family methods.
Window $^1$H 0–10 $\times$ $^{13}$C 0–165 ppm; STCC blur $\sigma_{^1\mathrm{H}} = 0.05$,
$\sigma_{^{13}\mathrm{C}} = 0.5$ ppm.

**Table S5.** Aggregate scores, sparse small-molecule $^1$H–$^{13}$C benchmark (5 same-compound
pairs, 40 different). All methods self-score exactly 1. The un-centred cosine is STCC without
mean-centring (§S6b), included here as the ablation baseline (bug-fix: it is now a reproducible
benchmark row, not an asserted number).

| method | mean same | mean diff | separation | margin |
| --- | :---: | :---: | :---: | :---: |
| **Local Contrast (experimental)** | **0.824** | 0.004 | **0.820** | **0.434** |
| Bin (Bodis 2009) | 0.748 | 0.073 | 0.675 | 0.212 |
| Bin + 45° | 0.768 | 0.096 | 0.672 | 0.194 |
| Quad-tree (Castillo 2013) | 0.914 | 0.502 | 0.413 | 0.069 |
| Nearest neighbour (Pierens 2012) | 0.997 | 0.898 | 0.099 | 0.041 |
| Cosine, un-centred (STCC ablation) | 0.747 | 0.002 | 0.744 | 0.373 |
| **STCC (default)** | 0.746 | 0.002 | 0.745 | 0.374 |

**Table S6.** Per-pair same-compound scores. Two pairs are genuinely shifted (menthol in two
solvents; a rotenone re-measurement); on these the shift-tolerant tree/NN keep the same compound
high while the position-sensitive bin/STCC penalize the shift. The catch is that tree/NN tolerate
*everything*, which is why their different-compound scores (Table S5) are also high. The un-centred
cosine column is identical to STCC to two decimals — in this sparse regime mean-centring changes
nothing (see the ablation, §S10).

| pair | Bin | Tree | NN | Cosine | STCC | Local Contrast |
| --- | :---: | :---: | :---: | :---: | :---: | :---: |
| menthol (two solvents) | 0.734 | 0.967 | 0.988 | 0.421 | 0.421 | **0.467** |
| rotenone (re-measured) | 0.323 | 0.839 | 1.000 | 0.389 | 0.388 | **0.727** |
| santonin | 0.973 | 0.993 | 1.000 | 1.000 | 1.000 | 0.999 |
| chartreusin | 0.813 | 0.852 | 0.999 | 0.930 | 0.930 | 0.932 |
| indanone | 0.897 | 0.923 | 1.000 | 0.992 | 0.992 | 0.997 |

Over a large common window a different small molecule still supplies *some* peak near every peak
and a similar mass-centre structure, so the NN distance and the tree overlap stay small — the same
saturation mechanism as in the protein case.

## S10 Robustness and ablation

**Ablation — mean-centring is the discriminating step, in the dense regime.** Replacing the
mean-centred correlation (Eq. S7) with an un-centred cosine overlap on the *same* rendered/blurred
images drops the dense-protein separation from **0.75 to 0.59**, confirming mean-centring (§S6b) as
the operation that turns overlap into discrimination when the amide fingerprint is crowded. In the
**sparse** $^1$H–$^{13}$C regime, by contrast, the un-centred cosine and centred STCC nearly coincide
(separation **0.74445 vs 0.74466**; Table S5, Table S6): where crosspeaks rarely overlap there is
almost no non-co-located intensity for mean-centring to penalize, so its benefit is specific to the
dense fingerprint that is the hard case this measure targets. The cosine baseline is a real,
reproducible benchmark row in both harnesses (`cosine_similarity` in `hsqc_lcc.py`), not an
asserted number.

**Robustness — the gain is not a tuned parameter.** On the protein benchmark STCC beats the bin
method across the whole physical blur range: separation is $\Sigma = 0.71$–$0.77$ at
$\sigma = 0.02$–$0.04\,/\,0.20$–$0.40$ ppm, and still $\Sigma = 0.36$ (vs the bin method's 0.29)
when the grid is coarsened all the way to the bin method's own resolution ($0.10\,/\,1.0$ ppm).

**Experimental Local Contrast.** On the expanded dense benchmark, the candidate increases
separation/margin from STCC's 0.567/0.405 to 0.672/0.523. On the stick-input sparse benchmark it
increases 0.745/0.374 to 0.820/0.434, driven especially by the shifted rotenone pair
(0.388→0.727). These are composite results for square-root compression *together with* fixed-scale
local background subtraction; the present experiment does not assign the gain to either component
alone. Its representative runtime is 1.34$\times$ STCC. The candidate remains experimental because
the held-out operating metrics tie default STCC (§S13).

**Table S7.** Direct comparison of default STCC and its experimental Local-Contrast extension.

| regime | method | mean same | mean different | separation | margin |
| --- | --- | :---: | :---: | :---: | :---: |
| Dense, 23 same + 2 decoys | STCC | 0.871 | 0.305 | 0.567 | 0.405 |
|  | Local Contrast | 0.752 | 0.080 | **0.672** | **0.523** |
| Sparse, 5 same + 40 different | STCC | 0.746 | 0.002 | 0.745 | 0.374 |
|  | Local Contrast | 0.824 | 0.004 | **0.820** | **0.434** |

The former sparse blur-sweep table is omitted because it was produced from already Gaussian-blurred
peak-list images. The current benchmark starts from sticks so that each method owns its smoothing;
a future sweep must use that same preprocessing before it is compared with the current table.

**Caveats.** (i) In the sparse benchmark, spectra are **stick-rendered from peak lists**, not
measured 2D matrices. Each method owns its smoothing, but exact values remain specific to these
coordinates, intensities and grid; the *ranking* is the result, not a universal calibration.
(ii) The tree and NN methods here are the repository's representative
reimplementations run inside one common wide window; the original papers tune per-spectrum windows
and parameters, so this is not a verdict on those methods as their authors deployed them — only on
how they behave as drop-in global similarity scores. (iii) The between-compound contrast is strong
(40 pairs) but the same-compound side is only 5 pairs, of which two (menthol, rotenone) are the
genuinely-shifted, discriminating tests; a byte-identical sixth pair (olivetol) was excluded
(§S9).

## S11 Relationship to the source methods

**Table S8.** STCC keeps the physically meaningful ingredient of each predecessor. Local Contrast
is an experimental feature transform layered on the completed STCC render.

| ingredient | taken from |
| --- | --- |
| window + negative handling (`clip`/`abs`) + area weighting; normalization to unit mass | Bodis et al. 2009 (bins) |
| physical blur width as shift tolerance ($\sigma = \sqrt{\ell^2 + d^2}$) | Castillo et al. 2013 (tree shift-insensitivity) |
| a peak / lineshape picture of the spectrum | Pierens et al. 2012 (NN peaks) |
| mean-centred zero-lag correlation as the discriminating step | this work (synthesis) |
| square root + fixed-$3\sigma$ local background subtraction | this work (experimental extension) |

The mean-centred zero-lag correlation is itself the zero-mean normalized cross-correlation (ZNCC)
standard in image registration and template matching; the contribution of this work is the
*synthesis* — applying it to physically blurred 2D-HSQC images together with the ingredients above —
not the correlation operation itself. The un-centred cosine of the same rendered/blurred images
(`center=False`) is the natural ablation of the mean-centred-correlation row — it isolates what
mean-centring buys
(§S10) and is exactly the contrast-angle similarity of mass-spectral library search [15,16].
Local Contrast instead uses the signed feature in Eq. SLC1 with an un-centred cosine. The current
benchmark tests that composite transform, not its two ingredients separately.

## S12 Software, reproduction and test suite

**Implementation.** STCC and Local Contrast are pure NumPy [17] (`hsqc_lcc.py`); they read the processed Bruker data (the
`2rr` matrix and its `procs` parameters) directly through the repository's own reader — without a
general-purpose package such as nmrglue [18] — need no peak picking and no external solver. On a
representative dense pair, STCC runs in 12.8 ms and Local Contrast in 17.2 ms
(1.34$\times$; Apple Silicon, single core). The
processed input can come from any standard NMR processing pipeline [19–22]. Figures use
Matplotlib [23].

Local-Contrast JSON records the method, both blur widths, actual grid steps, grid shape, fixed ppm
ranges, `background_factor=3` and `intensity_transform="sqrt"`, so the experimental calculation is
self-describing. Non-finite spectra/ranges and non-positive or non-finite blur/step inputs raise
explicit errors.

**Reproduction.** All benchmark harnesses are in the repository:

- `python3.11 bench.py` — dense protein $^1$H–$^{15}$N benchmark (Table S2), writing
  `results/method_comparison.{csv,json}` and `results/lcc_comparison.png`.
- `python3.11 bench_nhsqc.py` — expanded dense benchmark with two decoy proteins (Tables S3, S4;
  main-text Table 2), reading the `Nhsqc/` PRL3–olsalazine titration and OAA/EphB3 spectra and
  writing `results/nhsqc_dense.json`.
- `python3.11 bench_13c.py` — sparse small-molecule $^1$H–$^{13}$C benchmark (Tables S5, S6),
  downloading the simpleNMR peak lists on first run, rasterizing sticks and writing
  `results/comparison_13c.json`. It
  excludes the byte-identical olivetol pair and raises on any same-pair that renders identically.
- `python3.11 bench_retrieval.py` — top-1/MRR, AUROC/AUPRC, leave-one-compound-out threshold error,
  rejection FPR and compound-cluster bootstrap comparisons on the sparse regime (Table S9),
  writing `results/retrieval_13c.json`.
- `python3.11 results/make_plots.py` — regenerates Figure 2 (and `comparison_13c.png`) from the
  JSON; `python3.11 results/make_fig3.py` — regenerates Figure 3 (bootstrap CIs + ablation) from
  `retrieval_13c.json`; `python3.11 results/make_schematic.py` — regenerates Figure 1;
  `python3.11 results/make_si_figs.py` — regenerates the ten SI method figures (S1–S10).

**Test suite.** A 46-test suite (all passing) checks the properties proved in §S5 and the
Local-Contrast contract: exact self-score, symmetry, $[0,1]$ range, monotonic drift response,
global intensity-scale invariance, baseline/noise related-above-unrelated ordering, half-step grid
invariance and explicit invalid-input errors. Existing regressions cover mean-centring, the `abs`
baseline and the duplicate-pair guard. Run with `python3.11 -m pytest`
(files under `tests/`: `test_hsqc_lcc.py`, `test_hsqc_methods.py`, `test_hsqc_similarity.py`,
`test_bench_13c.py`, `test_bench_nhsqc.py`, `test_bench_retrieval.py`, `test_spectrum_similarity.py`).

## S13 Statistical strength: bootstrap CIs and a retrieval test

These metrics run on the public sparse $^1$H–$^{13}$C data and are reproducible with
`python3.11 bench_retrieval.py` (writes `results/retrieval_13c.json`).

**Bootstrap confidence intervals and operating metrics.** The sparse-regime separation is
resampled $10^4$ times at the level of the five **compounds** (whole compounds rather than the 45
pairs), so the dependence induced by each spectrum recurring in eight different-compound pairs is
respected [24]. The cluster 95% intervals are wide: STCC 0.508–0.982 and Local Contrast
0.623–0.985. A paired compound-cluster bootstrap, which evaluates both methods on the same
resamples, gives Local Contrast − STCC $\Delta=0.0757$ (95% CI 0.0014–0.2057,
$P(\Delta\leq0)=0.0067$).

The promotion decision is based on held-out operation, not separation alone. For each held-out
compound, a threshold minimizing balanced error is learned from the remaining compounds and then
applied to every pair involving the held-out compound. Both STCC and Local Contrast give AUROC and
AUPRC of 1.00, top-1 and MRR of 1.00, zero leave-one-compound-out error and zero rejection FPR and
false-negative rate. The candidate therefore improves descriptive separation on these five
compounds without improving the held-out decision boundary. It remains experimental.

**Table S9.** Current stick-input sparse $^1$H–$^{13}$C evaluation. The CI is the
compound-cluster interval; top-1 and MRR [25] use the 10-spectrum retrieval library.

| method | separation | cluster 95% CI | AUROC | AUPRC |
| --- | :---: | :---: | :---: | :---: |
| **STCC (default)** | 0.745 | 0.508–0.982 | 1.00 | 1.00 |
| **Local Contrast (experimental)** | **0.820** | 0.623–0.985 | 1.00 | 1.00 |

\newpage

**Table S9 (continued).** Held-out and retrieval operating metrics.

| method | LOO error | rejection FPR | top-1 | MRR |
| --- | :---: | :---: | :---: | :---: |
| **STCC (default)** | 0.00 | 0.00 | 1.00 | 1.00 |
| **Local Contrast (experimental)** | 0.00 | 0.00 | 1.00 | 1.00 |

**Retrieval (dereplication).** Cast as a library search — for each of the 10 spectra, rank the
other 9 by similarity and ask whether the true same-compound partner is retrieved first — **every**
method achieves top-1 = 1.00 and MRR = 1.00. This is an honest null result: on a small, clean
library the correct partner is always the single most-similar spectrum, even for the saturating
tree/NN methods whose *different*-compound scores are also high. Retrieval accuracy therefore does
**not** discriminate the methods at this scale; the discriminating quantity is the thresholded
same/different decision (separation and margin), which is what would matter against a large decoy
library or with a rejection threshold. The held-out threshold test above is also tied for STCC and
Local Contrast; a larger, harder library is required to test whether the descriptive gap becomes
an operating advantage. We report the null result rather than omit an experiment
that did not favour the proposed method. A prospective, at-scale dereplication study [26,27] against
a community library (BMRB [28], HMDB [29], NMRShiftDB [30]) is the natural way to turn this into a
discriminating chemical application; STCC provides the scoring primitive, complementary to the
multivariate workflows of NMR-based profiling [31]. A second, equally natural application is the
**higher-order-structure (HOS) comparability** of protein biotherapeutics — biosimilar-versus-
originator, or before/after a formulation or process change — where the fold is confirmed by
matching 2D-HSQC fingerprints and where 2D-NMR is an established, high-precision readout
[32–37]. The dense-protein benchmark (§S8) exercises the distinct-fold end of that task — cleanly
separating fingerprints of different folds and tracking a titration smoothly — but not the subtle
originator-versus-biosimilar discrimination, which falls in the intermediate-score regime and is
left to future work (main text).

## S14 Supporting references

1. W. Willker, D. Leibfritz, R. Kerssebaum, W. Bermel, *Gradient selection in inverse
   heteronuclear correlation spectroscopy*, Magn. Reson. Chem. **1993**, *31*, 287–292.
   doi:10.1002/mrc.1260310315.
2. L. Bodis, A. Ross, J. Bodis, E. Pretsch, *Automatic compatibility tests of HSQC NMR spectra
   with proposed structures of chemical compounds*, Talanta **2009**, *79*, 1379–1386.
   doi:10.1016/j.talanta.2009.06.017.
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
7. H. Samet, *The quadtree and related hierarchical data structures*, ACM Comput. Surv. **1984**,
   *16*, 187–260. doi:10.1145/356924.356930.
8. G. K. Pierens, S. Brossi, Z. Yang, D. C. Reutens, V. Vegh, *HSQC spectral based similarity
   matching of compounds using nearest neighbours and a fast discrete genetic algorithm*,
   J. Cheminform. **2012**, *4*, 25. doi:10.1186/1758-2946-4-25.
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
18. J. J. Helmus, C. P. Jaroniec, *Nmrglue: an open source Python package for the analysis of
    multidimensional NMR data*, J. Biomol. NMR **2013**, *55*, 355–367.
    doi:10.1007/s10858-013-9718-x.
19. F. Delaglio, S. Grzesiek, G. W. Vuister, G. Zhu, J. Pfeifer, A. Bax, *NMRPipe: a multidimensional
    spectral processing system based on UNIX pipes*, J. Biomol. NMR **1995**, *6*, 277–293.
    doi:10.1007/BF00197809.
20. W. F. Vranken, W. Boucher, T. J. Stevens, R. H. Fogh, A. Pajon, M. Llinas, E. L. Ulrich,
    J. L. Markley, J. Ionides, E. D. Laue, *The CCPN data model for NMR spectroscopy: development of
    a software pipeline*, Proteins **2005**, *59*, 687–696. doi:10.1002/prot.20449.
21. S. P. Skinner, R. H. Fogh, W. Boucher, T. J. Ragan, L. G. Mureddu, G. W. Vuister, *CcpNmr
    AnalysisAssign: a flexible platform for integrated NMR analysis*, J. Biomol. NMR **2016**, *66*,
    111–124. doi:10.1007/s10858-016-0060-y.
22. C. Ludwig, U. L. Günther, *MetaboLab — advanced NMR data processing and analysis for
    metabolomics*, BMC Bioinformatics **2011**, *12*, 366. doi:10.1186/1471-2105-12-366.
23. J. D. Hunter, *Matplotlib: a 2D graphics environment*, Comput. Sci. Eng. **2007**, *9*, 90–95.
    doi:10.1109/MCSE.2007.55.
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
37. R. G. Brinson, J. P. Marino, F. Delaglio, L. W. Arbogast, R. M. Evans, A. Kearsley, G. Gingras,
    H. Ghasriani, Y. Aubin, G. K. Pierens, X. Jia, M. Mobli, H. G. Grant, D. W. Keizer, K. Schweimer,
    J. Ståhle, G. Widmalm, E. R. Zartler, C. W. Lawrence, P. N. Reardon, J. R. Cort, P. Xu, F. Ni,
    S. Yanaka, K. Kato, S. R. Parnham, D. Tsao, A. Blomgren, T. Rundlöf, N. Trieloff, P. Schmieder,
    A. Ross, K. Skidmore, K. Chen, D. Keire, D. I. Freedberg, T. Suter-Stahel, G. Wider, G. Ilc,
    J. Plavec, S. A. Bradley, D. M. Baldisseri, M. L. Sforça, A. C. de M. Zeri, J. Y. Wei,
    C. M. Szabo, C. A. Amezcua, J. B. Jordan, M. Wikström, *Enabling adoption of 2D-NMR
    for the higher order structure assessment of monoclonal antibody therapeutics*, mAbs **2019**,
    *11*, 94–105. doi:10.1080/19420862.2018.1544454.
38. C. E. Hart, T. Kind, P. C. Dorrestein, D. Healey, D. Domingo-Fernández, *Weighting
    low-intensity MS/MS ions and m/z frequency for spectral library annotation*, J. Am. Soc. Mass
    Spectrom. **2024**, *35*, 266–274. doi:10.1021/jasms.3c00353.
