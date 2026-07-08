# Lineshape Correlation Coefficient (LCC)

A 2D-HSQC similarity measure synthesized from the three methods already in this
repository — the Bodis *et al.* bin method, the Castillo *et al.* quad-tree, and the
Pierens *et al.* nearest-neighbour matcher — designed to keep the bin method's
discrimination on dense protein `1H-15N` fingerprints while removing its
shift-brittleness and without the tree/NN saturation. Implemented in
[`../hsqc_lcc.py`](../hsqc_lcc.py).

**One line:** render each spectrum to a common grid, blur it by the physical NMR
linewidth, and take the mean-centred normalized cross-correlation (Pearson / ZNCC)
of the two images at zero lag.

---

## 1. Setup and preprocessing

Two spectra $a \in \{x, y\}$ are compared inside a common window
$\Omega = [\delta^{F2}_\text{lo}, \delta^{F2}_\text{hi}] \times
[\delta^{F1}_\text{lo}, \delta^{F1}_\text{hi}]$
(direct dimension $F2$, e.g. $^1\mathrm{H}$; indirect dimension $F1$, e.g.
$^{15}\mathrm{N}$). Each processed point $p$ has shifts $\delta_{F2}(p),
\delta_{F1}(p)$ and intensity $I_a(p)$.

With the default `clip` baseline negative intensities are zeroed, and every point
is weighted by its local integration area so spectra of different digital resolution
integrate consistently:

$$
w_a(p) \;=\; \max\!\big(I_a(p),\,0\big)\;\cdot\;
\lvert \Delta\delta_{F1}(p)\rvert \;\cdot\; \lvert \Delta\delta_{F2}(p)\rvert .
$$

(This is exactly the area weighting of the bin method. The mass is normalized to
$\sum_p w_a(p) = 1$; the final ZNCC is invariant to this scale, so normalization
only guards against overflow.)

## 2. Rendering to a shared grid

The window is discretized once into a fixed grid shared by **both** spectra: $K$
bins along $F1$ of width $\Delta_1$ (default `step_f1` $=0.10$ ppm) and $L$ bins
along $F2$ of width $\Delta_2$ (default `step_f2` $=0.01$ ppm). Each spectrum is
rendered by summing its weighted points into cells $(k,l)$:

$$
R_a[k,l] \;=\!\!\sum_{p \,:\, \delta(p)\,\in\,\text{cell}(k,l)}\!\! w_a(p),
\qquad k = 1,\dots,K,\;\; l = 1,\dots,L.
$$

Rendering both spectra onto the *same* edges is what makes the two images directly
comparable regardless of their native resolutions.

**Optional intensity compression.** To stop a single dominant amide from dominating
the correlation, the rendered intensities may be raised to a power
$\rho$ (`--p`, default $\rho = 1$, i.e. off):

$$
R_a[k,l] \;\leftarrow\; R_a[k,l]^{\,\rho}.
$$

## 3. Lineshape blur

The discrete image is convolved with a separable Gaussian whose width is the
physical NMR linewidth *plus* the expected chemical-shift drift, one $\sigma$ per
axis:

$$
\sigma_{F2} = \sqrt{\ell_{F2}^{2} + d_{F2}^{2}},
\qquad
\sigma_{F1} = \sqrt{\ell_{F1}^{2} + d_{F1}^{2}},
$$

with linewidth $\ell$ and expected drift $d$ (defaults $\sigma_{F2}=0.03$,
$\sigma_{F1}=0.30$ ppm). The 1-D kernel, with $s = \sigma/\Delta$ the width in
pixels and radius $r = \lceil 3s \rceil$, is

$$
g_\sigma(m) \;=\; \frac{\exp\!\big(-m^{2}/2s^{2}\big)}
{\sum_{m'=-r}^{r}\exp\!\big(-m'^{2}/2s^{2}\big)},
\qquad m = -r,\dots,r,
$$

and the blurred ("lineshape-rendered") image is the separable convolution

$$
G_a \;=\; g_{\sigma_{F1}} \;*_{F1}\; \big(\, g_{\sigma_{F2}} \;*_{F2}\; R_a \,\big).
$$

The continuous Gaussian is precisely what replaces the bin method's *hard* bin
edges with a *soft* response (see §6).

## 4. Mean-centring and the similarity score

Each image is centred by subtracting its global mean

$$
\bar G_a = \frac{1}{KL}\sum_{k,l} G_a[k,l],
\qquad
\widehat G_a[k,l] = G_a[k,l] - \bar G_a,
$$

and the similarity is the mean-centred normalized cross-correlation at **zero lag**
(the Pearson correlation coefficient of the two flattened images), clamped to
$[0,1]$:

$$
\boxed{\;
S(x,y) \;=\; \max\!\left(0,\;
\frac{\displaystyle\sum_{k,l} \widehat G_x[k,l]\,\widehat G_y[k,l]}
{\sqrt{\displaystyle\sum_{k,l}\widehat G_x[k,l]^{2}}\;
 \sqrt{\displaystyle\sum_{k,l}\widehat G_y[k,l]^{2}}}
\right) \;}
$$

No shift search / alignment is performed — this is deliberate (§6).

## 5. Properties

**Self-similarity is exactly 1.** For $y = x$ the numerator and each factor of the
denominator are the same sum $\lVert\widehat G_x\rVert^{2}$, so

$$
S(x,x) = \frac{\lVert\widehat G_x\rVert^{2}}
{\sqrt{\lVert\widehat G_x\rVert^{2}}\,\sqrt{\lVert\widehat G_x\rVert^{2}}} = 1,
$$

and the $\min(1,\cdot)$ clamp in the code absorbs the last-ULP floating-point wobble
so an identical spectrum returns `1.0` bit-exactly.

**Symmetry.** $S(x,y) = S(y,x)$ (the score is an inner product).

**Range.** As a Pearson coefficient $S \in [-1,1]$ before clamping; the clamp maps
anticorrelation to $0$, giving $S \in [0,1]$.

**Degenerate guard.** A flat image ($\widehat G_a \equiv 0$) has zero denominator
and returns $0$; the `_window` step already excludes empty windows upstream.

## 6. Why it beats the existing methods

The three source methods fail in *opposite* ways: the bin method is shift-**brittle**
(a peak straddling a bin edge splits its mass, so a small titration drift drops the
score in a discontinuous jump), while the tree and NN methods are shift-**blind** (in
a dense fingerprint every peak has a near neighbour, so both saturate near 1). LCC
sits between them via three mechanisms:

**(a) Graded shift tolerance.** Model a peak as a point mass; after the blur it is a
Gaussian of width $\sigma$. Two copies of one peak displaced by
$(\Delta_H, \Delta_N)$ contribute to the (centred, normalized) correlation a factor

$$
\exp\!\left(-\frac{\Delta_H^{2}}{4\sigma_{F2}^{2}}
              -\frac{\Delta_N^{2}}{4\sigma_{F1}^{2}}\right),
$$

which is **smooth and monotone** in the drift — no bin-edge discontinuity. A
physical titration drift ($\Delta_H \!\sim\! 0.01$–$0.05$, $\Delta_N \!\sim\! 0.2$–$0.8$
ppm) is a fraction of $\sigma$ and costs little; a random relocation
($\Delta \gg \sigma$) drives the factor toward $0$. The single knob $\sigma$ sets the
tolerance scale.

**(b) Coverage penalty from mean-centring.** Writing the score as a covariance,

$$
S \propto \operatorname{cov}(G_x, G_y)
= \frac{1}{KL}\sum_{k,l}\big(G_x[k,l]-\bar G_x\big)\big(G_y[k,l]-\bar G_y\big),
$$

a cell where one spectrum has a peak (value above its mean) and the other is empty
(value below its mean) contributes a **negative** product and *lowers* $S$. So
intensity that is not co-located is actively penalised — a differently scattered
protein cannot hide behind coincidental proximity the way it does under NN's
one-way nearest-neighbour distance. This mean-centring step is worth $\approx +0.17$
of separation over an un-centred overlap.

**(c) Zero lag on purpose.** Adding a global shift search (FFT/block matching) would
let a different protein slide into registration and re-saturate — the same failure as
tree/NN. Scoring at zero lag keeps the position information that carries the
discrimination.

## 7. Parameters

| symbol | flag | default | meaning |
| --- | --- | --- | --- |
| $\sigma_{F2}$ | `--sigma-f2` | 0.03 ppm | $^1$H linewidth + drift blur width (main $F2$ tolerance) |
| $\sigma_{F1}$ | `--sigma-f1` | 0.30 ppm | $^{15}$N/$^{13}$C linewidth + drift (main tuning lever) |
| $\Delta_2$ | `--step-f2` | 0.01 ppm | $F2$ render pixel width ($\geq 3$ px per $\sigma$; score grid-invariant) |
| $\Delta_1$ | `--step-f1` | 0.10 ppm | $F1$ render pixel width |
| $\rho$ | `--p` | 1.0 | intensity compression exponent (off) |

Fixed by design: `clip` negatives, centring **on**, shift search **off**.

## 8. Benchmark

Reference PRL3 experiment 2 (`1H-15N` HSQC) versus the same protein plus ligand
(experiments 4–14, should score high) and a different protein (OAA 103, should score
low). Window $F2\,6.5$–$10$, $F1\,105$–$130$ ppm. Separation
$\Sigma = \overline{S}_\text{same} - S_\text{diff}$; margin
$= \min_e S_\text{same} - S_\text{diff}$.

| method | self | mean same | different | $\Sigma$ | margin |
| --- | --- | --- | --- | --- | --- |
| Bin (Bodis 2009) | 1.00 | 0.79 | 0.49 | 0.29 | 0.24 |
| Bin + 45° | 1.00 | 0.82 | 0.57 | 0.25 | 0.20 |
| Quad-tree (Castillo 2013) | 1.00 | 0.90 | 0.87 | 0.03 | −0.01 |
| Nearest-neighbour (Pierens 2012) | 1.00 | 0.99 | 0.96 | 0.04 | 0.03 |
| **LCC (this work)** | 1.00 | **0.94** | **0.18** | **0.75** | **0.71** |

LCC separates same-protein from different-protein spectra $\approx 2.6\times$ better
than the previous best, pushing the different protein below *every* same-protein
score. The gain holds across the whole physical blur range ($\Sigma = 0.71$–$0.77$
at $\sigma = 0.02$–$0.04\,/\,0.20$–$0.40$ ppm; still $0.36$ when coarsened to the bin
method's own resolution $0.10\,/\,1.0$), so it is not an artefact of one tuned
parameter. Numbers and plot in [`../results/`](../results/README.md).

## 9. Relationship to the source methods

| ingredient | taken from |
| --- | --- |
| window + negative clip + area weighting; normalization to unit mass | Bodis *et al.* 2009 (bins) |
| physical blur width as a shift-tolerance ($\sigma = \sqrt{\ell^2 + d^2}$) | Castillo *et al.* 2013 (tree shift-insensitivity) |
| a peak/lineshape picture of the spectrum | Pierens *et al.* 2012 (NN peaks) |
| mean-centred zero-lag correlation (the discriminating step) | this work |

## References

1. L. Bodis, A. Ross, J. Bodis, E. Pretsch, *Automatic compatibility tests of HSQC
   NMR spectra with proposed structures of chemical compounds*, Talanta 79 (2009)
   1379–1386. doi:10.1016/j.talanta.2009.06.017.
2. A.M. Castillo, L. Uribe, L. Patiny, J. Wist, *Fast and shift-insensitive
   similarity comparisons of NMR using a tree-representation of spectra*,
   Chemometrics and Intelligent Laboratory Systems 127 (2013) 1–6.
   doi:10.1016/j.chemolab.2013.05.009.
3. G.K. Pierens, S. Brossi, Z. Yang, D.C. Reutens, V. Vegh, *HSQC spectral based
   similarity matching of compounds using nearest neighbours and a fast discrete
   genetic algorithm*, Journal of Cheminformatics 4 (2012) 25.
   doi:10.1186/1758-2946-4-25.
