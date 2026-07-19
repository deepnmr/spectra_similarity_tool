---
title: "Robust Sameness Scoring Using a Shift-Tolerant Correlation Coefficient"
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
---

**NMR Spectroscopy | Research Article**

Korea Basic Science Institute (KBSI), Ochang, Republic of Korea. \* Correspondence:
kbsi.bionmr@gmail.com

---

**Abstract.** Testing whether two two-dimensional $^1$H–X HSQC spectra are the *same* — a
sameness test that underlies molecular fingerprinting, higher-order-structure comparability of
protein biotherapeutics, ligand screening and chemical-shift-perturbation mapping — reduces to a
single quantitative comparison, yet
the established similarity measures each fail in one of two opposite ways: histogram-binning
scores are brittle to the small peak drifts a titration induces, whereas tree- and
nearest-neighbour scores saturate because in any crowded window every peak finds a
coincidental near neighbour. Here we introduce the Shift-Tolerant Correlation Coefficient (STCC),
which renders both spectra onto a single shared grid, blurs each by the physical NMR
linewidth combined with the expected shift drift, and scores the pair by the mean-centred,
zero-lag normalized cross-correlation of the two images. The Gaussian render supplies graded,
physically-tunable shift tolerance while mean-centring actively penalizes intensity that is
not co-located. On two deliberately opposite benchmarks — a dense protein $^1$H–$^{15}$N amide
fingerprint and a sparse small-molecule $^1$H–$^{13}$C set — STCC gives the widest same/different
separation among the primary reference methods (0.75 dense, 0.74 sparse on a $[0,1]$ scale), far
outscoring the saturating tree and nearest-neighbour measures. Its dense-fingerprint advantage is
large and, on an expanded test, keeps both of two decoy proteins below every same-protein titration
point; removing the mean-centring step collapses the original dense-protein separation from 0.75 to
0.59. An experimental, post hoc local-contrast extension raises separation from 0.57 to 0.67 on
the expanded dense test and from 0.74 to 0.82 on the sparse set. The sparse compound-cluster paired
gain is 0.076 (95% CI 0.001–0.206), but both variants already give perfect ranking metrics and zero
leave-one-compound-out threshold error and rejection false-positive rate. STCC therefore remains
the default. It is a two-knob, peak-picking-free, solver-free measure with both knobs set by
spectroscopy rather than by tuning.

**Keywords:** NMR spectroscopy · HSQC · sameness testing · spectral similarity · higher-order structure · chemical-shift perturbation

---

Nuclear magnetic resonance (NMR) spectroscopy is one of chemistry's most powerful probes of
matter, reading out the structure and local environment of individual nuclei non-destructively
and at atomic resolution — from small organic molecules and natural products to metabolites and
folded proteins. Much of that power is exercised not by solving a spectrum *de novo* but by
**comparison**: a great deal of what NMR establishes about an unknown follows from how closely the
spectral fingerprint it produces matches a reference, or from how that fingerprint changes across a
series. This fingerprint-comparison axis — is this the same compound as the catalogued one, is the
protein still folded the same way, has a binding site been perturbed — is at root a **sameness
test**, and it is one of the pillars of
NMR in chemistry, and a two-dimensional $^1$H–X HSQC spectrum is close to an ideal fingerprint for
it: each crosspeak reports the local environment of a single $^1$H–X pair, so the pattern of
crosspeaks encodes molecular structure at atomic resolution.

The comparisons that matter most divide into two families — small molecules and proteins — and
both turn on the same quantitative question. For small molecules, metabolites and natural products,
a $^1$H–$^{13}$C spectrum is matched against a reference, the operation at the heart of
2D-NMR-based dereplication$^{[1-3]}$ and of searches against community libraries such as the BMRB,
NMRShiftDB and HMDB.$^{[4-6]}$ For proteins, an amide $^1$H–$^{15}$N fingerprint is tracked across
a ligand titration to read out binding through chemical-shift perturbation in fragment- and
structure-based drug discovery,$^{[7-10]}$ and — increasingly — is compared between samples to
confirm that a protein's **higher-order structure (HOS)** is preserved. Because a protein's
biological activity, stability and immunogenicity are governed by its folded secondary, tertiary
and quaternary structure and not by its amino-acid sequence alone, establishing that this
higher-order structure is unchanged — between a biosimilar and its originator, across a
manufacturing or formulation change, or on storage — has become a central requirement of
biopharmaceutical characterization and of regulatory comparability.$^{[11,12]}$ Two-dimensional
HSQC is a premier atomic-resolution probe of precisely that structure: two fingerprints coincide
only when the two samples share the same fold, and 2D $^{13}$C- and $^{15}$N-HSQC fingerprinting is
now an established, high-precision method for the HOS comparability of therapeutic proteins and
monoclonal antibodies.$^{[13-15]}$ Resolution is why HSQC, rather than a lower-resolution
technique, carries this burden: methods such as circular dichroism, infrared spectroscopy and
size-exclusion chromatography report structure largely as a global average,$^{[11]}$ whereas an HSQC
fingerprint resolves the local environment of essentially every $^1$H–X site at once and therefore
registers the subtle, localized change in fold or dynamics that an averaged readout can
miss$^{[13,14]}$ — the very differences a sameness test must catch.

Across all of these settings the same primitive is required, and it is deceptively basic: a single
number that is high when two HSQC spectra come from the same species and low when they do not — the
scoring core of a sameness test, exactly the measurement that multivariate and chemometric HOS
workflows already carry out,$^{[16]}$ and the primitive the present work makes robust. Getting that
one number right is deceptively hard, because a sameness test has to resolve differences that are
both small and consequential from the ordinary run-to-run variability of the measurement: a genuine
change in fold, a low-level degradant or a single perturbed residue must register above analytical
noise, solvent- and temperature-dependent peak drift, and the slow drift of the reference material
itself, while an innocuous change must not.$^{[11,12]}$ Posed
precisely, "are these two the same?" is a question of *equivalence* rather than of difference — the
province of the two-one-sided-tests procedure$^{[17]}$ and of interval-based equivalence
testing$^{[18]}$ — and it is the same logic that underpins the regulatory *comparability* of a
biologic before and after a manufacturing change$^{[19]}$ and the
analytical-similarity assessment of a biosimilar against its reference product.$^{[20]}$ Thresholding
that number is what turns a similarity score into a same/different verdict, so the useful measure of
a sameness test is not the raw score but how cleanly it separates same-species from
different-species pairs.

Three families of measure are in current use, and they fail in two opposite ways. The bin
method of Bodis et al. integrates each spectrum over an $n \times n$ grid at a range of
resolutions and reports a weighted Jaccard (Ružička) index,$^{[21,22]}$ in the lineage of the
binning and normalization schemes long used in NMR metabonomics.$^{[23,24]}$ Because its bin
edges are hard, a peak straddling an edge splits its intensity across neighbouring bins; a small
titration drift therefore erodes the score in a discontinuous jump. The measure is
shift-**brittle**. The quad-tree method of Castillo et al. encodes a spectrum as a recursive
centre-of-mass tree and compares trees node-by-node with a shift-tolerant node score,$^{[25]}$
and the peak method of Pierens et al. picks peaks and matches each to its nearest neighbour in
the other spectrum.$^{[26]}$ Both build in genuine shift tolerance — but in a crowded window
every peak has *some* near neighbour and every spectrum a similar mass-centre structure, so both
saturate near unity and retain almost no contrast between related and unrelated spectra. They
are shift-**blind**. A measure fit for real HSQC comparison must tolerate the small physical
drifts a titration or a solvent change induces, yet still penalize intensity that is genuinely
relocated.

We resolve this tension with the Shift-Tolerant Correlation Coefficient (STCC), a measure that keeps
the physically meaningful ingredients of all three predecessors — area-weighted integration and
unit normalization (Bodis), a physical shift tolerance (Castillo), and a lineshape picture of
the spectrum (Pierens) — and adds a mean-centred correlation as the discriminating step. Given
two spectra compared inside a common window $\Omega$ (processed Bruker data from any standard
pipeline$^{[27-31]}$ or peak lists rendered to a grid), STCC proceeds in three steps (Figure 1a).

**Render.** Negative intensities are clipped, each point is weighted by its local integration
area, and *both* spectra are histogrammed onto a **single shared grid**, so they are directly
comparable irrespective of their native digital resolution. The mass is normalized to unit sum;
the final score is invariant to this scale.

**Blur.** Each image is convolved with a separable Gaussian whose width per axis is the physical
linewidth combined with the expected drift, $\sigma = \sqrt{\ell^2 + d^2}$ — an assumed,
effective lineshape (a Gaussian standing in for the true peak shape and deliberately widened by the
expected drift), not a fitted one, so "lineshape" here names this modelled kernel rather than a
measured profile. This replaces the
bin method's hard edges with a continuous response. For an isolated peak on a large grid, two
copies displaced by $(\Delta_2, \Delta_1)$ correlate as
$\exp(-\Delta_2^2/4\sigma_2^2 - \Delta_1^2/4\sigma_1^2)$ — smooth and monotone in the drift
(Figure 1b), so a physical shift costs little while a random relocation ($\Delta \gg \sigma$)
costs almost everything. The blur width, not a bin count, is the tolerance knob, and it is set
by two spectroscopic quantities: the linewidth (the floor) and the expected shift perturbation
(the tolerance).

**Score.** Each image is mean-centred and the similarity is the zero-lag normalized
cross-correlation — the Pearson coefficient — of the two centred images, clamped to $[0,1]$.
This is the mean-centred (Pearson) counterpart of the cosine / contrast-angle similarity long used
in mass-spectral library search$^{[32,33]}$ — equivalently, the zero-mean normalized
cross-correlation (ZNCC) standard in image registration and template matching$^{[34]}$ — applied to
blurred 2D-HSQC images; the un-centred cosine we report as an ablation is exactly that contrast
angle, and mean-centring is what separates STCC from it. Mean-centring is what turns overlap into
discrimination: a cell where one spectrum has a peak (a value above its mean) and the other is
empty (a value below its mean) contributes a **negative** product and lowers the score, so
intensity that is not co-located is actively penalized rather than merely ignored, as it is by
the one-way nearest-neighbour distance. No shift search is performed. Aligning the images — as
in NMR peak-alignment tools such as icoshift$^{[35]}$ — would let a genuinely different spectrum
slide into registration and re-saturate, reintroducing exactly the failure mode of the tree and
nearest-neighbour methods. Scoring at zero lag keeps the position information that carries the
discrimination.

*Experimental local-contrast extension.* We also tested a post hoc candidate that changes only the
feature used for scoring, not the STCC render or its physical shift tolerance. For the existing
Gaussian-blurred image $G$, it forms

$$
H=\sqrt{\max(G,0)}, \qquad
F=H-\operatorname{GaussianBlur}(H,3\sigma),
$$

then takes the clipped zero-lag cosine of $F_x$ and $F_y$. The square root reduces domination by a
few intense peaks, an intensity-weighting rationale used in mass-spectral library search but not
direct evidence of NMR performance;$^{[39]}$ the fixed $3\sigma$ background subtraction removes
broad baseline and common
low-frequency density. The factor three is fixed by design rather than exposed as another tuning
parameter. The candidate adds no dependency, alignment, peak picking or model and preserves the
same fixed-window comparison. Because it was developed after inspecting these benchmarks, we label
it **experimental** and keep STCC as the default unless held-out threshold error or rejection
false-positive rate improves without degrading the other regime.

The measure has the properties one expects of a similarity: self-similarity is exactly one, the
score is exactly symmetric, and it is monotone in the drift of an isolated peak (all verified in
the test suite; Supporting Information). The implementation is pure NumPy,$^{[36]}$ reads
processed Bruker data directly, and requires neither peak picking nor an external solver. On a
representative dense pair, STCC takes 12.8 ms and the experimental local-contrast extension takes
17.2 ms (1.34$\times$; Apple Silicon, single core). A 46-test suite (all passing) checks
self-similarity, symmetry, boundedness, related-beats-unrelated ordering, monotonicity in drift,
intensity-scale and render-grid invariance, explicit invalid-input errors, and the existing
regressions; figures use Matplotlib.$^{[37]}$ A full derivation with proofs is given in the
Supporting Information.

We benchmarked STCC against the three reference methods on two deliberately **opposite** regimes,
so that any measure claiming general utility must succeed on both (Figure 2; Table 1). We
summarize discrimination by the **separation** (mean same-class score − mean different-class
score) and the **margin** (worst same-class − best different-class score, which is what a single
classification threshold — the operating point of a sameness test — actually sees); both are
reported on each method's own $[0,1]$ scale. Every method self-scores exactly 1.00.
Within each benchmark, every method and candidate is evaluated in the same fixed ppm window; no
method-specific automatic overlap or alignment is used.

*Dense protein $^1$H–$^{15}$N.* The reference is the apo spectrum of the phosphatase PRL3,
compared with six spectra of the same protein along a ligand titration (which should score high)
and with a **different** protein (which should score low), inside a
$^1$H 6.5–10 $\times$ $^{15}$N 105–130 ppm window. STCC scores the same protein 0.94 (range
0.89–0.97) and the different protein 0.18, a margin of 0.71 — three times the bin method's margin of 0.24,
the strongest baseline. Mean-centring earns that gap: replacing it with an un-centred cosine
overlap on the *same* images drops the separation from 0.75 to 0.59. This gain is not a tuned
artefact; it holds across the entire physical blur range and even survives coarsening the grid
to the bin method's own resolution (Supporting Information).

*A second decoy protein and an independent titration.* Because that margin is measured against a
single decoy protein, we tested whether STCC's advantage survives a harder, self-contained design:
an independent PRL3 apo→olsalazine titration (23 points, from a 500 µM sample) scored against
**two** distinct decoy folds — OAA and the kinase-domain EphB3 — inside the same window (Table 2;
Supporting Information). As ligand is added the chemical-shift perturbation grows and STCC decays
**smoothly** from 0.98 to a plateau near 0.82 — its graded shift tolerance acting on a real
titration — while both decoys sit far below every titration point (OAA 0.19, EphB3 0.41). STCC keeps
a separation of 0.57 and margin of 0.41: both decoy proteins clear *every*
same-protein score by a comfortable gap, whereas the tree and nearest-neighbour margins are a
knife-edge (0.12 and 0.02, a single noisy titration point from crossing). The un-centred cosine
trails (0.46 / 0.35), so mean-centring still earns the gap against a genuine negative distribution. This test is harder than the single-decoy benchmark on both axes —
a longer titration with larger shifts and a second, tougher decoy — so STCC's dense-regime margin
is not an artefact of one negative spectrum. The experimental
Local-Contrast extension reduces the two decoy scores to 0.02 and 0.14 and increases separation
and margin to 0.67 and 0.52, respectively (Table 2).

*Sparse small-molecule $^1$H–$^{13}$C.* To probe the opposite regime — and to replace the single
negative spectrum above with a large negative set — we used public HSQC peak lists from the
simpleNMR example set$^{[38]}$ for five compounds each recorded twice (one nominal sixth pair was
excluded on inspection as a byte-identical duplicate; Supporting Information), and asked each
method to score the 5 same-compound pairs above the 40 different-compound pairs, inside a
$^1$H 0–10 $\times$ $^{13}$C 0–165 ppm window. The peak lists are rasterized as sticks, after which
each method applies its own shift-tolerance processing; this removes the earlier double blur.
Among the primary methods STCC again leads (mean same 0.746, mean different 0.0016, separation
0.745, margin 0.374). The bin method is far stronger here than on the dense
fingerprint (0.29): with few, well-separated peaks its hard bins rarely straddle a peak. The
tree and nearest-neighbour methods, however, **leave the two classes barely separable in this
regime too** (margins 0.07 and 0.04, despite separations of 0.41 and 0.10): over a wide window a
different molecule still supplies a near neighbour for every peak. Their shift tolerance is
real — per pair, it keeps two shifted recordings of the same compound high (0.83–1.00,
where STCC scores 0.39–0.42) — but it tolerates
*everything*, so it never becomes discrimination. STCC recovers that same tolerance
controllably, through its blur width, without surrendering the between-compound contrast
(Supporting Information). The experimental Local-Contrast extension lifts the rotenone pair from
0.388 to 0.727 while keeping the different-compound mean at 0.0041; overall it reaches separation
0.820 and margin 0.434. Its square-root transform and local background subtraction therefore
improve the descriptive sparse separation without changing peak positions or aligning spectra.

*Statistical strength and held-out operation.* Compound-level resampling keeps all pairs involving
the same molecule together. On the current stick-input sparse benchmark, the cluster 95% interval
for separation is 0.508–0.982 for STCC and 0.623–0.985 for Local Contrast. A paired
compound-cluster bootstrap gives Local Contrast − STCC $\Delta=0.0757$
($[0.0014,0.2057]$, $P(\Delta\leq0)=0.0067$), supporting a descriptive gain on these five
compounds. Ranking and threshold operation provide the promotion check: both variants have
AUROC = AUPRC = top-1 = MRR = 1.00, and both have zero leave-one-compound-out classification error,
false-negative rate and rejection false-positive rate. Thus the candidate does not improve the
held-out decision boundary even though its separation is larger. The expanded dense result is used
as regression evidence, not as an independent held-out estimate because its same-class points are
one autocorrelated titration. Default promotion is reserved for lower held-out threshold error or
rejection FPR without a loss in the other regime; that condition is not met here.

Across both regimes, default STCC beats the saturating tree and nearest-neighbour methods while
retaining graded shift tolerance, and it matches or exceeds the bin method. The post hoc
Local-Contrast extension improves descriptive separation and margin on the expanded dense and
sparse benchmarks, but it does not improve held-out threshold error or rejection FPR and is
therefore not promoted. STCC's two knobs, the per-axis blur widths, are both
physical: the linewidth sets the floor and the expected shift perturbation sets how much drift is
tolerated. Two caveats bound the claim honestly. First, STCC
measures position-and-intensity coincidence, not molecular identity, so a close structural
analogue — a near-identical dereplication candidate, or a biosimilar that differs from its
originator only in a subtle structural detail — will score intermediate; this is the correct
behaviour for a *similarity*, but it is also why the subtle end of higher-order-structure
comparability is a task STCC scores rather than one it yet decides. Second, the three reference methods are our own representative reimplementations
compared inside one common window, and the $^1$H–$^{13}$C spectra are rendered from peak lists;
accordingly the robust **ranking**, not the exact numerical values, is the result, and this is
not a verdict on the original methods as their authors deployed them with per-spectrum windows.

In summary, a single physically-motivated construction — render both spectra onto a shared grid,
blur by the NMR linewidth, and take the mean-centred zero-lag correlation — reconciles shift
tolerance with discrimination, the two requirements that the incumbent binning, tree and
nearest-neighbour measures satisfy only one at a time. Both of STCC's knobs are spectroscopic
quantities rather than fitting parameters, its self-similarity is exact (as for every method
tested), and it needs no peak picking, no alignment search and no external solver. A square-root
local-contrast feature is a promising experimental extension, but the default remains the simpler
validated STCC until it improves held-out operating error rather than descriptive separation alone. We are candid
about scope: this is a **methodological** contribution — a robust, reproducible scoring primitive
for sameness testing, not a new chemical discovery — and although its most decisive number, the dense-protein
margin, was first measured against a single decoy, we now show it holds against two decoy folds on
an independent titration (Table 2), alongside the fully public sparse benchmark and retrieval test.
And the register is deliberately that of a scoring primitive, not a finished assay: STCC supplies the
similarity score a sameness test thresholds, while fixing the equivalence bound and characterizing
the false-same and false-different rates it implies — the inferential half of a formal comparability
test$^{[17,19]}$ — is a task the present benchmarks scope rather than settle.

Made robust, that primitive feeds straight back into the comparisons that motivated it. For small
molecules it is the scoring step for library and dereplication searches and for the quality control
of repeated measurements; for proteins it is the read-out for automated titration tracking and
chemical-shift-perturbation mapping, and the scoring step that higher-order-structure comparability
workflows require.$^{[15,16]}$ On that last application we claim only what the present data establish:
the dense benchmark demonstrates two of the properties an HOS read-out needs — clean separation of
distinct folds (both decoys below every same-protein titration point) and a smooth, near-monotone
response to a genuinely perturbed fingerprint (STCC decaying 0.98→0.82 along the titration rather than
saturating) — but not the harder task at the centre of biosimilar comparability, telling a close
variant from its originator, which falls in the intermediate-score regime noted above and is not
demonstrated here. Characterizing where STCC places true originator/biosimilar or point-variant pairs
is therefore the natural next test on the protein side, as a prospective, at-scale dereplication
study against a community library is on the small-molecule side; for both, STCC now supplies the
tolerant yet discriminating scoring primitive.

More broadly, deciding whether two spectral fingerprints are the same is one of the central operations by
which NMR turns a spectrum into an identification or a structural verdict, and a measure that is at
once tolerant of the small shifts real chemistry induces and unforgiving of genuine relocation makes
that operation more trustworthy wherever it is used. Insofar as NMR remains one of chemistry's
primary windows onto matter — on molecules large and small, in solution, at atomic resolution —
sharpening the way we compare its fingerprints sharpens the readout itself. The full derivation, all
three benchmark protocols, every numerical table with confidence intervals, and the 46-test suite
are provided as Supporting Information; the software is freely available under the MIT license.

## Figures

![**Figure 1.** The Shift-Tolerant Correlation Coefficient. (a) The three-step construction:
area-weighted rendering of both spectra onto one shared grid, a separable Gaussian blur of
physical width $\sigma = \sqrt{\ell^2 + d^2}$ per axis, and the mean-centred zero-lag normalized
cross-correlation of the two images (red/blue = positive/negative contribution to the score).
(b) Graded shift tolerance: two copies of an isolated peak displaced by $\Delta$ correlate as
$\exp(-\Delta^2/4\sigma^2)$, smooth and monotone in the drift, so a physical titration drift
(shaded) costs little while a random relocation costs almost
everything.](../results/lcc_schematic.png)

![**Figure 2.** Same/different separation (mean same-class − mean different-class similarity) for
the primary methods on the two benchmarks: dense protein $^1$H–$^{15}$N (left of each pair) and sparse
small-molecule $^1$H–$^{13}$C (right). STCC leads the primary comparison in both regimes; the tree and nearest-neighbour
methods saturate in both. The un-centred cosine (an ablation of STCC without mean-centring) matches
STCC closely on the sparse regime (0.74 vs 0.74) but falls from 0.75 to 0.59 on the dense fingerprint, isolating
mean-centring as the step that rescues discrimination when peaks
crowd. The experimental Local-Contrast extension is reported in Table 2 and the Supporting
Information rather than mixed into this original dense comparison.](../results/comparison_all.png)

![**Figure 3.** Statistical strength and the mean-centring ablation. **(a)** Sparse
$^1$H–$^{13}$C separation with compound-cluster bootstrap 95% CIs ($10^4$ resamples of the five
compounds), including experimental Local Contrast (purple). The paired Local Contrast − STCC gain
is $\Delta=0.076$ ($[0.001,0.206]$), while both variants have identical perfect ranking and zero
held-out threshold error (Table S9). **(b)** Replacing mean-centred STCC with the un-centred cosine
costs 0.16 of separation on the dense protein fingerprint (0.75→0.59) but only 0.0002 on the sparse
stick set (0.7447→0.7445), showing where global mean-centring matters.](../results/fig3_stats_ablation.png)

**Table 1.** Separation and worst-case margin for the primary methods on both regimes, computed from
unrounded scores. Higher is better; every method self-scores 1.00.

| Method | $^1$H–$^{15}$N sep | margin | $^1$H–$^{13}$C sep | margin |
| --- | :---: | :---: | :---: | :---: |
| Bin (Bodis, 2009) | 0.29 | 0.24 | 0.67 | 0.21 |
| Bin + 45° rotation | 0.25 | 0.20 | 0.67 | 0.19 |
| Quad-tree (Castillo, 2013) | 0.03 | −0.01 | 0.41 | 0.07 |
| Nearest neighbour (Pierens, 2012) | 0.04 | 0.03 | 0.10 | 0.04 |
| Cosine, un-centred (STCC ablation) | 0.59 | — | 0.74 | 0.37 |
| **STCC (this work)** | **0.75** | **0.71** | **0.74** | **0.37** |

The un-centred cosine is STCC without mean-centring, on the same rendered/blurred images; it is an
ablation, not an independent method; its dense-protein margin is not tabulated here (the code
computes a margin for every method, but the stored dense ablation entry carries the separation
only). The dense separations and margins in this table are `bench.py` output, regenerable once the
dense spectra are deposited (Data Availability Statement). The $^1$H–$^{13}$C benchmark uses 5 same-compound pairs and 40
different-compound pairs; a nominal sixth pair (olivetol) was excluded as a byte-identical
duplicate. On the current stick-input sparse benchmark, compound-cluster 95% CIs are 0.508–0.982
for STCC and 0.623–0.985 for experimental Local Contrast. Their paired separation difference is
$\Delta=0.0757$ ($[0.0014,0.2057]$, $P(\Delta\leq0)=0.0067$). Both nevertheless have
AUROC/AUPRC/top-1/MRR = 1.00 and zero leave-one-compound-out error and rejection FPR. The candidate
therefore improves descriptive separation without improving the held-out operating point.

**Table 2.** Expanded dense $^1$H–$^{15}$N benchmark addressing the single-decoy limitation: a PRL3
apo–olsalazine titration (23 same-protein points, exps 4–48) scored against **two** decoy proteins
(OAA, EphB3) in the $^1$H 6.5–10 $\times$ $^{15}$N 105–130 ppm window. Separation = mean same − mean
decoy; margin = worst same − worst decoy. Higher is better; every method self-scores 1.00.
Computed by `python3 bench_nhsqc.py` (harness and derived scores `results/nhsqc_dense.json` are in
the repository; raw spectra per the Data Availability Statement).

| Method | mean same | min same | OAA | EphB3 |
| --- | :---: | :---: | :---: | :---: |
| Local Contrast (experimental) | 0.75 | 0.66 | 0.02 | 0.14 |
| **STCC (default)** | **0.87** | **0.82** | 0.19 | 0.41 |
| Cosine, un-centred (STCC ablation) | 0.90 | 0.85 | 0.36 | 0.51 |
| Bin (Bodis, 2009) | 0.82 | 0.75 | 0.46 | 0.50 |
| Bin + 45° rotation | 0.86 | 0.81 | 0.52 | 0.56 |
| Quad-tree (Castillo, 2013) | 0.95 | 0.91 | 0.79 | 0.75 |
| Nearest neighbour (Pierens, 2012) | 0.99 | 0.98 | 0.96 | 0.96 |

| Method | separation | margin |
| --- | :---: | :---: |
| Local Contrast (experimental) | **0.67** | **0.52** |
| **STCC (default)** | 0.57 | 0.41 |
| Cosine, un-centred (STCC ablation) | 0.46 | 0.35 |
| Bin (Bodis, 2009) | 0.34 | 0.25 |
| Bin + 45° rotation | 0.32 | 0.25 |
| Quad-tree (Castillo, 2013) | 0.18 | 0.12 |
| Nearest neighbour (Pierens, 2012) | 0.03 | 0.02 |

Both decoy proteins score below every one of the 23 titration points for STCC (worst same 0.82 >
worst decoy 0.41) and Local Contrast (0.66 > 0.14), so each margin is positive across two distinct
folds. The candidate suppresses both decoys more strongly, raising separation from 0.57 to 0.67 and
margin from 0.41 to 0.52. It is nevertheless reported as an experimental, post hoc extension: the
held-out sparse operating metrics below do not improve over default STCC. The STCC numbers are lower
than the single-decoy benchmark (Table 1) because this test is harder on both axes: the olsalazine
titration reaches larger shifts (min same 0.82 vs 0.89) and EphB3 is a tougher decoy than OAA
(0.41 vs 0.19).

## Supporting Information

The Supporting Information (full mathematical derivation with proofs, both benchmark protocols,
all per-comparison and per-pair numerical tables, robustness sweeps, parameter reference and
test-suite description) is available and referenced throughout the text.

## Acknowledgements

The author thanks the developers of the simpleNMR project for the public HSQC peak-list example
set used in the sparse-regime benchmark.

## Conflict of Interest

The author declares no conflict of interest.

## Data Availability Statement

The software (STCC, experimental Local Contrast, the three reference methods, the
un-centred-cosine ablation, all benchmark harnesses, the tabulated scores and the 46-test
suite) is freely available under the MIT
license at https://github.com/deepnmr/spectra_similarity_tool. The **sparse $^1$H–$^{13}$C
benchmark and the dereplication/retrieval experiment are fully reproducible from public data**:
the harnesses download the simpleNMR peak lists (https://github.com/EricHughesABC/simpleNMR)
automatically and contain no hard-coded paths. The dense $^1$H–$^{15}$N benchmark uses PRL3 and OAA
Bruker spectra available from the author on reasonable request (the harness reads their location
from `--prl3/--oaa` or `$PRL3_DIR/$OAA_DIR`); a public deposition of the processed spectra is in
preparation so that the dense benchmark also becomes independently reproducible. The expanded
dense benchmark of Table 2 (PRL3–olsalazine titration and the OAA and EphB3 decoy spectra) is
computed by `bench_nhsqc.py` from the same processed Bruker data and is included in that deposition.

## References

[1] S. L. Robinette, R. Brüschweiler, F. C. Schroeder, A. S. Edison, *Acc. Chem. Res.* **2012**,
*45*, 288–297.

[2] K. Bingol, R. Brüschweiler, *Anal. Chem.* **2014**, *86*, 47–57.

[3] B. Worley, R. Powers, *Curr. Metabolomics* **2013**, *1*, 92–107.

[4] E. L. Ulrich, H. Akutsu, J. F. Doreleijers, Y. Harano, Y. E. Ioannidis, J. Lin, M. Livny,
S. Mading, D. Maziuk, Z. Miller, E. Nakatani, C. F. Schulte, D. E. Tolmie, R. K. Wenger, H. Yao,
J. L. Markley, *Nucleic Acids Res.* **2008**, *36*, D402–D408.

[5] C. Steinbeck, S. Kuhn, *Phytochemistry* **2004**, *65*, 2711–2717.

[6] D. S. Wishart, A. Guo, E. Oler, F. Wang, A. Anjum, H. Peters, R. Dizon, Z. Sayeeda, S. Tian,
B. L. Lee, M. Berjanskii, R. Mah, M. Yamamoto, J. Jovel, C. Torres-Calzada, M. Hiebert-Giesbrecht,
V. W. Lui, D. Varshavi, D. Varshavi, D. Allen, D. Arndt, N. Khetarpal, A. Sivakumaran, K. Harford,
S. Sanford, K. Yee, X. Cao, Z. Budinski, J. Liigand, L. Zhang, J. Zheng, R. Mandal, N. Karu,
M. Dambrova, H. B. Schiöth, R. Greiner, V. Gautam, *Nucleic Acids Res.* **2022**, *50*, D622–D631.

[7] S. B. Shuker, P. J. Hajduk, R. P. Meadows, S. W. Fesik, *Science* **1996**, *274*, 1531–1534.

[8] M. Pellecchia, I. Bertini, D. Cowburn, C. Dalvit, E. Giralt, W. Jahnke, T. L. James,
S. W. Homans, H. Kessler, C. Luchinat, B. Meyer, H. Oschkinat, J. Peng, H. Schwalbe, G. Siegal,
*Nat. Rev. Drug Discov.* **2008**, *7*, 738–745.

[9] D. A. Erlanson, S. W. Fesik, R. E. Hubbard, W. Jahnke, H. Jhoti, *Nat. Rev. Drug Discov.*
**2016**, *15*, 605–619.

[10] M. P. Williamson, *Prog. Nucl. Magn. Reson. Spectrosc.* **2013**, *73*, 1–16.

[11] S. A. Berkowitz, J. R. Engen, J. R. Mazzeo, G. B. Jones, *Nat. Rev. Drug Discov.* **2012**,
*11*, 527–540.

[12] J. P. Gabrielson, W. F. Weiss IV, *J. Pharm. Sci.* **2015**, *104*, 1240–1245.

[13] L. W. Arbogast, R. G. Brinson, J. P. Marino, *Anal. Chem.* **2015**, *87*, 3556–3561.

[14] H. Ghasriani, D. J. Hodgson, R. G. Brinson, I. McEwen, L. F. Buhse, S. Kozlowski, J. P. Marino,
Y. Aubin, D. A. Keire, *Nat. Biotechnol.* **2016**, *34*, 139–141.

[15] R. G. Brinson, J. P. Marino, F. Delaglio, L. W. Arbogast, R. M. Evans, A. Kearsley, G. Gingras,
H. Ghasriani, Y. Aubin, G. K. Pierens, X. Jia, M. Mobli, H. G. Grant, D. W. Keizer, K. Schweimer,
J. Ståhle, G. Widmalm, E. R. Zartler, C. W. Lawrence, P. N. Reardon, J. R. Cort, P. Xu, F. Ni,
S. Yanaka, K. Kato, S. R. Parnham, D. Tsao, A. Blomgren, T. Rundlöf, N. Trieloff, P. Schmieder,
A. Ross, K. Skidmore, K. Chen, D. Keire, D. I. Freedberg, T. Suter-Stahel, G. Wider, G. Ilc,
J. Plavec, S. A. Bradley, D. M. Baldisseri, M. L. Sforça, A. C. de M. Zeri, J. Y. Wei, C. M. Szabo,
C. A. Amezcua, J. B. Jordan, M. Wikström, *mAbs* **2019**, *11*, 94–105.

[16] L. W. Arbogast, F. Delaglio, J. E. Schiel, J. P. Marino, *Anal. Chem.* **2017**, *89*,
11839–11845.

[17] D. J. Schuirmann, *J. Pharmacokinet. Biopharm.* **1987**, *15*, 657–680.

[18] S. Wellek, *Testing Statistical Hypotheses of Equivalence and Noninferiority*, 2nd ed.,
Chapman & Hall/CRC, Boca Raton, **2010**.

[19] *ICH Harmonised Tripartite Guideline Q5E: Comparability of Biotechnological/Biological Products
Subject to Changes in Their Manufacturing Process*, International Council for Harmonisation of
Technical Requirements for Pharmaceuticals for Human Use, **2004**.

[20] J. Woodcock, J. Griffin, R. Behrman, B. Cherney, T. Crescenzi, B. Fraser, D. Hixon,
C. Joneckis, S. Kozlowski, A. Rosenberg, L. Schrager, E. Shacter, R. Temple, K. Webber, H. Winkle,
*Nat. Rev. Drug Discov.* **2007**, *6*, 437–442.

[21] L. Bodis, A. Ross, E. Pretsch, *Chemometr. Intell. Lab. Syst.* **2007**, *85*, 1–8.

[22] L. Bodis, A. Ross, J. Bodis, E. Pretsch, *Talanta* **2009**, *79*, 1379–1386.

[23] A. Craig, O. Cloarec, E. Holmes, J. K. Nicholson, J. C. Lindon, *Anal. Chem.* **2006**,
*78*, 2262–2267.

[24] F. Dieterle, A. Ross, G. Schlotterbeck, H. Senn, *Anal. Chem.* **2006**, *78*, 4281–4290.

[25] A. M. Castillo, L. Uribe, L. Patiny, J. Wist, *Chemometr. Intell. Lab. Syst.* **2013**,
*127*, 1–6.

[26] G. K. Pierens, S. Brossi, Z. Yang, D. C. Reutens, V. Vegh, *J. Cheminform.* **2012**, *4*, 25.

[27] F. Delaglio, S. Grzesiek, G. W. Vuister, G. Zhu, J. Pfeifer, A. Bax, *J. Biomol. NMR*
**1995**, *6*, 277–293.

[28] W. F. Vranken, W. Boucher, T. J. Stevens, R. H. Fogh, A. Pajon, M. Llinas, E. L. Ulrich,
J. L. Markley, J. Ionides, E. D. Laue, *Proteins* **2005**, *59*, 687–696.

[29] S. P. Skinner, R. H. Fogh, W. Boucher, T. J. Ragan, L. G. Mureddu, G. W. Vuister,
*J. Biomol. NMR* **2016**, *66*, 111–124.

[30] C. Ludwig, U. L. Günther, *BMC Bioinformatics* **2011**, *12*, 366.

[31] J. J. Helmus, C. P. Jaroniec, *J. Biomol. NMR* **2013**, *55*, 355–367.

[32] S. E. Stein, D. R. Scott, *J. Am. Soc. Mass Spectrom.* **1994**, *5*, 859–866.

[33] K. X. Wan, I. Vidavsky, M. L. Gross, *J. Am. Soc. Mass Spectrom.* **2002**, *13*, 85–88.

[34] J. P. Lewis, *Vision Interface* **1995**, *10*, 120–123.

[35] F. Savorani, G. Tomasi, S. B. Engelsen, *J. Magn. Reson.* **2010**, *202*, 190–202.

[36] C. R. Harris, K. J. Millman, S. J. van der Walt, R. Gommers, P. Virtanen, D. Cournapeau,
E. Wieser, J. Taylor, S. Berg, N. J. Smith, R. Kern, M. Picus, S. Hoyer, M. H. van Kerkwijk,
M. Brett, A. Haldane, J. F. del Río, M. Wiebe, P. Peterson, P. Gérard-Marchant, K. Sheppard,
T. Reddy, W. Weckesser, H. Abbasi, C. Gohlke, T. E. Oliphant, *Nature* **2020**, *585*, 357–362.

[37] J. D. Hunter, *Comput. Sci. Eng.* **2007**, *9*, 90–95.

[38] E. Hughes, A. M. Kenwright, *Magn. Reson. Chem.* **2024**, *62*, 556–565.

[39] C. E. Hart, T. Kind, P. C. Dorrestein, D. Healey, D. Domingo-Fernández,
*J. Am. Soc. Mass Spectrom.* **2024**, *35*, 266–274.

---

## Table of Contents Entry

**Tolerant yet discriminating.** Testing whether two 2D-HSQC spectra are the same demands tolerance
of small chemical-shift drifts *and* penalization of genuinely relocated intensity — requirements the
incumbent binning, tree and nearest-neighbour measures satisfy only one at a time. The Shift-Tolerant
Correlation Coefficient renders both spectra onto one grid, blurs them by the physical NMR
linewidth, and takes the mean-centred zero-lag correlation. An experimental local-contrast
extension further improves descriptive separation across dense protein and sparse small-molecule
examples, while the validated STCC remains the default.

*Suggested TOC graphic: the separation bar chart (Figure 2) with the three-step STCC schematic
inset.*
