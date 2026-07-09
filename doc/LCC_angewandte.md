---
title: "A Lineshape Correlation Coefficient for Robust Similarity of Two-Dimensional HSQC NMR Spectra"
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
---

**NMR Spectroscopy | Research Article**

Korea Basic Science Institute (KBSI), Ochang, Republic of Korea. \* Correspondence:
kbsi.bionmr@gmail.com

---

**Abstract.** Quantitative comparison of two-dimensional $^1$H–X HSQC spectra underlies
ligand screening, chemical-shift-perturbation mapping, higher-order-structure comparability of
protein biotherapeutics, and molecular fingerprinting, yet
the established similarity measures each fail in one of two opposite ways: histogram-binning
scores are brittle to the small peak drifts a titration induces, whereas tree- and
nearest-neighbour scores saturate because in any crowded window every peak finds a
coincidental near neighbour. Here we introduce the Lineshape Correlation Coefficient (LCC),
which renders both spectra onto a single shared grid, blurs each by the physical NMR
linewidth combined with the expected shift drift, and scores the pair by the mean-centred,
zero-lag normalized cross-correlation of the two images. The Gaussian render supplies graded,
physically-tunable shift tolerance while mean-centring actively penalizes intensity that is
not co-located. On two deliberately opposite benchmarks — a dense protein $^1$H–$^{15}$N amide
fingerprint and a sparse small-molecule $^1$H–$^{13}$C set — LCC gives the widest same/different
separation in both regimes (0.75 dense, 0.78 sparse on a $[0,1]$ scale), significantly outscoring
the saturating tree and nearest-neighbour measures, which leave the classes overlapping (margin
$\le 0.09$) in both. A bootstrap shows LCC's sparse-regime edge over the bin method is within noise
at five compounds, while its dense-fingerprint advantage is large and, on an expanded test, keeps
both of two decoy proteins below every same-protein titration point; removing the mean-centring step
collapses the dense-protein separation from 0.75 to 0.59, yet in the sparse regime — where
peaks rarely overlap — centred and un-centred scores coincide (0.78), so mean-centring is the
operation that rescues discrimination specifically when peaks crowd. LCC is a two-knob,
peak-picking-free, solver-free measure with both knobs set by spectroscopy rather than by tuning.

**Keywords:** NMR spectroscopy · HSQC · spectral similarity · higher-order structure · chemical-shift perturbation

---

The comparison of two two-dimensional $^1$H–X HSQC spectra is a deceptively basic operation
that sits at the centre of several NMR workflows. In fragment-based and structure-based drug
discovery, an amide $^1$H–$^{15}$N fingerprint is tracked across a ligand titration to read out
binding through chemical-shift perturbation,$^{[1-4]}$ and the quality of that read-out depends
on being able to say, quantitatively, how different two fingerprints are. In metabolomics and
natural-product chemistry, a small-molecule $^1$H–$^{13}$C spectrum is matched against a
reference — an operation central to 2D-NMR-based dereplication$^{[5-7]}$ and to searches against
community libraries such as the BMRB, NMRShiftDB and HMDB.$^{[8-10]}$

A third and increasingly prominent setting is the **higher-order-structure (HOS)** assessment of
protein biopharmaceuticals. A protein's biological activity, stability and immunogenicity are
governed not by its amino-acid sequence alone but by its folded secondary, tertiary and quaternary
structure, so confirming that this higher-order structure is preserved — between a biosimilar and
its originator, across a manufacturing or formulation change, or on storage — has become a central
requirement of biochemical characterization and of regulatory comparability.$^{[11,12]}$
Two-dimensional HSQC is a premier atomic-resolution probe of that structure: each crosspeak reports
the local environment of one $^1$H–X pair, so two fingerprints coincide only when the two samples
share the same fold, and 2D $^{13}$C- and $^{15}$N-HSQC fingerprinting is now an established,
high-precision method for the HOS comparability of therapeutic proteins and monoclonal
antibodies.$^{[13,14,16]}$ Deciding *how similar* two such fingerprints are is exactly the
measurement that multivariate and chemometric HOS workflows carry out,$^{[15]}$ and exactly the
primitive the present work makes robust.

In every case the same
primitive is required: a single number that is high when two spectra come from the same species
and low when they do not.

Three families of measure are in current use, and they fail in two opposite ways. The bin
method of Bodis et al. integrates each spectrum over an $n \times n$ grid at a range of
resolutions and reports a weighted Jaccard (Ružička) index,$^{[17,18]}$ in the lineage of the
binning and normalization schemes long used in NMR metabonomics.$^{[19,20]}$ Because its bin
edges are hard, a peak straddling an edge splits its intensity across neighbouring bins; a small
titration drift therefore erodes the score in a discontinuous jump. The measure is
shift-**brittle**. The quad-tree method of Castillo et al. encodes a spectrum as a recursive
centre-of-mass tree and compares trees node-by-node with a shift-tolerant node score,$^{[21]}$
and the peak method of Pierens et al. picks peaks and matches each to its nearest neighbour in
the other spectrum.$^{[22]}$ Both build in genuine shift tolerance — but in a crowded window
every peak has *some* near neighbour and every spectrum a similar mass-centre structure, so both
saturate near unity and retain almost no contrast between related and unrelated spectra. They
are shift-**blind**. A measure fit for real HSQC comparison must tolerate the small physical
drifts a titration or a solvent change induces, yet still penalize intensity that is genuinely
relocated.

We resolve this tension with the Lineshape Correlation Coefficient (LCC), a measure that keeps
the physically meaningful ingredients of all three predecessors — area-weighted integration and
unit normalization (Bodis), a physical shift tolerance (Castillo), and a lineshape picture of
the spectrum (Pierens) — and adds a mean-centred correlation as the discriminating step. Given
two spectra compared inside a common window $\Omega$ (processed Bruker data from any standard
pipeline$^{[23-27]}$ or peak lists rendered to a grid), LCC proceeds in three steps (Figure 1a).

**Render.** Negative intensities are clipped, each point is weighted by its local integration
area, and *both* spectra are histogrammed onto a **single shared grid**, so they are directly
comparable irrespective of their native digital resolution. The mass is normalized to unit sum;
the final score is invariant to this scale.

**Blur.** Each image is convolved with a separable Gaussian whose width per axis is the physical
linewidth combined with the expected drift, $\sigma = \sqrt{\ell^2 + d^2}$. This replaces the
bin method's hard edges with a continuous response. For an isolated peak on a large grid, two
copies displaced by $(\Delta_2, \Delta_1)$ correlate as
$\exp(-\Delta_2^2/4\sigma_2^2 - \Delta_1^2/4\sigma_1^2)$ — smooth and monotone in the drift
(Figure 1b), so a physical shift costs little while a random relocation ($\Delta \gg \sigma$)
costs almost everything. The blur width, not a bin count, is the tolerance knob, and it is set
by two spectroscopic quantities: the linewidth (the floor) and the expected shift perturbation
(the tolerance).

**Score.** Each image is mean-centred and the similarity is the zero-lag normalized
cross-correlation — the Pearson coefficient — of the two centred images, clamped to $[0,1]$.
This is the spectral-imaging analogue of the cosine / contrast-angle similarity long used in
mass-spectral library search.$^{[28,29]}$ Mean-centring is what turns overlap into
discrimination: a cell where one spectrum has a peak (a value above its mean) and the other is
empty (a value below its mean) contributes a **negative** product and lowers the score, so
intensity that is not co-located is actively penalized rather than merely ignored, as it is by
the one-way nearest-neighbour distance. No shift search is performed. Aligning the images — as
in NMR peak-alignment tools such as icoshift$^{[30]}$ — would let a genuinely different spectrum
slide into registration and re-saturate, reintroducing exactly the failure mode of the tree and
nearest-neighbour methods. Scoring at zero lag keeps the position information that carries the
discrimination.

The measure has the properties one expects of a similarity: self-similarity is exactly one, the
score is exactly symmetric, and it is monotone in the drift of an isolated peak (all verified in
the test suite; Supporting Information). The implementation is pure NumPy,$^{[31]}$ reads
processed Bruker data directly, requires neither peak picking nor an external solver, and runs
in about 6 ms per spectrum pair on a $250 \times 350$ grid (Apple Silicon, single core). A
26-test suite (all passing) checks self-similarity, related-beats-unrelated ordering, and
monotonicity in drift, plus regressions for the mean-centring ablation, the `abs` baseline and the
duplicate-pair guard; figures use Matplotlib.$^{[32]}$ A full derivation with proofs is given in
the Supporting Information.

We benchmarked LCC against the three reference methods on two deliberately **opposite** regimes,
so that any measure claiming general utility must succeed on both (Figure 2; Table 1). We
summarize discrimination by the **separation** (mean same-class score − mean different-class
score) and the **margin** (worst same-class − best different-class score, which is what a single
classification threshold actually sees); both are reported on each method's own $[0,1]$ scale.
Every method self-scores exactly 1.00.

*Dense protein $^1$H–$^{15}$N.* The reference is the apo spectrum of the phosphatase PRL3,
compared with six spectra of the same protein along a ligand titration (which should score high)
and with a **different** protein (which should score low), inside a
$^1$H 6.5–10 $\times$ $^{15}$N 105–130 ppm window. LCC scores the same protein 0.94 (range
0.89–0.97) and the different protein 0.18, a margin of 0.71 — three times the bin method's 0.24,
the strongest baseline. Mean-centring earns that gap: replacing it with an un-centred cosine
overlap on the *same* images drops the separation from 0.75 to 0.59. This gain is not a tuned
artefact; it holds across the entire physical blur range and even survives coarsening the grid
to the bin method's own resolution (Supporting Information).

*A second decoy protein and an independent titration.* Because that margin is measured against a
single decoy protein, we tested whether LCC's advantage survives a harder, self-contained design:
an independent PRL3 apo→olsalazine titration (23 points, from a 500 µM sample) scored against
**two** distinct decoy folds — OAA and the kinase-domain EphB3 — inside the same window (Table 2;
Supporting Information). As ligand is added the chemical-shift perturbation grows and LCC decays
**smoothly** from 0.98 to a plateau near 0.82 — its graded shift tolerance acting on a real
titration — while both decoys sit far below every titration point (OAA 0.19, EphB3 0.41). LCC keeps
the widest separation (0.57) and the widest margin (0.41): both decoy proteins clear *every*
same-protein score by a comfortable gap, whereas the tree and nearest-neighbour margins are a
knife-edge (0.12 and 0.02, a single noisy titration point from crossing). The un-centred cosine
trails (0.46 / 0.35), so mean-centring still earns the gap against a genuine negative distribution. This test is harder than the single-decoy benchmark on both axes —
a longer titration with larger shifts and a second, tougher decoy — yet the ordering is unchanged,
so LCC's dense-regime margin is not an artefact of one negative spectrum.

*Sparse small-molecule $^1$H–$^{13}$C.* To probe the opposite regime — and to replace the single
negative spectrum above with a large negative set — we used public HSQC peak lists from the
simpleNMR example set$^{[33]}$ for five compounds each recorded twice (one nominal sixth pair was
excluded on inspection as a byte-identical duplicate; Supporting Information), and asked each
method to score the 5 same-compound pairs above the 40 different-compound pairs, inside a
$^1$H 0–10 $\times$ $^{13}$C 0–165 ppm window. LCC again leads (mean same 0.78, mean different
0.003, separation 0.78). The bin method is far stronger here (0.65) than on the dense
fingerprint (0.29): with few, well-separated peaks its hard bins rarely straddle a peak. The
tree and nearest-neighbour methods, however, **leave the two classes barely separable in this
regime too** (margins 0.09 and 0.04, despite separations of 0.41 and 0.10): over a wide window a
different molecule still supplies a near neighbour for every peak. Their shift tolerance is
real — per pair, it keeps two solvent-shifted recordings of the same compound high (0.83–1.00,
where the position-sensitive LCC and bin methods drop to 0.30–0.60) — but it tolerates
*everything*, so it never becomes discrimination. LCC recovers that same tolerance
controllably, through its blur width, without surrendering the between-compound contrast
(Supporting Information). Notably, in this sparse regime the un-centred cosine matches LCC exactly
(separation 0.78): where crosspeaks rarely overlap there is no coincidental non-co-located
intensity for mean-centring to penalize, so its benefit is specific to the dense fingerprint.

*Statistical strength and a retrieval test.* Two additions probe how much of this the small
samples actually support (Figure 3). First, bootstrapping the sparse-regime separation (resampling
the 5 same- and 40 different-compound pairs, $10^4$ draws) gives LCC 0.78 (95% CI 0.55–0.98),
**overlapping** the bin method's 0.65 (0.42–0.84) but cleanly above the saturating tree (0.41,
0.34–0.49) and nearest-neighbour (0.10, 0.09–0.11) measures (Figure 3a): on this five-compound set
LCC's edge over the bin method is within noise, while its advantage over the two shift-tolerant
methods is significant. Mean-centring is what separates LCC from the un-centred cosine, and only
where peaks crowd: the ablation gains $+0.16$ on the dense fingerprint but nothing on the sparse set
(Figure 3b). Second, cast as a library search — rank the other nine spectra for each of the ten and
ask whether the true same-compound partner ranks first — **every** method scores top-1 = 1.00 (mean
reciprocal rank 1.00): on a small, clean library the correct partner is always the single closest,
even for the saturating methods, so retrieval becomes discriminating only under a rejection
threshold or at library scale — which is exactly what the separation and margin measure. The
dense-protein margin (0.71) remains LCC's strongest single number, but rather than rest it on one
decoy we add an expanded dense test with a second decoy protein (EphB3) and an independent PRL3
titration (Table 2), on which LCC keeps *both* decoys below *every* same-protein point; the sparse
regime and this retrieval test further broaden the evidence with a larger, public,
fully-reproducible negative set.

Across both regimes LCC is thus the only measure that beats the saturating tree and
nearest-neighbour methods while retaining their shift tolerance, and it matches or exceeds the bin
method — decisively on the dense amide fingerprint (margin 0.71 vs 0.24), comparably on the sparse
set (within the bootstrap CI) (Table 1). Its two knobs, the per-axis blur widths, are both
physical: the linewidth sets the floor and the expected shift perturbation sets how much drift is
tolerated. Two caveats bound the claim honestly. First, LCC
measures position-and-intensity coincidence, not molecular identity, so a close structural
analogue will score intermediate — the correct behaviour for a *similarity*, but a caveat for
hard classification. Second, the three reference methods are our own representative reimplementations
compared inside one common window, and the $^1$H–$^{13}$C spectra are rendered from peak lists;
accordingly the robust **ranking**, not the exact numerical values, is the result, and this is
not a verdict on the original methods as their authors deployed them with per-spectrum windows.

In summary, a single, physically-motivated construction — render onto a shared grid, blur by the
NMR linewidth, and take the mean-centred zero-lag correlation — reconciles shift tolerance with
discrimination, the two requirements that the incumbent measures satisfy only one at a time. Both
of LCC's knobs are spectroscopic quantities rather than fitting parameters (its self-similarity is
exact, as for every method tested), and it needs no peak picking, no alignment search, and no
external solver. We are candid about scope: this is a **methodological** contribution — a robust,
reproducible similarity primitive, not a new chemical discovery — and although its most decisive
number, the dense-protein margin, was first measured against one decoy, we now show it holds against
two decoy proteins on an independent titration (Table 2) and add the larger public sparse benchmark
and the retrieval test. Within those bounds LCC is directly useful wherever a
robust 2D-HSQC similarity is the bottleneck: automated titration tracking and CSP read-out, library
and dereplication searches, quality control of repeated measurements, and the higher-order-structure
comparability assessment of protein biotherapeutics and biosimilars, where 2D-NMR fingerprint
comparison is already an established readout.$^{[15,16]}$ The natural next steps are a prospective,
at-scale dereplication study against a community library and a formulation/biosimilar HOS
comparison, for both of which LCC supplies the scoring primitive. The full derivation, all three benchmark protocols, every numerical table with
confidence intervals, and the 26-test suite are provided as Supporting Information; the software is
freely available under the MIT license.

## Figures

![**Figure 1.** The Lineshape Correlation Coefficient. (a) The three-step construction:
area-weighted rendering of both spectra onto one shared grid, a separable Gaussian blur of
physical width $\sigma = \sqrt{\ell^2 + d^2}$ per axis, and the mean-centred zero-lag normalized
cross-correlation of the two images (red/blue = positive/negative contribution to the score).
(b) Graded shift tolerance: two copies of an isolated peak displaced by $\Delta$ correlate as
$\exp(-\Delta^2/4\sigma^2)$, smooth and monotone in the drift, so a physical titration drift
(shaded) costs little while a random relocation costs almost
everything.](../results/lcc_schematic.png)

![**Figure 2.** Same/different separation (mean same-class − mean different-class similarity) for
every method on the two benchmarks: dense protein $^1$H–$^{15}$N (left of each pair) and sparse
small-molecule $^1$H–$^{13}$C (right). LCC leads in both regimes; the tree and nearest-neighbour
methods saturate in both. The un-centred cosine (an ablation of LCC without mean-centring) matches
LCC on the sparse regime but falls from 0.75 to 0.59 on the dense fingerprint, isolating
mean-centring as the step that rescues discrimination when peaks
crowd.](../results/comparison_all.png)

![**Figure 3.** Statistical strength and the mean-centring ablation. **(a)** Sparse
$^1$H–$^{13}$C separation with a bootstrap 95% CI ($10^4$ resamples of the 5 same- and 40
different-compound pairs): LCC (blue) and the bin method (orange) have **overlapping** intervals —
their difference is within noise at five compounds — while both sit cleanly above the saturating
quad-tree and nearest-neighbour measures (grey), and the un-centred cosine (light blue, the LCC
ablation) coincides with LCC. **(b)** Mean-centring ablation by regime: replacing the mean-centred
correlation with the un-centred cosine costs 0.16 of separation on the dense protein fingerprint
(0.75 → 0.59) but nothing on the sparse set (0.78 → 0.78) — mean-centring rescues discrimination
only where peaks crowd.](../results/fig3_stats_ablation.png)

## Table 1

**Table 1.** Separation and worst-case margin for all methods on both regimes, computed from
unrounded scores. Higher is better; every method self-scores 1.00.

| Method | $^1$H–$^{15}$N sep | margin | $^1$H–$^{13}$C sep | margin |
| --- | :---: | :---: | :---: | :---: |
| Bin (Bodis, 2009) | 0.29 | 0.24 | 0.65 | 0.20 |
| Bin + 45° rotation | 0.25 | 0.20 | 0.65 | 0.18 |
| Quad-tree (Castillo, 2013) | 0.03 | −0.01 | 0.41 | 0.09 |
| Nearest neighbour (Pierens, 2012) | 0.04 | 0.03 | 0.10 | 0.04 |
| Cosine, un-centred (LCC ablation) | 0.59 | — | 0.78 | 0.37 |
| **LCC (this work)** | **0.75** | **0.71** | **0.78** | **0.37** |

The un-centred cosine is LCC without mean-centring, on the same rendered/blurred images; it is an
ablation, not an independent method. Its dense-protein margin was not recomputed here (the
ablation reports separation only). The $^1$H–$^{13}$C benchmark uses 5 same-compound pairs and 40
different-compound pairs; a nominal sixth pair (olivetol) was excluded as a byte-identical
duplicate. Bootstrap 95% CIs on the sparse separation ($10^4$ resamples): LCC 0.55–0.98, bin
0.42–0.84 (overlapping LCC), tree 0.34–0.49, NN 0.09–0.11 (both below LCC). Cast as a library
search over the 10 spectra, all six methods retrieve the correct same-compound partner at rank 1
(top-1 = 1.00, MRR = 1.00) — retrieval on a small clean library does not separate the methods; the
thresholded separation/margin does.

## Table 2

**Table 2.** Expanded dense $^1$H–$^{15}$N benchmark addressing the single-decoy limitation: a PRL3
apo–olsalazine titration (23 same-protein points, exps 4–48) scored against **two** decoy proteins
(OAA, EphB3) in the $^1$H 6.5–10 $\times$ $^{15}$N 105–130 ppm window. Separation = mean same − mean
decoy; margin = worst same − worst decoy. Higher is better; every method self-scores 1.00.
Reproducible with `python3 bench_nhsqc.py`.

| Method | mean same | min same | OAA | EphB3 | separation | margin |
| --- | :---: | :---: | :---: | :---: | :---: | :---: |
| **LCC (this work)** | **0.87** | **0.82** | 0.19 | 0.41 | **0.57** | **0.41** |
| Cosine, un-centred (LCC ablation) | 0.90 | 0.85 | 0.36 | 0.51 | 0.46 | 0.35 |
| Bin (Bodis, 2009) | 0.82 | 0.75 | 0.46 | 0.50 | 0.34 | 0.25 |
| Bin + 45° rotation | 0.86 | 0.81 | 0.52 | 0.56 | 0.32 | 0.25 |
| Quad-tree (Castillo, 2013) | 0.95 | 0.91 | 0.79 | 0.75 | 0.18 | 0.12 |
| Nearest neighbour (Pierens, 2012) | 0.99 | 0.98 | 0.96 | 0.96 | 0.03 | 0.02 |

Both decoy proteins score below every one of the 23 titration points for LCC (worst same 0.82 >
worst decoy 0.41), so the margin is positive across two distinct folds — the negative side is now a
distribution, not a single spectrum. The numbers are lower than the single-decoy benchmark (Table 1)
because this test is harder on both axes: the olsalazine titration reaches larger shifts (min same
0.82 vs 0.89) and EphB3 is a tougher decoy than OAA (0.41 vs 0.19); LCC's ordering is nevertheless
unchanged.

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

The software (LCC, the three reference methods, the un-centred-cosine ablation, all three
benchmark harnesses, the tabulated scores and the 26-test suite) is freely available under the MIT
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

[1] S. B. Shuker, P. J. Hajduk, R. P. Meadows, S. W. Fesik, *Science* **1996**, *274*, 1531–1534.

[2] M. Pellecchia, I. Bertini, D. Cowburn, C. Dalvit, E. Giralt, W. Jahnke, T. L. James,
S. W. Homans, H. Kessler, C. Luchinat, B. Meyer, H. Oschkinat, J. Peng, H. Schwalbe, G. Siegal,
*Nat. Rev. Drug Discov.* **2008**, *7*, 738–745.

[3] D. A. Erlanson, S. W. Fesik, R. E. Hubbard, W. Jahnke, H. Jhoti, *Nat. Rev. Drug Discov.*
**2016**, *15*, 605–619.

[4] M. P. Williamson, *Prog. Nucl. Magn. Reson. Spectrosc.* **2013**, *73*, 1–16.

[5] S. L. Robinette, R. Brüschweiler, F. C. Schroeder, A. S. Edison, *Acc. Chem. Res.* **2012**,
*45*, 288–297.

[6] K. Bingol, R. Brüschweiler, *Anal. Chem.* **2014**, *86*, 47–57.

[7] B. Worley, R. Powers, *Curr. Metabolomics* **2013**, *1*, 92–107.

[8] E. L. Ulrich, H. Akutsu, J. F. Doreleijers, Y. Harano, Y. E. Ioannidis, J. Lin, M. Livny,
S. Mading, D. Maziuk, Z. Miller, E. Nakatani, C. F. Schulte, D. E. Tolmie, R. K. Wenger, H. Yao,
J. L. Markley, *Nucleic Acids Res.* **2008**, *36*, D402–D408.

[9] C. Steinbeck, S. Kuhn, *Phytochemistry* **2004**, *65*, 2711–2717.

[10] D. S. Wishart, A. Guo, E. Oler, F. Wang, A. Anjum, H. Peters, R. Dizon, Z. Sayeeda, S. Tian,
B. L. Lee, M. Berjanskii, R. Mah, M. Yamamoto, J. Jovel, C. Torres-Calzada, M. Hiebert-Giesbrecht,
V. W. Lui, D. Varshavi, D. Varshavi, D. Allen, D. Arndt, N. Khetarpal, A. Sivakumaran, K. Harford,
S. Sanford, K. Yee, X. Cao, Z. Budinski, J. Liigand, L. Zhang, J. Zheng, R. Mandal, N. Karu,
M. Dambrova, H. B. Schiöth, R. Greiner, V. Gautam, *Nucleic Acids Res.* **2022**, *50*, D622–D631.

[11] S. A. Berkowitz, J. R. Engen, J. R. Mazzeo, G. B. Jones, *Nat. Rev. Drug Discov.* **2012**,
*11*, 527–540.

[12] J. P. Gabrielson, W. F. Weiss IV, *J. Pharm. Sci.* **2015**, *104*, 1240–1245.

[13] L. W. Arbogast, R. G. Brinson, J. P. Marino, *Anal. Chem.* **2015**, *87*, 3556–3561.

[14] H. Ghasriani, D. J. Hodgson, R. G. Brinson, I. McEwen, L. F. Buhse, S. Kozlowski, J. P. Marino,
Y. Aubin, D. A. Keire, *Nat. Biotechnol.* **2016**, *34*, 139–141.

[15] L. W. Arbogast, F. Delaglio, J. E. Schiel, J. P. Marino, *Anal. Chem.* **2017**, *89*,
11839–11845.

[16] R. G. Brinson, J. P. Marino, F. Delaglio, L. W. Arbogast, et al., *mAbs* **2019**, *11*, 94–105.

[17] L. Bodis, A. Ross, E. Pretsch, *Chemometr. Intell. Lab. Syst.* **2007**, *85*, 1–8.

[18] L. Bodis, A. Ross, J. Bodis, E. Pretsch, *Talanta* **2009**, *79*, 1379–1386.

[19] A. Craig, O. Cloarec, E. Holmes, J. K. Nicholson, J. C. Lindon, *Anal. Chem.* **2006**,
*78*, 2262–2267.

[20] F. Dieterle, A. Ross, G. Schlotterbeck, H. Senn, *Anal. Chem.* **2006**, *78*, 4281–4290.

[21] A. M. Castillo, L. Uribe, L. Patiny, J. Wist, *Chemometr. Intell. Lab. Syst.* **2013**,
*127*, 1–6.

[22] G. K. Pierens, S. Brossi, Z. Yang, D. C. Reutens, V. Vegh, *J. Cheminform.* **2012**, *4*, 25.

[23] F. Delaglio, S. Grzesiek, G. W. Vuister, G. Zhu, J. Pfeifer, A. Bax, *J. Biomol. NMR*
**1995**, *6*, 277–293.

[24] W. F. Vranken, W. Boucher, T. J. Stevens, R. H. Fogh, A. Pajon, M. Llinas, E. L. Ulrich,
J. L. Markley, J. Ionides, E. D. Laue, *Proteins* **2005**, *59*, 687–696.

[25] S. P. Skinner, R. H. Fogh, W. Boucher, T. J. Ragan, L. G. Mureddu, G. W. Vuister,
*J. Biomol. NMR* **2016**, *66*, 111–124.

[26] C. Ludwig, U. L. Günther, *BMC Bioinformatics* **2011**, *12*, 366.

[27] J. J. Helmus, C. P. Jaroniec, *J. Biomol. NMR* **2013**, *55*, 355–367.

[28] S. E. Stein, D. R. Scott, *J. Am. Soc. Mass Spectrom.* **1994**, *5*, 859–866.

[29] K. X. Wan, I. Vidavsky, M. L. Gross, *J. Am. Soc. Mass Spectrom.* **2002**, *13*, 85–88.

[30] F. Savorani, G. Tomasi, S. B. Engelsen, *J. Magn. Reson.* **2010**, *202*, 190–202.

[31] C. R. Harris, K. J. Millman, S. J. van der Walt, R. Gommers, P. Virtanen, D. Cournapeau,
E. Wieser, J. Taylor, S. Berg, N. J. Smith, R. Kern, M. Picus, S. Hoyer, M. H. van Kerkwijk,
M. Brett, A. Haldane, J. F. del Río, M. Wiebe, P. Peterson, P. Gérard-Marchant, K. Sheppard,
T. Reddy, W. Weckesser, H. Abbasi, C. Gohlke, T. E. Oliphant, *Nature* **2020**, *585*, 357–362.

[32] J. D. Hunter, *Comput. Sci. Eng.* **2007**, *9*, 90–95.

[33] E. Hughes, A. M. Kenwright, *Magn. Reson. Chem.* **2024**, *62*, 556–565.

---

## Table of Contents Entry

**Tolerant yet discriminating.** Comparing two 2D-HSQC spectra demands tolerance of small
chemical-shift drifts *and* penalization of genuinely relocated intensity — requirements the
incumbent binning, tree and nearest-neighbour measures satisfy only one at a time. The Lineshape
Correlation Coefficient renders both spectra onto one grid, blurs them by the physical NMR
linewidth, and takes the mean-centred zero-lag correlation, giving the widest same/different
separation across both dense protein and sparse small-molecule regimes.

*Suggested TOC graphic: the separation bar chart (Figure 2) with the three-step LCC schematic
inset.*
