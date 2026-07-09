---
title: "LCC: a lineshape-correlation coefficient for 2D HSQC spectral similarity"
author: "Donghan Lee"
date: " "
geometry: margin=1in
colorlinks: true
---

Korea Basic Science Institute (KBSI), Ochang, Republic of Korea. Correspondence:
kbsi.bionmr@gmail.com

**Subject area:** Structural bioinformatics / NMR spectroscopy

## Abstract

**Motivation:** Quantitative comparison of two-dimensional HSQC spectra underlies ligand
screening and molecular fingerprinting, yet the established similarity measures each fail in
one of two opposite ways. On a dense protein $^1$H-$^{15}$N amide fingerprint, bin-histogram
methods lose margin under the small peak drifts a titration induces, while tree and
nearest-neighbour methods saturate because in a crowded fingerprint every peak finds a
coincidental near neighbour. A measure is needed that tolerates small physical drifts yet
still penalizes genuinely relocated intensity.

**Results:** We introduce the Lineshape Correlation Coefficient (LCC), which renders each
spectrum onto a shared grid, blurs it by the physical linewidth plus the expected shift drift,
and scores the two images by their mean-centred normalized cross-correlation at zero lag. The
Gaussian render gives graded shift tolerance while mean-centring penalizes intensity that is
not co-located. On two independent benchmarks — a dense protein $^1$H-$^{15}$N titration
(same protein vs a different protein) and a sparse small-molecule $^1$H-$^{13}$C set (six
compounds each recorded twice, 6 same-compound vs 60 different-compound pairs) — LCC gives the
widest same/different separation of all methods tested (0.75 and 0.81); the tree and
nearest-neighbour methods leave the two classes overlapping (margin $\le 0.04$) in both regimes, and
the bin method trails (separation 0.29 and 0.69). Removing the mean-centring drops the protein
separation from 0.75 to 0.59,
confirming it as the discriminating step.

**Availability and implementation:** Freely available under the MIT license at
https://github.com/deepnmr/spectra_similarity_tool; implemented in Python 3 (NumPy; Matplotlib
optional for figures); platform independent. The repository includes the three reference
methods, both benchmark harnesses and a test suite.

**Contact:** kbsi.bionmr@gmail.com

## 1 Introduction

A recurring task in NMR is to score how similar two 2D HSQC spectra are: tracking a
$^1$H-$^{15}$N amide fingerprint across a ligand titration, mapping a chemical shift
perturbation on binding (Shuker *et al.*, 1996; Pellecchia *et al.*, 2008; Erlanson *et al.*,
2016; Williamson, 2013), or matching a small-molecule $^1$H-$^{13}$C spectrum against a
database — a comparison central to 2D-NMR metabolomics and natural-product work (Robinette
*et al.*, 2012; Bingol and Brüschweiler, 2014; Worley and Powers, 2013), where spectra are also
matched against community libraries (Ulrich *et al.*, 2008; Steinbeck and Kuhn, 2004; Wishart
*et al.*, 2022). Three families of measure are in use. The bin method of Bodis *et al.*
integrates each spectrum over an $n \times n$ grid at a range of resolutions and reports a
weighted Jaccard (Ružička) index (Bodis *et al.*, 2007, 2009), in the lineage of the binning and
normalization schemes used across NMR metabonomics (Craig *et al.*, 2006; Dieterle *et al.*,
2006). The tree method of Castillo *et al.* encodes a spectrum as a
recursive centre-of-mass quad-tree and compares the trees node-by-node with a shift-tolerant
node score (Castillo *et al.*, 2013). The peak method of Pierens *et al.* picks peaks and
matches each to its nearest neighbour in the other spectrum (Pierens *et al.*, 2012).

These fail in two opposite ways. Bin histograms are shift-**brittle**: a peak straddling a bin
edge splits its intensity across neighbouring bins, so a small drift erodes the score. The tree
and nearest-neighbour measures are shift-**blind**: when many peaks share a region every peak
has some near neighbour and every spectrum a similar mass-centre structure, so both saturate and
retain little contrast. A measure for real HSQC comparison must tolerate small physical drifts
yet still penalize genuinely relocated intensity.

## 2 Implementation

LCC keeps the physically meaningful ingredients of the three methods — area-weighted integration
and unit normalization (Bodis), a physical shift tolerance (Castillo), a lineshape picture of the
spectrum (Pierens) — and adds a mean-centred correlation as the discriminating step. Given two
spectra (processed Bruker data from any standard pipeline — Delaglio *et al.*, 1995; Vranken *et
al.*, 2005; Skinner *et al.*, 2016; Ludwig and Günther, 2011 — or peak lists) compared inside a
common window, it proceeds in three steps.

**Render.** Negative intensities are clipped, each point is weighted by its local integration
area, and both spectra are histogrammed onto a **single shared grid** so they are directly
comparable irrespective of their native digital resolution.

**Blur.** Each image is convolved with a separable Gaussian whose width per axis is the physical
linewidth combined with the expected drift, $\sigma = \sqrt{\ell^2 + d^2}$. This replaces the bin
method's hard edges with a continuous response: for an isolated peak on a large grid (where
mean-centring is negligible) two copies displaced by $(\Delta_1, \Delta_2)$ correlate as
$\exp(-\Delta_1^2/4\sigma_1^2 - \Delta_2^2/4\sigma_2^2)$, smooth and monotone in the drift, so a
physical shift costs little and a random relocation costs a lot.

**Score.** Each image is mean-centred and the similarity is the zero-lag normalized
cross-correlation (the Pearson coefficient) of the two centred images, clamped to $[0,1]$ — the
spectral-imaging analogue of the cosine/contrast-angle similarity long used in mass-spectral
library search (Stein and Scott, 1994; Wan *et al.*, 2002). Mean-centring is the discriminating
step: a cell where one spectrum has a peak and the other is empty contributes a negative product
and lowers the score, so intensity that is not co-located is actively penalized rather than, as in
the nearest-neighbour distance, ignored. No shift search is performed — aligning the images (as in
NMR peak-alignment tools such as icoshift; Savorani *et al.*, 2010) would let a different spectrum
slide into registration and re-saturate, the failure mode of the tree and nearest-neighbour
methods.

Self-similarity is exactly one (as it is for all methods tested); the score is symmetric exactly,
and monotone in the drift of an isolated peak (verified in the unit tests). The implementation is
pure NumPy (Harris *et al.*, 2020), reads processed Bruker data directly (cf. Helmus and Jaroniec,
2013), needs no peak picking and no external solver, and runs in about 6 ms per spectrum pair on a
$250 \times 350$ grid (Apple Silicon, single core); a 19-test suite (all passing) checks
self-similarity, related-beats-unrelated ordering and monotonicity in drift. Figures use
Matplotlib (Hunter, 2007). A full derivation with proofs is provided as Supplementary Material.

## 3 Results

We benchmarked LCC against the three reference methods on **two independent, opposite regimes**.
Discrimination is summarized by the **separation** (mean same-class score − mean different-class
score) and the **margin** (worst same-class − best different-class score, what a classification
threshold sees); both are on each method's own $[0,1]$ scale (Table 1; Fig. 1). All methods
self-score exactly 1.

**Dense protein $^1$H-$^{15}$N.** The reference is the apo spectrum of the phosphatase PRL3,
compared with six spectra of the same protein along a ligand titration (should score high) and
with a **different** protein (should score low), inside $^1$H 6.5–10 $\times\ ^{15}$N 105–130 ppm.
LCC scores the same protein 0.94 (range 0.89–0.97) and the different protein 0.18, a margin of
0.71 — three times the bin method's 0.24, the strongest baseline. Mean-centring earns that gap:
an un-centred cosine correlation on the same images drops the separation from 0.75 to 0.59. The
gain is not a tuned parameter (it holds across the physical blur range; Supplementary).

**Sparse small-molecule $^1$H-$^{13}$C.** To test the opposite regime — and to replace the single
negative spectrum above with a large negative set — we used public HSQC peak lists (simpleNMR;
Data availability) for six compounds each recorded twice, and asked each method to score the 6
same-compound pairs above the 60 different-compound pairs, inside $^1$H 0–10 $\times\ ^{13}$C
0–165 ppm. LCC again leads (mean same 0.82, mean different 0.003, separation 0.81). The bin method
is far stronger here (0.69) than on the dense fingerprint (0.29) — with few, well-separated peaks
its hard bins rarely straddle a peak — while the tree and nearest-neighbour methods **leave the two
classes barely separable in this regime too** (margins 0.00 and 0.04, despite separations of 0.39
and 0.10): over a wide window a different molecule still has a near neighbour for every peak. Their
shift tolerance is genuine (per-pair, it keeps two solvent-shifted same-compound recordings high,
0.83–1.00, where the position-sensitive LCC and bin methods drop to 0.30–0.60) but it tolerates
*everything*, so it does not become discrimination; LCC recovers that tolerance controllably
through its blur width.

Table: Separation and worst-case margin for all methods on both regimes, computed from unrounded
scores. Every method self-scores 1.00.

| Method | $^1$H-$^{15}$N sep | margin | $^1$H-$^{13}$C sep | margin |
| --- | --- | --- | --- | --- |
| Bin (Bodis, 2009) | 0.29 | 0.24 | 0.69 | 0.20 |
| Bin + 45° rotation | 0.25 | 0.20 | 0.69 | 0.18 |
| Quad-tree (Castillo, 2013) | 0.03 | −0.01 | 0.39 | 0.00 |
| Nearest neighbour (Pierens, 2012) | 0.04 | 0.03 | 0.10 | 0.04 |
| **LCC (this work)** | **0.75** | **0.71** | **0.81** | **0.37** |

**LCC leads in both regimes** — the only method that keeps same-class similarity high while
pushing different-class similarity down, whether the fingerprint is a dense protein amide region
or a handful of scattered small-molecule crosspeaks. Its two knobs are the per-axis blur widths,
both physical: the linewidth sets the floor and the expected shift perturbation sets how much
drift is tolerated. LCC measures position-and-intensity coincidence, not molecular identity, so a
close analogue will score intermediate — correct for a similarity, a caveat for hard
classification. The three reference methods are our own reimplementations, and the $^1$H-$^{13}$C
spectra are rendered from peak lists inside one common window, so the *ranking*, not the exact
values, is the result.

![**Figure 1.** Same/different separation (mean same-class − mean different-class similarity) for
all five methods on the two benchmarks: dense protein $^1$H-$^{15}$N (blue) and sparse
small-molecule $^1$H-$^{13}$C (orange). LCC is the only method with a large separation in both
regimes; the tree and nearest-neighbour methods saturate in both.](../results/comparison_all.png)

## Data availability

The dense $^1$H-$^{15}$N benchmark uses PRL3 and OAA Bruker spectra available from the author on
reasonable request. The sparse $^1$H-$^{13}$C benchmark uses public HSQC peak lists from the
simpleNMR example set (Hughes and Kenwright, 2024;
https://github.com/EricHughesABC/simpleNMR), downloaded automatically by `bench_13c.py`. Both harnesses (`bench.py`, `bench_13c.py`) and the tabulated scores are in the
repository.

## Funding

None declared.

## Conflict of interest

None declared.

## References

Bingol,K. and Brüschweiler,R. (2014) Multidimensional approaches to NMR-based metabolomics.
*Anal. Chem.*, **86**, 47–57.

Bodis,L. *et al.* (2007) A novel spectra similarity measure. *Chemometr. Intell. Lab. Syst.*,
**85**, 1–8.

Bodis,L. *et al.* (2009) Automatic compatibility tests of HSQC NMR spectra with proposed
structures of chemical compounds. *Talanta*, **79**, 1379–1386.

Castillo,A.M. *et al.* (2013) Fast and shift-insensitive similarity comparisons of NMR using a
tree-representation of spectra. *Chemometr. Intell. Lab. Syst.*, **127**, 1–6.

Craig,A. *et al.* (2006) Scaling and normalization effects in NMR spectroscopic metabonomic
data sets. *Anal. Chem.*, **78**, 2262–2267.

Delaglio,F. *et al.* (1995) NMRPipe: a multidimensional spectral processing system based on
UNIX pipes. *J. Biomol. NMR*, **6**, 277–293.

Dieterle,F. *et al.* (2006) Probabilistic quotient normalization as robust method to account for
dilution of complex biological mixtures. Application in 1H NMR metabonomics. *Anal. Chem.*,
**78**, 4281–4290.

Erlanson,D.A. *et al.* (2016) Twenty years on: the impact of fragments on drug discovery.
*Nat. Rev. Drug Discov.*, **15**, 605–619.

Harris,C.R. *et al.* (2020) Array programming with NumPy. *Nature*, **585**, 357–362.

Helmus,J.J. and Jaroniec,C.P. (2013) Nmrglue: an open source Python package for the analysis of
multidimensional NMR data. *J. Biomol. NMR*, **55**, 355–367.

Hughes,E. and Kenwright,A.M. (2024) SimpleNMR: an interactive graph network approach to aid
constitutional isomer verification using standard 1D and 2D NMR experiments. *Magn. Reson. Chem.*,
**62**, 556–565.

Hunter,J.D. (2007) Matplotlib: a 2D graphics environment. *Comput. Sci. Eng.*, **9**, 90–95.

Ludwig,C. and Günther,U.L. (2011) MetaboLab — advanced NMR data processing and analysis for
metabolomics. *BMC Bioinformatics*, **12**, 366.

Pellecchia,M. *et al.* (2008) Perspectives on NMR in drug discovery: a technique comes of age.
*Nat. Rev. Drug Discov.*, **7**, 738–745.

Pierens,G.K. *et al.* (2012) HSQC spectral based similarity matching of compounds using nearest
neighbours and a fast discrete genetic algorithm. *J. Cheminform.*, **4**, 25.

Robinette,S.L. *et al.* (2012) NMR in metabolomics and natural products research: two sides of
the same coin. *Acc. Chem. Res.*, **45**, 288–297.

Savorani,F. *et al.* (2010) icoshift: a versatile tool for the rapid alignment of 1D NMR spectra.
*J. Magn. Reson.*, **202**, 190–202.

Shuker,S.B. *et al.* (1996) Discovering high-affinity ligands for proteins: SAR by NMR.
*Science*, **274**, 1531–1534.

Skinner,S.P. *et al.* (2016) CcpNmr AnalysisAssign: a flexible platform for integrated NMR
analysis. *J. Biomol. NMR*, **66**, 111–124.

Stein,S.E. and Scott,D.R. (1994) Optimization and testing of mass spectral library search
algorithms for compound identification. *J. Am. Soc. Mass Spectrom.*, **5**, 859–866.

Steinbeck,C. and Kuhn,S. (2004) NMRShiftDB — compound identification and structure elucidation
support through a free community-built web database. *Phytochemistry*, **65**, 2711–2717.

Ulrich,E.L. *et al.* (2008) BioMagResBank. *Nucleic Acids Res.*, **36**, D402–D408.

Vranken,W.F. *et al.* (2005) The CCPN data model for NMR spectroscopy: development of a software
pipeline. *Proteins*, **59**, 687–696.

Wan,K.X. *et al.* (2002) Comparing similar spectra: from similarity index to spectral contrast
angle. *J. Am. Soc. Mass Spectrom.*, **13**, 85–88.

Williamson,M.P. (2013) Using chemical shift perturbation to characterise ligand binding.
*Prog. Nucl. Magn. Reson. Spectrosc.*, **73**, 1–16.

Wishart,D.S. *et al.* (2022) HMDB 5.0: the Human Metabolome Database for 2022. *Nucleic Acids
Res.*, **50**, D622–D631.

Worley,B. and Powers,R. (2013) Multivariate analysis in metabolomics. *Curr. Metabolomics*,
**1**, 92–107.
