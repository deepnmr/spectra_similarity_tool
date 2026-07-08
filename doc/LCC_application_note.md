---
title: "LCC: a lineshape-correlation coefficient for discriminating dense protein $^{1}$H-$^{15}$N HSQC fingerprints"
author: "Donghan Lee"
date: " "
geometry: margin=1in
colorlinks: true
---

Korea Basic Science Institute (KBSI), Ochang, Republic of Korea. Correspondence:
kbsi.bionmr@gmail.com

**Subject area:** Structural bioinformatics / NMR spectroscopy

## Abstract

**Motivation:** Quantitative comparison of two-dimensional $^1$H-$^{15}$N HSQC spectra
underlies ligand screening and protein fingerprinting, yet the established similarity
measures were designed for sparse small-molecule $^1$H-$^{13}$C spectra. On the dense
amide fingerprint of a protein they fail in two opposite ways: bin-histogram methods
lose margin under the small peak drifts that a titration induces, while tree and
nearest-neighbour methods saturate high (different-protein score 0.87–0.96) because in
a crowded fingerprint every peak finds a coincidental near neighbour, leaving almost no
contrast between the same protein and a different one.

**Results:** We introduce the Lineshape Correlation Coefficient (LCC), which renders
each spectrum onto a shared grid, blurs it by the physical linewidth plus the expected
shift drift, and scores the two images by their mean-centred normalized cross-correlation
at zero lag. The Gaussian render supplies graded shift tolerance while mean-centring
penalizes intensity that is not co-located, so a titrated ligand series stays high while a
different protein decorrelates; removing the mean-centring drops the separation from 0.75
to 0.59, confirming it as the discriminating step. On a $^1$H-$^{15}$N HSQC benchmark of
six same-protein titration points and one different-protein spectrum, LCC widens the
worst-case margin between the two classes to 0.71 (versus 0.24 for the best prior method);
the between-protein contrast rests on a single negative spectrum and should be read as
indicative pending a larger negative set.

**Availability and implementation:** Freely available under the MIT license at
https://github.com/deepnmr/spectra_similarity_tool; implemented in Python 3 (NumPy;
Matplotlib optional for figures); platform independent. The repository includes the three
reference methods, a benchmark harness and a test suite.

**Contact:** kbsi.bionmr@gmail.com

## 1 Introduction

A recurring task in biomolecular NMR is to score how similar two 2D HSQC spectra are:
tracking a $^1$H-$^{15}$N amide fingerprint across a ligand titration, mapping a chemical
shift perturbation on binding (Shuker *et al.*, 1996; Pellecchia *et al.*, 2008;
Williamson, 2013), or matching a spectrum against a spectral library (Ulrich *et al.*,
2008). Three families of measure are in use. The bin method of Bodis
*et al.* integrates each spectrum over an $n \times n$ grid at a range of resolutions and
reports a weighted Jaccard (Ružička) index (Bodis *et al.*, 2007, 2009). The tree method of
Castillo *et al.* encodes a spectrum as a recursive centre-of-mass quad-tree and compares
the trees node-by-node with a shift-tolerant node score (Castillo *et al.*, 2013). The peak
method of Pierens *et al.* picks peaks and matches each to its nearest neighbour in the
other spectrum (Pierens *et al.*, 2012).

All three were developed for **sparse** small-molecule $^1$H-$^{13}$C HSQC spectra. The
amide fingerprint of a folded protein is the opposite regime: on the order of a hundred
crosspeaks (roughly one per residue) crowded into a single $^1$H-$^{15}$N region. Here the
measures fail in two opposite ways. Bin histograms are shift-**brittle**: a peak straddling
a bin edge splits its intensity across neighbouring bins, so the small drifts a titration
induces move mass between bins and erode the score. The tree and nearest-neighbour measures
are shift-**blind**: in a crowded fingerprint every peak has some near neighbour and every
spectrum has a similar mass-centre structure, so both saturate high and retain almost no
contrast (separation 0.03–0.04 below). A measure for this regime must tolerate small physical
drifts yet still penalize genuinely relocated intensity.

## 2 Implementation

LCC keeps the physically meaningful ingredients of the three methods — area-weighted
integration and unit normalization (Bodis), a physical shift tolerance (Castillo), a
lineshape picture of the spectrum (Pierens) — and adds a mean-centred correlation as the
discriminating step. Given two processed Bruker spectra — produced by any standard pipeline
or analysis suite (Delaglio *et al.*, 1995; Vranken *et al.*, 2005; Skinner *et al.*, 2016) —
compared inside a common window, it proceeds in three steps.

**Render.** Negative intensities are clipped, each point is weighted by its local integration
area, and both spectra are histogrammed onto a **single shared grid** so they are directly
comparable irrespective of their native digital resolution.

**Blur.** Each image is convolved with a separable Gaussian whose width per axis is the
physical linewidth combined with the expected drift, $\sigma = \sqrt{\ell^2 + d^2}$ (defaults
$\sigma_{^1\mathrm{H}} = 0.03$, $\sigma_{^{15}\mathrm{N}} = 0.30$ ppm). This replaces the bin
method's hard edges with a continuous response: for an isolated peak on a large grid (where
mean-centring is negligible) two copies displaced by $(\Delta_H, \Delta_N)$ correlate as
$\exp(-\Delta_H^2/4\sigma_H^2 - \Delta_N^2/4\sigma_N^2)$, smooth and monotone in the drift, so
a physical titration shift costs little and a random relocation costs a lot.

**Score.** Each image is mean-centred and the similarity is the zero-lag normalized
cross-correlation (the Pearson coefficient) of the two centred images, clamped to $[0,1]$.
Mean-centring is the discriminating step: a cell where one spectrum has a peak and the other
is empty contributes a negative product and lowers the score, so intensity that is not
co-located is actively penalized rather than, as in the nearest-neighbour distance, ignored.
No shift search is performed — aligning the images would let a different protein slide into
registration and re-saturate, the very failure mode of the tree and nearest-neighbour methods.

Self-similarity is exactly one (as it is for all methods tested); the score is symmetric
exactly, and monotone in the drift of an isolated peak (verified in the unit tests). The
implementation is pure NumPy (Harris *et al.*, 2020), reads processed Bruker data directly
(cf. Helmus and Jaroniec, 2013), needs no peak picking and no external solver, and runs in
about 6 ms per spectrum pair on a $250 \times 350$ grid (Apple Silicon, single core); a
17-test suite (all passing) checks self-similarity, related-beats-unrelated ordering and
monotonicity in drift. Figures use Matplotlib (Hunter, 2007). A full derivation with proofs
is provided as Supplementary Material.

## 3 Results

We benchmarked LCC against the three reference methods on a protein $^1$H-$^{15}$N HSQC data
set. The reference is the apo spectrum of the phosphatase PRL3, compared with six spectra of
the same protein along a ligand titration (which should score high) and with the spectrum of
a **different** protein (which should score low). Because the six positives are titration
points of one protein–ligand system they form a correlated series, not independent replicates,
and the negative class is a single spectrum; the numbers below are therefore descriptive of
this set rather than population estimates. All spectra were compared inside $^1$H 6.5–10 ppm
$\times\ ^{15}$N 105–130 ppm. We summarize discrimination by the separation (mean same-protein
score minus different-protein score) and, because it is what a classification threshold sees,
by the worst-case margin (minimum same-protein score minus different-protein score); both are
on each method's own $[0,1]$ score scale (Table 1; Fig. 1).

Table: Similarity of the reference PRL3 spectrum to the same protein (mean over six titration
points) and to a different protein. Separation and margin are computed from unrounded scores.
All methods self-score 1.00.

| Method | mean same | different | separation | margin |
| --- | --- | --- | --- | --- |
| Bin (Bodis, 2009) | 0.79 | 0.49 | 0.29 | 0.24 |
| Bin + 45° rotation | 0.82 | 0.57 | 0.25 | 0.20 |
| Quad-tree (Castillo, 2013) | 0.90 | 0.87 | 0.03 | −0.01 |
| Nearest neighbour (Pierens, 2012) | 0.99 | 0.96 | 0.04 | 0.03 |
| **LCC (this work)** | **0.94** | **0.18** | **0.75** | **0.71** |

The six same-protein scores span 0.89–0.97 (mean 0.94) and the different protein scores 0.18,
so LCC opens a worst-case margin of 0.71 — three times the 0.24 of the bin method, the
strongest prior baseline, and well clear of the tree and nearest-neighbour methods, which
saturate as expected for this regime. Both LCC and the bin method separate this particular
seven-spectrum set with a single threshold; the practical difference is the width of the gap,
which governs robustness to an unseen spectrum. The mean-centring is what earns that gap:
replacing it with an un-centred cosine correlation on the same images drops the separation
from 0.75 to 0.59 and the margin from 0.71 to 0.55. The gain is also not an artefact of one
tuned parameter — LCC beats the bin method across the whole physical blur range (separation
0.71–0.77 for $\sigma_{^1\mathrm{H}}$ 0.02–0.04 and $\sigma_{^{15}\mathrm{N}}$ 0.20–0.40 ppm)
and still wins (0.36) when coarsened to the bin method's own resolution.

The two knobs are the per-axis blur widths, both physical: the NMR linewidth sets the floor
and the expected chemical shift perturbation sets how much drift is tolerated; widening
$\sigma_{^{15}\mathrm{N}}$ trades separation for robustness when strong binding-site shifts are
anticipated. LCC measures position-and-intensity coincidence rather than protein identity, so a
homologue sharing the amide envelope will score intermediate — correct behaviour for a
similarity but a caveat for classification. LCC is intended for dense protein fingerprints and
titration tracking; the tree and nearest-neighbour methods remain preferable for sparse
small-molecule spectra where shift insensitivity is the goal.

Table 1 is regenerated by `bench.py --prl3 <dir> --oaa <dir> --json` (Python 3.11, NumPy
2.4.6); the full per-comparison scores and Fig. 1 are committed under `results/`.

![**Figure 1.** Per-comparison LCC similarity of the reference PRL3 $^1$H-$^{15}$N HSQC
spectrum. Left: similarity along the six-point ligand titration of the same protein (green,
should score high, 0.89–0.97) and to a different protein (OAA, red, 0.18), for all five
methods; LCC (purple) is the only measure that keeps the titration high while dropping the
different protein far below it. Right: the resulting separation (mean same − different) per
method.](../results/lcc_comparison.png)

## Data availability

The processed PRL3 and OAA Bruker spectra used for the benchmark are available from the author
on reasonable request. The benchmark harness accepts their locations via `--prl3`/`--oaa`.

## Funding

None declared.

## Conflict of interest

None declared.

## References

Bodis,L. *et al.* (2007) A novel spectra similarity measure. *Chemometr. Intell. Lab. Syst.*,
**85**, 1–8.

Bodis,L. *et al.* (2009) Automatic compatibility tests of HSQC NMR spectra with proposed
structures of chemical compounds. *Talanta*, **79**, 1379–1386.

Castillo,A.M. *et al.* (2013) Fast and shift-insensitive similarity comparisons of NMR using a
tree-representation of spectra. *Chemometr. Intell. Lab. Syst.*, **127**, 1–6.

Delaglio,F. *et al.* (1995) NMRPipe: a multidimensional spectral processing system based on
UNIX pipes. *J. Biomol. NMR*, **6**, 277–293.

Harris,C.R. *et al.* (2020) Array programming with NumPy. *Nature*, **585**, 357–362.

Helmus,J.J. and Jaroniec,C.P. (2013) Nmrglue: an open source Python package for the analysis of
multidimensional NMR data. *J. Biomol. NMR*, **55**, 355–367.

Hunter,J.D. (2007) Matplotlib: a 2D graphics environment. *Comput. Sci. Eng.*, **9**, 90–95.

Pellecchia,M. *et al.* (2008) Perspectives on NMR in drug discovery: a technique comes of age.
*Nat. Rev. Drug Discov.*, **7**, 738–745.

Pierens,G.K. *et al.* (2012) HSQC spectral based similarity matching of compounds using nearest
neighbours and a fast discrete genetic algorithm. *J. Cheminform.*, **4**, 25.

Shuker,S.B. *et al.* (1996) Discovering high-affinity ligands for proteins: SAR by NMR.
*Science*, **274**, 1531–1534.

Skinner,S.P. *et al.* (2016) CcpNmr AnalysisAssign: a flexible platform for integrated NMR
analysis. *J. Biomol. NMR*, **66**, 111–124.

Ulrich,E.L. *et al.* (2008) BioMagResBank. *Nucleic Acids Res.*, **36**, D402–D408.

Vranken,W.F. *et al.* (2005) The CCPN data model for NMR spectroscopy: development of a software
pipeline. *Proteins*, **59**, 687–696.

Williamson,M.P. (2013) Using chemical shift perturbation to characterise ligand binding.
*Prog. Nucl. Magn. Reson. Spectrosc.*, **73**, 1–16.
