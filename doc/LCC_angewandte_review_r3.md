---
title: "Simulated Peer Review — Round 3 (two-decoy dense benchmark): A Lineshape Correlation Coefficient for 2D HSQC Similarity"
author: "Editorial simulation (4 referees + adversarial verification + handling editor)"
date: " "
geometry: margin=1in
colorlinks: true
linkcolor: RoyalBlue
urlcolor: RoyalBlue
fontsize: 11pt
---

# Simulated peer review — Round 3

**Manuscript:** *A Lineshape Correlation Coefficient for Robust Similarity of Two-Dimensional HSQC
NMR Spectra* (D. Lee, KBSI), re-reviewed after the addition of an **expanded dense $^1$H–$^{15}$N
benchmark with two decoy proteins** (main-text Table 2; SI "Expanded benchmark — a second decoy
protein", Tables S8–S9). Method: four independent referees, adversarial verification of every MAJOR
against the actual files, then a handling-editor re-decision. This round specifically asks whether
the new data changes the Round-2 disposition (reject-and-transfer).

---

## Round-3 decision: **REJECT (reject-and-transfer)** — unchanged outcome, one round-2 objection now closed

| Referee | Lens | Round 2 | Round 3 |
| --- | --- | :---: | :---: |
| Referee 1 | Novelty, significance, journal fit | Reject | **Reject** |
| Referee 2 | Technical correctness & statistics | Major revision | **Reject-and-transfer** |
| Referee 3 | Benchmarking fairness & reproducibility | Major revision | **Major revision** |
| Referee 4 | Presentation & format | Reject | **Reject** |
| **Handling editor** | **Consolidated** | **Reject-and-transfer** | **Reject-and-transfer** |

**Verification outcome:** of **14** MAJOR findings raised this round, adversarial verification
CONFIRMED **5** and judged the rest OVERSTATED/REFUTED. Crucially, **three of the five confirmed
findings were newly introduced by this revision** and are fixable defects (now corrected — see
below); the other two are the standing significance/reproducibility items.

---

## What the new data changed (verified)

The two-decoy benchmark is a **genuine, reproducible strengthening** that closes one round-2 item:

- **Round-2 cap "the mean-centring novelty and the dense margin rest on a single negative spectrum
  (n=1)" is materially resolved.** The negative side is now a distribution over **two distinct
  folds** — OAA (277 K) and the kinase domain EphB3 (900 MHz) — scored against an **independent
  23-point PRL3 apo→olsalazine titration**. LCC keeps both decoys below *all 23* same-protein points
  (worst same 0.82 > worst decoy 0.41; margin +0.41), mean-centring still earns the gap against a
  real negative distribution (LCC separation 0.57 vs un-centred cosine 0.46), and the titration
  decays smoothly 0.98→0.82 — a clean demonstration of the graded shift tolerance on real data.
- **Referees 2 and 3 independently re-ran `bench_nhsqc.py`**: Table 2 and Tables S8/S9 reproduce to
  the digit, the separation/margin definitions are correct and reduce to the Table-1 definitions at
  one decoy, the reference is genuinely apo, and the decoys carry real overlapping amide structure
  (tree/NN saturate at 0.75–0.96) that LCC nonetheless separates — so it is not an empty-window
  artefact.

**Why the outcome is unchanged.** Round 2 had *already* down-ranked the single-negative concern to
MINOR; its **decisive** gate was *broad chemical significance for a general journal*, closeable only
by "at least one real chemical result where LCC changes a conclusion an incumbent gets wrong." Table 2 adds
more same/different **discrimination scores** (23 comparisons + 2 decoys) — the same evidence type at
larger $n$ — not an application: no binder identified, no CSP site mapped, no unknown dereplicated,
no QC/comparability failure caught. The HOS framing (refs 11–16) is well-cited motivation, not a
demonstrated comparability result. So rigor is up (helping the specialist-venue case) but the
general-journal gate is untouched.

---

## Confirmed findings (survived adversarial verification)

**Newly introduced by this revision — now fixed:**

1. **[MAJOR, fixed] Factual overclaim "the only margin."** Main text called LCC's the "only margin
   that pushes both decoy proteins below every same-protein score," yet all six methods have positive
   Table S8 margins (Bin 0.25, Bin+45 0.25, tree 0.12, NN 0.02, cosine 0.35, LCC 0.41) and the same
   sentence cites the positive tree/NN margins. **Fixed:** reworded to "the widest margin (0.41)…
   whereas the tree and nearest-neighbour margins are a knife-edge (0.12 and 0.02)."
2. **[MAJOR, fixed] Reproducibility regression in `bench_nhsqc.py`.** The new harness hard-coded the
   private `Nhsqc/` paths with no flag/env interface and a raw `FileNotFoundError`, undoing the
   round-2 `--prl3/--oaa` + `$PRL3_DIR` fix. **Fixed:** added `--titration/--oaa/--ephb3` flags and
   `$NHSQC_*` env overrides with a clean `p.error()` when the data are absent (numbers unchanged).
3. **[MAJOR, fixed] SI §S13 de-synchronised.** §S13 still called the dense margin "a single negative
   spectrum (§S8)" — contradicting §S8, where the second decoy was added. **Fixed:** §S13 now names
   the two-decoy expansion and the 0.41 margin.

**Standing items (decisive / editorial — not fixable by editing):**

4. **[MAJOR → decisive] No application-grade chemical result.** The reject-driving gate; unchanged
   since round 2. Needs new experimental chemistry, not more discrimination numbers.
5. **[MAJOR → minor at a methods venue] The decisive dense evidence is not yet public.**
   `bench_nhsqc.py`, the `Nhsqc/` raw data and `results/nhsqc_dense.json` are untracked in git and
   the deposition is "in preparation" (no DOI), so Table 2's "reproducible" claim is false on a
   public checkout. Fix by committing the harness + depositing the matrices with a DOI.

---

## Remaining blockers (ranked)

1. **Chemical significance for a general journal** — add at least one application where LCC corrects an
   incumbent's chemical conclusion (a blind BMRB/HMDB dereplication hit at library scale under a
   rejection threshold, or a CSP/biosimilar-HOS screen where bin/tree/NN misrank and LCC ranks
   correctly). New application-grade work → resubmit, not held revision.
2. **Deposit the dense data** (both single- and two-decoy) with a DOI; commit `bench_nhsqc.py` so the
   in-repo claim is true on a clean checkout.
3. **Carryover minors** (unchanged from round 2): unlabelled 3× (margin) vs 2.6× (separation)
   multiplier; the per-regime grid step as a de-facto third parameter despite "set by spectroscopy";
   the retrieval null (top-1 = 1.00 tie) as the only task-level test.

---

## Bottom line

The two-decoy benchmark is a **real success on the axis round 2 asked about**: the dense result no
longer rests on a single negative spectrum, mean-centring's benefit is demonstrated against two
distinct folds and a 23-point real titration, and everything reproduces locally to the digit. It is
**not** a success on the axis that decides venue: no chemical result was added, so the
general-chemistry significance gate is untouched, and the revision briefly introduced three fixable
defects (all now corrected). The intellectually honest disposition is unchanged and, if anything,
sharper: **transfer to a specialized methods journal** (J. Magn. Reson., Anal. Chem., or
J. Cheminform.), where this is now an even stronger, near-acceptable methods paper, with a path back
to Angewandte only on a genuine, at-scale chemical demonstration.
