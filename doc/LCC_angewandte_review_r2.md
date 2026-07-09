---
title: "Simulated Peer Review — Round 2 (revised manuscript): A Lineshape-Correlation Coefficient for 2D HSQC Similarity"
author: "Editorial simulation (4 referees + adversarial verification + handling editor)"
date: " "
geometry: margin=1in
colorlinks: true
linkcolor: RoyalBlue
urlcolor: RoyalBlue
fontsize: 11pt
---

# Simulated peer review — Round 2

**Manuscript:** revised *A Lineshape-Correlation Coefficient for Robust Similarity of
Two-Dimensional HSQC NMR Spectra* (D. Lee, KBSI). This round re-reviews the manuscript after the
revisions made in response to Round 1 (see `LCC_angewandte_review.md`). Same method: four
independent referees (each carrying their Round-1 lens), adversarial verification of every
FATAL/MAJOR finding against the revised text, then a handling-editor re-decision.

---

## Round-2 decision: **REJECT (reject-and-transfer)** — unchanged outcome, materially stronger manuscript

| Referee | Lens | Round 1 | Round 2 |
| --- | --- | :---: | :---: |
| Referee 1 | Novelty, significance, journal fit | Reject | **Reject** |
| Referee 2 | Technical correctness & statistics | Major revision | **Major revision** |
| Referee 3 | Benchmarking fairness & reproducibility | Major revision | **Major revision** |
| Referee 4 | Presentation & Angewandte format | Minor revision | **Reject** |
| **Handling editor** | **Consolidated** | **Reject** | **Reject-and-transfer** |

**Verification outcome:** of the **10** MAJOR findings raised this round, adversarial verification
judged **all 10 OVERSTATED** and none CONFIRMED as a rebuttable defect (Round 1: 2 confirmed / 14
overstated). **The revised manuscript is technically clean** — no MAJOR survives scrutiny. The
outcome is unchanged only because the decisive gate is not technical: it is **chemical
significance for a general journal**, which the revision *conceded* rather than closed.

### The honest paradox of this revision

The revision did exactly what Round 1 asked on every achievable axis — and in doing so it made its
own ceiling explicit. Two of the paper's *new*, honest results are the ones that now cap it:

- **The added dereplication/retrieval "application" is a null result:** every method (LCC, bin,
  tree, NN, cosine) scores top-1 = 1.00 / MRR = 1.00 on the clean 5-compound library. It documents
  the *absence* of a task-level advantage, so it cannot lift significance.
- **The added bootstrap + cosine ablation show the sole novel ingredient is a no-op on the
  reproducible public benchmark:** un-centred cosine 0.777 ≈ LCC 0.778 on sparse $^1$H–$^{13}$C, and
  LCC's edge over the bin method there is *within noise* (overlapping CIs). Mean-centring's only
  demonstrated benefit (0.59 → 0.75) lives entirely on the **dense** benchmark, whose separation
  rests on a **single** negative spectrum.

So the more rigorously and honestly the paper was revised, the clearer it became that its
conceptual advance is demonstrable only where the evidence is weakest (n = 1). This is not a defect
to fix by editing — it needs new, application-grade chemistry.

### Editor's rationale

For a top general-chemistry journal the binding constraint is broad chemical significance, and the
revision answers it by *agreeing* with the first-round objection rather than overturning it: no
binder identified, no CSP site mapped, no unknown dereplicated, no QC failure caught; the result
set is still 7 protein comparisons plus 45 small-molecule pairs. The one experiment added as a
chemical application is an honest null. Meanwhile the sole novel ingredient (mean-centring) is a
no-op on the reproducible public benchmark, so its only demonstrated benefit lives against a single
dense negative spectrum, and the flagship number is still not depositable (a DOI-less "deposition
in preparation"). Nothing was refuted this round, so the three technical MAJORs (dense
reproducibility, mean-centring-rests-on-n=1, retrieval null) are down-ranked to MINOR as
disclosed-and-bounded, and the significance finding is kept as the one that decides venue fit. A
candidly-scoped similarity primitive with no chemical result, whose decisive number is n = 1 and
whose novelty is a no-op on the only public benchmark, is below Angewandte's bar. It is not fixable
by editing what is here; it needs a different, application-grade result — a **reject** for this
journal, **not** a held major revision. In fairness to the genuine improvement, the verifier
consensus is that at a specialist venue this manuscript is **already at minor-revision-accept
quality**, so the correct disposition is **reject-and-transfer** (J. Magn. Reson., Anal. Chem., or
J. Cheminform.), with a stated path back to Angewandte only if a real chemical demonstration is
added.

---

## Round-1 issues now resolved (verified)

The revision cleared, and verifiers independently re-ran/reproduced, the following Round-1 items:

- **Reproducibility of the code:** `bench.py` no longer hard-codes private absolute paths — it reads
  `--prl3/--oaa` or `$PRL3_DIR/$OAA_DIR` and errors cleanly when absent (Referee 3, verified).
- **Cosine ablation reproducible:** the mean-centring ablation is now a benchmark row re-run and
  reproduced to the digit by Referees 2 and 3 (sparse 0.777), not an asserted number.
- **Uncertainty quantification:** bootstrap 95% CIs are computed correctly and used *honestly* —
  overlapping CIs decline a claim ("within noise") rather than assert one.
- **Data-integrity bug fixed:** the byte-identical olivetol duplicate is removed, guarded by a
  raising `_assert_distinct` + regression test, and its removal correctly *lowered* the sparse
  separation (0.81 → 0.78) rather than flattering the result.
- **Fairness wording:** "faithful reimplementations" softened to "representative"; self-similarity
  = 1 no longer sold as an LCC advantage ("as for every method tested").
- **Sparse benchmark fully public & reproducible end-to-end** (auto-downloads the peak lists), plus
  the retrieval/MRR test; **26/26 tests pass** (was 19).
- **Presentation cluster largely cleared** and significance overclaiming replaced by honest
  self-scoping as a methodological contribution.

---

## What remains (ranked)

1. **[MAJOR → decisive] No demonstrated chemical significance for a general journal** — *conceded,
   not closed.* The only added "application" is a retrieval null (all methods tie at top-1 = 1.00).
   This is the reject-driving gate and an editorial call, not a technical flaw.
2. **[MAJOR → MINOR by verification] Flagship dense benchmark still not independently reproducible**
   — the most decisive number (margin 0.71) rests on request-only data with a DOI-less "deposition
   in preparation"; should be depositable *now*.
3. **[MAJOR → MINOR by verification] The novel step is unpowered** — mean-centring is a no-op on the
   reproducible public benchmark (0.778 vs 0.777); its benefit (0.59 → 0.75) rests on a single dense
   negative spectrum.
4. **[MINOR]** Dense cosine 0.59 is a hard-coded literal in `make_plots.py`/`comparison_all.json`,
   not read from a committed harness output (the stale `method_comparison.json` lacks the row) — an
   artefact-hygiene gap.
5. **[MINOR]** Abstract juxtaposes "decisive" (the n = 1 dense result) with "within noise" (the
   quantified sparse result) and omits the n = 1 caveat present in the body; "widest separation in
   both regimes" sits inside the bin CI on the sparse set.
6. **[MINOR carryovers]** Inconsistent 3× (margin) vs 2.6× (separation) multiplier, neither
   labelled; the per-regime grid step (0.01/0.10 → 0.02/0.2) is a de-facto third, non-spectroscopic
   parameter despite the "grid-invariant" framing; the bootstrap treats the 40 different-compound
   pairs as independent (they are C(10,2) of 10 spectra) — conservative but undisclosed;
   method-name hyphenation inconsistent; Table 1 leaves the dense cosine margin as "—".

---

## Revision plan (editor)

1. **Gate for Angewandte (only thing that can change the outcome):** add ≥1 real chemical result
   where LCC changes a conclusion an incumbent gets wrong — a blind BMRB/HMDB dereplication hit run
   at **library scale under a rejection threshold** (so separation converts into a retrieval
   difference the clean-library null cannot show), or a CSP/titration binding-rank screen where
   bin/tree/NN misrank and LCC ranks correctly. New application-grade work → reject-and-resubmit,
   not held major revision.
2. **Populate the dense benchmark with multiple decoy proteins** (varied folds/sizes) so margin 0.71
   and the 0.59 → 0.75 ablation gain a negative distribution + CI/ROC-AUC, validating the novel step
   where it actually operates instead of against one spectrum.
3. **Deposit** the processed/windowed PRL3+OAA matrices (Zenodo/BMRB) and cite the DOI; regenerate
   and commit `method_comparison.json` with the `cosine_uncentred` row so the dense 0.59 is read from
   a harness output, not a literal.
4. **Presentation fixes regardless of venue:** attach the n = 1 caveat to "decisive" in the abstract
   and reword the sparse lead as "ties the bin method within bootstrap noise, leads over tree/NN";
   label each multiplier's basis; fix the grid step by a stated rule or disclose it as a third
   parameter; note the different-pair non-independence; unify hyphenation; complete the Table 1 dense
   ablation margin.
5. **Strong transfer recommendation:** as a candidly-scoped, reproducible similarity primitive this
   manuscript is already at **minor-revision-accept** quality for a specialized methods venue, where
   the n = 1 / deposition items become ordinary minor revisions.

---

## Bottom line

The revision is a success as *engineering and scientific integrity* and a **non-success for the
Angewandte target**: every technical MAJOR collapsed under verification (10/10 overstated), the
paper is now honest, reproducible on its public legs, and statistically correct — but the same
honesty exposed that LCC's novelty is a no-op on the reproducible benchmark and its decisive result
is n = 1, so the general-chemistry significance gate cannot be closed by editing. The intellectually
honest disposition, echoed by the verifiers, is **transfer to a specialized methods journal**, where
this is already a strong, near-acceptable paper, and return to Angewandte only with a genuine,
at-scale chemical demonstration.

---

# Appendix: full per-referee Round-2 reports

### Referee 1 — REJECT

*Focus: Whether the revised manuscript now demonstrates broad chemical significance sufficient for Angewandte's general-chemistry bar, or whether the added retrieval test and honest methodological scoping instead confirm it as a specialized-methods contribution. First-round view: REJECT.*

**Assessment.** My first-round verdict was REJECT because the paper produced no chemical result, only discrimination scores, and its novelty was a recombination of three prior methods plus a textbook Pearson/contrast-angle score. The revision is genuinely improved — more honest, better-powered on the sparse side, cleaner code — but on my lens it moves in the wrong direction for Angewandte. The author now explicitly concedes 'this is a methodological contribution ... not a new chemical discovery' (line 191-192): my central first-round finding is not resolved, it is agreed with. The one new experiment framed as a chemical application — library retrieval/dereplication — is an honest null result in which every method, including the incumbents LCC claims to beat, scores top-1 = 1.00, so it demonstrates no advantage on a real task and cannot lift significance. Worse for the significance case, the paper's own new bootstrap shows LCC's sparse-regime edge over the bin method is 'within noise', and its own ablation shows that on the reproducible public benchmark the novel mean-centring step is a no-op (0.778 vs 0.777) — LCC there reduces exactly to the established mass-spec cosine. The decisive advantage now rests, by the author's own admission, on a single negative spectrum. This is a technically clean, well-scoped methods paper for J. Magn. Reson., Anal. Chem. or J. Cheminform.; it is below Angewandte's breadth/significance bar, and the honest self-scoping the revision adds makes that home even clearer.

**Prior concerns:** Neither of my two first-round MAJOR concerns is resolved. (1) 'No demonstrated chemical discovery' REMAINS and is now explicitly conceded by the author (lines 191-194); the added retrieval test does not supply a discovery — it is a null result. (2) 'Incremental recombination / textbook scoring step' REMAINS and is arguably sharpened: the new reproducible ablation shows the sole novel ingredient (mean-centring) contributes nothing on the public sparse benchmark (0.778 vs 0.777), so LCC there is identically the mass-spec contrast-angle cosine. Reproducibility (a secondary point I raised via benchmark thinness) is PARTIALLY improved — private paths de-hardcoded, sparse benchmark now public and duplicate-cleaned, bootstrap CIs added — but the flagship dense result is still 'from the author on reasonable request'.

**Findings:**

- **[MAJOR · REMAINING] Still no chemical result; the significance gap is now conceded, not closed** — *verification: OVERSTATED → MINOR*
  - *Location:* MS penultimate paragraph (lines 191-194); Conclusion; Table 1; SI S8, S9
  - The manuscript now states outright: 'this is a **methodological** contribution — a robust, reproducible similarity primitive, not a new chemical discovery — and its most decisive number, the dense-protein margin, rests on a single negative spectrum' (lines 191-193). No binder is identified, no CSP site mapped, no unknown dereplicated, no QC failure caught anywhere in MS or SI. The entire result set remains separation/margin numbers on 7 protein comparisons and 45 small-molecule pairs (Table 1; Tables S2, S3).
  - *Suggested fix:* For Angewandte, add at least one real chemical result where LCC changes a conclusion an incumbent gets wrong (a blind BMRB/HMDB dereplication hit, or a titration/CSP screen where LCC correctly ranks binding while bins/tree/NN fail). Absent that, the author's own scoping points to the correct venue: transfer directly to J. Magn. Reson., Anal. Chem., or J. Cheminform., where the current evidence base is already sufficient for a minor-revision accept.
  - *Verifier:* The referee quotes lines 191-193 verbatim and correctly: the manuscript states "this is a **methodological** contribution — a robust, reproducible similarity primitive, not a new chemical discovery — and its most decisive number, the dense-protein margin, rests on a single negative spectrum". So the referee did NOT miss a fix — the concession is real and the factual observations (no binder identified, no CSP mapped, no unknown dereplicated, no QC failure caught; result set = 7 protein comparisons + 45 small-molecule pairs in Table 1) are all accurate.

However, characterizing this as a MAJOR unrebuttable defect overstates it, for three reasons grounded in the revised text:

1. The manuscript does not merely concede — it adds compensating, quantitative, PUBLIC evidence precisely to offset the single-negative-spectrum weakness. Lines 169-171: "The dense-protein margin (0.71) remains LCC's strongest result, but it rests on a single negative spectrum; the sparse regime and this retrieval test broaden the evidence with a larger, public, fully-reproducible negative set." Lines 192-194 continue: "...which is why we add the larger public sparse benchmark and the retrieval test." The 40-pair different-compound negative set, the bootstrap CIs (lines 160-163), and the retrieval/MRR test (lines 164-168) are the closure the referee says is absent.

2. A "methodological contribution / similarity primitive" is a legitimate paper category. The referee's demand — identify a binder, map a CSP site, dereplicate an unknown, catch a QC failure — is a demand for a different paper (an application study), not correction of a flaw in this one. The manuscript honestly scopes itself and names that application as future work: lines 194-198 "...the natural next step is a prospective, at-scale dereplication study against a community library, for which LCC supplies the scoring primitive."

3. The internal consistency and honesty are intact: the caveats at lines 178-184 bound the claim (LCC measures position-and-intensity coincidence, not molecular identity; reference methods are reimplementations; 13C spectra rendered from peak lists). Nothing is overclaimed.

Thus the referee's underlying point is a significance/venue-fit judgment (does a candidly-scoped methods primitive clear the impact bar for this journal?), which is an editorial call, not a rebuttable technical defect. A competent author CAN rebut it: the paper is accurate, honestly bounded, and already broadens its evidence base with a public reproducible benchmark and a retrieval test rather than leaning solely on the single negative spectrum.
- **[MAJOR · NEW] The one added 'chemical application' is a null result that ties LCC with the methods it claims to beat** — *verification: OVERSTATED → MINOR*
  - *Location:* MS 'Statistical strength and a retrieval test' (lines 158-171); Table 1 note (lines 240-243); SI S13 (lines 416-426)
  - Cast as library search/dereplication, 'rank the other nine spectra for each of the ten and ask whether the true same-compound partner ranks first — **every** method scores top-1 = 1.00 (mean reciprocal rank 1.00)' (lines 164-166; SI Table S7 confirms top-1=1.00, MRR=1.00 for all six methods, verified in results/retrieval_13c.json). On the concrete chemical task the paper adds to answer my significance objection, LCC has zero measured advantage over the saturating tree/NN and the bin method. An honest null result is commendable, but it does not demonstrate significance — it documents its absence on the only real-task test in the paper.
  - *Suggested fix:* Replace the small clean-library retrieval with an at-scale test against a real decoy library (BMRB/HMDB) or under a rejection threshold, where separation actually converts to a retrieval difference — i.e. run the prospective dereplication study the paper names as 'the natural next step' (line 196), and report it as the result rather than as future work. Only such a demonstration would lift significance to Angewandte's level.
  - *Verifier:* The referee's factual core is confirmed: all six methods score top-1 = 1.00, MRR = 1.00 (results/retrieval_13c.json — every row "top1": 1.0, "mrr": 1.0; SI Table S7 lines 409-414; MS lines 164-166). But the finding's MAJOR framing overstates in three verifiable ways.

(1) "LCC has zero measured advantage over the saturating tree/NN" is false. In the SAME paragraph (MS lines 160-163) and Table 1 note (lines 239-240) the revision reports disjoint bootstrap CIs on sparse separation: LCC 0.78 (0.55-0.98), tree 0.41 (0.34-0.49), NN 0.10 (0.09-0.11) — confirmed in JSON ci95 fields (lcc [0.547,0.981], tree [0.338,0.488], nn [0.089,0.112]). Non-overlapping CIs = statistically significant advantage. LCC ties tree/NN only on the retrieval metric, not the discriminating one.

(2) The paper never claims retrieval demonstrates significance; it says the opposite. MS lines 166-168: "retrieval becomes discriminating only under a rejection threshold or at library scale — which is exactly what the separation and margin measure." SI lines 421-423: "Retrieval accuracy therefore does not discriminate the methods at this scale; the discriminating quantity is the thresholded same/different decision (separation and margin)." The significance claim rests on separation/margin + bootstrap, not retrieval.

(3) "Ties LCC with the methods it claims to beat" mischaracterizes the calibrated claims. MS lines 173-176 claim LCC "beats the saturating tree and nearest-neighbour methods... and it matches or exceeds the bin method," and explicitly downgrades the sparse bin comparison to "comparably... within the bootstrap CI" (line 175), with "LCC's edge over the bin method is within noise" (line 162-163). The methods it claims to BEAT are tree/NN, which it does beat on separation via disjoint CIs; it never claims to beat the bin method on the sparse set.

The honest residual (no significant edge over the bin method on sparse; retrieval ties everyone) is real but fully and repeatedly disclosed (MS 162-171, 239-243; SI 418-426), and the referee concedes the null is honestly reported. The revision already contains the caveat that settles the charge.
- **[MAJOR · REMAINING] On the reproducible benchmark the sole novel ingredient does nothing; LCC reduces to the textbook contrast-angle cosine** — *verification: OVERSTATED → MINOR*
  - *Location:* MS Score step (lines 99-110); MS lines 154-156; Table 1 (lines 232-233); SI S10, S11
  - The revision's own numbers show the un-centred cosine ablation equals LCC on the public sparse set: Table 1 row 'Cosine, un-centred (LCC ablation) | 0.59 | — | 0.78 | 0.37' vs 'LCC | 0.75 | 0.71 | 0.78 | 0.37', and the text confirms 'in this sparse regime the un-centred cosine matches LCC exactly (separation 0.78)' (lines 154-156). Mean-centring — the paper's one 'this work' ingredient (SI Table S6) and the analogue of 'the cosine / contrast-angle similarity long used in mass-spectral library search' (lines 101-102) — therefore adds value only on the dense benchmark, which the author concedes rests on a single negative spectrum. The claimed conceptual advance is thus demonstrable only where the evidence is weakest (n=1 negative).
  - *Suggested fix:* Either establish the dense-regime advantage on a population of unrelated proteins (multiple decoy folds/sizes with ROC/AUC), so the novel step is validated where it actually operates, or reframe the contribution honestly as an engineering synthesis for a methods venue rather than a conceptual advance for a general journal.
  - *Verifier:* Every factual sub-claim the referee makes is true and is drawn verbatim from the paper's own disclosures: Table 1 (lines 232-233) shows Cosine-un-centred sparse sep 0.78 / margin 0.37 identical to LCC 0.78 / 0.37; the text at 154-156 states "in this sparse regime the un-centred cosine matches LCC exactly (separation 0.78)... so its benefit is specific to the dense fingerprint"; the abstract (42-43) and lines 192-193 concede "its most decisive number, the dense-protein margin, rests on a single negative spectrum"; SI Table S6 attributes mean-centring alone to "this work"; and Benchmark 1 (S8) has exactly one different-protein spectrum (OAA exp 103). So the underlying observation — the one novel ingredient's demonstrable value is confined to the n=1-negative dense benchmark — is real and is already openly disclosed, not hidden. That makes it a known, conceded limitation, not a newly-exposed fatal flaw.

But the referee's MAJOR framing is OVERSTATED on two counts a competent author can rebut. (1) "LCC reduces to the textbook contrast-angle cosine" is false. The code shows cosine_similarity = lcc_similarity(center=False) (hsqc_lcc.py:158-159; _zncc line 103-104): the ablation is a cosine of the SAME shared-grid, area-weighted, physically-Gaussian-blurred images — the render+blur pipeline that constitutes the rest of LCC — not a raw textbook cosine. The paper states this explicitly ("it is an ablation, not an independent method," lines 235-236; S10). On the reproducible sparse benchmark that pipeline still delivers sep 0.78, decisively beating the actual textbook-lineage methods (tree 0.41, NN 0.10; Table 1 / Table S3) and matching bin. So "the sole novel ingredient does nothing" conflates mean-centring with LCC's whole reproducible advantage; LCC as a measure is validated on public data. (2) The paper's claimed advance is the unified measure that wins BOTH regimes, with mean-centring the discriminating step specifically in the crowded regime "that is the hard case this measure targets" (S10) — the single-negative issue is a data-availability limit (public deposition "in preparation," lines 269-270), which the sparse benchmark + retrieval test are explicitly added to offset (S13). Because the paper neither overclaims nor hides anything, and the referee mischaracterizes the ablation and ignores LCC's real reproducible win over the textbook methods, this is a disclosed evidence-strength limitation (genuine, MINOR), not an unrebuttable MAJOR defect.
- **[MINOR · NEW] Abstract's 'widest separation in both regimes' is undercut by the paper's own new bootstrap**
  - *Location:* MS abstract (lines 36-40); MS lines 173-176; SI S13 (lines 392-401)
  - The abstract asserts 'LCC gives the widest same/different separation in both regimes (0.75 dense, 0.78 sparse)' (lines 36-37), yet two sentences later concedes 'A bootstrap shows LCC's sparse-regime edge over the bin method is within noise at five compounds' (lines 39-40) — LCC 0.78 (CI 0.55-0.98) vs bin 0.65 (CI 0.42-0.84), overlapping (SI Table S7). So on the sparse regime 'widest' is a point estimate the paper's own statistics cannot resolve from the incumbent. Disclosed, but the lead claim still reads stronger than the data support.
  - *Suggested fix:* Qualify the sparse-regime lead in the headline sentence (e.g. 'ties the bin method within bootstrap noise on the sparse set, and decisively leads on the dense fingerprint'), so the abstract and the CI tell the same story.
- **[MINOR · NEW] Inconsistent ablation-margin reporting in Table 1**
  - *Location:* Table 1 (lines 232-237)
  - Table 1 gives the un-centred cosine a dense-protein margin of '—' but a sparse margin of 0.37 (line 232), and the note explains 'Its dense-protein margin was not recomputed here (the ablation reports separation only)' (lines 235-237). Reporting the ablation margin in one regime but not the other is an avoidable asymmetry in the table that most directly isolates the paper's novelty claim.
  - *Suggested fix:* Compute and report the un-centred-cosine dense margin (bench.py already produces per-comparison scores), so the ablation row is complete in both regimes.

### Referee 2 — MAJOR REVISION

*Focus: Correctness of the new bootstrap CIs and cosine ablation; whether the honest sparse-regime disclosure and retrieval null are computed and reported correctly; whether the n=1 dense negative still undercuts the headline margin 0.71 and the central mean-centring novelty claim.*

**Assessment.** The revision is a genuine and largely commendable improvement in statistical honesty against my first-round asks. The bootstrap CIs are computed correctly and reproduce exactly (I re-ran the sparse pipeline: LCC sep 0.7778, CI [0.5467,0.9812]; bin 0.6457, CI [0.4183,0.8374]), match results/retrieval_13c.json and the manuscript, and the paper now openly states LCC's sparse edge over the bin method is within noise — the direction of caution is correct (overlapping 95% CIs are used only to decline a claim, not to assert one). The olivetol byte-identical duplicate is verifiably removed (same MD5), 26 tests pass, and the retrieval top-1=1.00 null result is reported rather than hidden. However, two load-bearing issues remain that my first round flagged and the revision did not close. First, the paper's single novel ingredient — mean-centring — is now shown by the authors' own honest ablation to do NOTHING in the only regime with a real negative set (sparse, 40 negatives: centred 0.778 vs un-centred 0.777); its entire demonstrated benefit (0.59→0.75) exists solely against the single dense OAA negative, so the central novelty rests on n=1, and the abstract now leans harder on this than before. Second, that dense 0.59 ablation value is a hard-coded literal (comparison_all.json, make_plots.py) that is absent from the primary dense results JSON and cannot be regenerated without the still-non-public PRL3/OAA data. The math is sound and the sparse story is fully reproducible, but the paper's most decisive numbers — margin 0.71 and the 0.75→0.59 ablation — both still hang on one non-public negative spectrum. Major revision from my lens.

**Prior concerns:** RESOLVED: (a) missing/non-reproducible cosine ablation — now reproducible in the SPARSE regime (0.777 in comparison_13c.json, I reproduced it); (b) uncertainty quantification — bootstrap 95% CIs added, correctly computed, and honestly interpreted; (c) olivetol duplicate removed and guarded by a raising check + regression test; (d) 'faithful' softened to 'representative reimplementations'; (e) self-similarity=1 no longer sold as an LCC advantage ('as for every method tested'). PARTIALLY_RESOLVED: the n=1 dense negative is now thoroughly DISCLOSED (abstract, body, conclusion) but not REMEDIED — no additional decoy proteins were added, so margin 0.71 and the 0.59 ablation still rest on one spectrum, and the bare 'three times the bin method' multiplier survives in the dense regime that still has no CI. REMAINS: dense cosine 0.59 not in the primary results JSON / not publicly reproducible; 3x vs 2.6x basis inconsistency; grid step still framed away as non-tuning; 'strictly decreasing' monotonicity only scoped, not fully qualified.

**Findings:**

- **[MAJOR · NEW] Central novelty (mean-centring) has zero demonstrated benefit in the only well-powered regime; its entire measured value rests on n=1** — *verification: OVERSTATED → MINOR*
  - *Location:* MS abstract lines 39-43; MS lines 134-135, 154-156; SI S10 lines 307-314
  - Sparse regime (40 negatives): 'the un-centred cosine matches LCC exactly (separation 0.78)' (MS line 154-156); comparison_13c.json gives lcc 0.7778 vs cosine 0.7773. Dense regime (1 negative): 'replacing it with an un-centred cosine ... drops the separation from 0.75 to 0.59' (MS line 134-135). So mean-centring's only demonstrated benefit is 0.59->0.75 against the single OAA spectrum.
  - *Suggested fix:* State explicitly that the mean-centring benefit is demonstrated only against one dense negative and vanishes on the 40-negative sparse set, so the novelty claim is not yet statistically powered. Add several decoy proteins to the dense benchmark so the 0.59->0.75 ablation is supported by a negative distribution, or reframe mean-centring as a mechanism argued (S6b) rather than empirically established.
  - *Verifier:* The referee's three factual claims are all TRUE and confirmed in code/SI: (1) sparse regime, cosine ≈ LCC — comparison_13c.json gives lcc_new.separation 0.7778 vs cosine_uncentred.separation 0.7773 (MS line 154-156, SI S10 "separation 0.777 vs 0.778"); (2) dense regime, mean-centring moves separation 0.59→0.75 (MS 134-135, SI S10); (3) the dense benchmark has a single negative spectrum (OAA). So the demonstrated benefit of the central novelty (mean-centring) does rest on n=1 negative.

BUT the referee frames this as a hidden/fatal flaw ("zero demonstrated benefit," "entire measured value rests on n=1") that the revision misses — and it does not miss it. The revised manuscript already discloses this openly and repeatedly, and bounds the claim to exactly what the data support:
- Abstract (41-43): "removing the mean-centring step collapses the dense-protein separation from 0.75 to 0.59, yet in the sparse regime — where peaks rarely overlap — centred and un-centred scores coincide (0.78), so mean-centring is the operation that rescues discrimination specifically when peaks crowd." The claim is deliberately dense-specific, not universal.
- Line 154-156 states the sparse coincidence itself and explains it: "where crosspeaks rarely overlap there is no coincidental non-co-located intensity for mean-centring to penalize, so its benefit is specific to the dense fingerprint."
- Line 192-193: "its most decisive number, the dense-protein margin, rests on a single negative spectrum, which is why we add the larger public sparse benchmark and the retrieval test." — explicit n=1 disclosure.
- SI S10 gives the mechanistic derivation (§S6b) predicting the effect is dense-specific, so the sparse null is theoretically expected, not a failure.

A competent author CAN rebut the "defect" framing: they do not claim a well-powered demonstration of mean-centring across many negatives; the sparse null is a mechanistic prediction they confirm, and the n=1 dense limitation is stated three times. What the author cannot do is add statistical power — the substantive concern (central novelty's empirical benefit validated against only one negative spectrum) is real and legitimate as a request for more dense negatives. That residual is a disclosed, mechanistically-bounded scope limitation, not a correctable defect or an overclaim. The referee's factual core survives but its severity/spin ("zero benefit / entire value rests on n=1" as if concealed) is overstated relative to a revision whose claims are calibrated to precisely this evidence.
- **[MAJOR · REMAINING] Dense cosine ablation (0.59) is a hard-coded literal, absent from the primary results JSON, and not publicly reproducible** — *verification: OVERSTATED → MINOR*
  - *Location:* results/method_comparison.json (no cosine row); results/comparison_all.json line 26-31; results/make_plots.py lines 27-33; MS abstract lines 39-43
  - bench.py now computes cosine_uncentred, but results/method_comparison.json has no cosine row (only bin/bin_rot45/tree/nn/lcc_new); the 0.59 lives only as a literal in comparison_all.json ('sep15': 0.59) and make_plots.py line 33, whose own comment says 'not regenerable here without the private Bruker data'. The abstract headlines this exact number (0.75->0.59).
  - *Suggested fix:* Regenerate method_comparison.json with the cosine_uncentred row from an actual bench.py run and deposit the dense PRL3/OAA data (or the windowed rendered matrices) so both 0.75 and 0.59 are reproducible from a public artefact. Until then, mark 0.59 in Table 1 as a value that cannot be independently regenerated.
  - *Verifier:* The referee's facts check out: results/method_comparison.json's separation block lists only bin_Bodis09, bin_rot45, tree_Castillo13, nn_Pierens12, lcc_new — no cosine row (file dated Jul 8, stale, predating the cosine addition). The dense 0.59 appears only as a literal in comparison_all.json ("cosine_uncentred":{"sep15":0.59}) and make_plots.py PROTEIN dict ("cosine_uncentred":(0.59, None)). make_plots.py docstring: "method_comparison.json is not regenerable here without the private Bruker data." Abstract headlines "0.75 to 0.59." All accurate.

But the MAJOR framing is overstated. (a) The dense 0.59 is not uniquely non-reproducible — the entire dense 1H-15N benchmark requires private PRL3/OAA data; the headline LCC 0.75 (sep15 in method_comparison.json) is from the same private run, as are bin 0.29 etc. The make_plots.py comment the referee quotes applies to ALL dense numbers, not cosine specifically. (b) The revision openly discloses this in the Data Availability Statement (MS 262-270): "The dense 1H-15N benchmark uses PRL3 and OAA Bruker spectra available from the author on reasonable request... a public deposition of the processed spectra is in preparation so that the dense benchmark also becomes independently reproducible." (c) bench.py:39 now computes cosine_uncentred via the same code path as every other method, and the SPARSE cosine ablation (0.777) is public, live and reproducible (bench_13c.py -> comparison_13c.json / Table S3), confirming the ablation mechanism. SI 314-315: "a real, reproducible benchmark row in both harnesses (cosine_similarity in hsqc_lcc.py), not an asserted number" — accurate about the code.

Residual genuine but MINOR issue: method_comparison.json is stale, so within the author's private workflow the dense 0.59 in the repo remains a hand-transcribed figure literal not read back from a regenerated harness output; the SI's "reproducible benchmark row in both harnesses" is slightly overstated w.r.t. the committed output artifacts (the dense output JSON lacks the row), though accurate about the harness code. Fixable artifact-hygiene gap, not fabrication or a reproducibility failure beyond the already-disclosed private-data limitation.
- **[MAJOR · PARTIALLY_RESOLVED] Dense benchmark still n=1 negative with no CI, yet quoted as 'three times' and 'decisive'** — *verification: OVERSTATED → MINOR*
  - *Location:* MS line 133; abstract line 41; MS lines 168-171, 191-194; bench.py diff = read_bruker_2d(oaa/'103')
  - 'a margin of 0.71 - three times the bin method's 0.24' (MS line 133) and abstract 'its dense-fingerprint advantage is decisive' (line 41), while the honest caveat that this 'rests on a single negative spectrum' appears at lines 168-171/191-194. The bootstrap was added only for the sparse regime; the dense n=1 admits no CI, so a bare multiplier and 'decisive' overstate precision the data cannot support.
  - *Suggested fix:* Drop or qualify the bare 'three times' and 'decisive' for the dense regime until backed by a distribution over multiple decoy proteins; state 'three times on worst-case margin against a single negative'.
  - *Verifier:* The referee's facts are accurate but the finding is overstated because the revision already contains prominent, repeated caveats the referee itself acknowledges. The "three times" at MS line 133 ("a margin of 0.71 — three times the bin method's 0.24") is a literal arithmetic ratio (0.71/0.24≈2.96), stated as fact with no implied CI. "Decisive" is an effect-size/ranking claim, not a precision claim, and is caveated adjacent to every use: line 169-171 "The dense-protein margin (0.71) remains LCC's strongest result, but it rests on a single negative spectrum; the sparse regime and this retrieval test broaden the evidence with a larger, public, fully-reproducible negative set"; line 191-194 "its most decisive number, the dense-protein margin, rests on a single negative spectrum, which is why we add the larger public sparse benchmark and the retrieval test"; and SI S13 line 398-401 "LCC's decisive, non-overlapping advantage is in the dense protein regime (margin 0.71), which unfortunately rests on a single negative spectrum (§S8) — the limitation the sparse benchmark and this retrieval test are added to offset." The abstract itself (line 39-41) draws the correct distinction — "within noise at five compounds, while its dense-fingerprint advantage is decisive" — where "decisive" means decisively above the bin baseline (0.71 vs 0.24, negative 0.18 below every same-protein score 0.89-0.97), versus sparse where LCC ties the bin method within the bootstrap CI. The manuscript explicitly frames ranking, not exact numbers, as the result (line 183-184). The only residual is a minor wording matter: the abstract's "decisive" is not re-qualified inside the abstract, so juxtaposition with "within noise" could momentarily read as a statistical rather than effect-size claim. That is a softening suggestion, not an unrebuttable MAJOR defect.
- **[MINOR · REMAINING] Bootstrap resamples 40 different-compound pairs as if independent (pseudoreplication)**
  - *Location:* bench_retrieval.py bootstrap_sep (diff resample); SI S13 lines 392-401
  - bench_retrieval.py bootstrap_sep resamples diff_vals (40 pairs) independently, but the 40 pairs come from 10 spectra each appearing in 9 pairs, so they are correlated; the percentile CI on the diff side is anti-conservative. In practice immaterial: mean_diff=0.003 with near-zero variance, so CI width is dominated by the 5 same-compound values (genuinely independent compounds). Conclusion ('within noise vs bin') is unaffected and, if anything, conservative.
  - *Suggested fix:* Add one sentence noting the different-pair non-independence and that the CI is dominated by the 5 same-compound draws; optionally report a compound-clustered bootstrap for completeness.
- **[MINOR · REMAINING] Inconsistent superiority multiplier (3x vs 2.6x) with basis not flagged**
  - *Location:* MS line 133 vs SI line 252
  - MS line 133 'three times the bin method's 0.24' (margin, 0.71/0.24=2.96) vs SI S8 line 252 'approx 2.6x better than the previous best' (separation, 0.75/0.29=2.59). Both correct on their basis, neither states which.
  - *Suggested fix:* State the basis at each mention ('3x on worst-case margin', '2.6x on separation') or unify to one figure.
- **[MINOR · REMAINING] Grid step still framed as non-tuning; 'grid-invariant' vs measured 0.36 when coarsened**
  - *Location:* MS abstract line 44; SI Table S1 line 217; SI S10 line 320; bench_13c.py line 53
  - Abstract line 44 'both knobs set by spectroscopy rather than by tuning'; Table S1 line 217 labels the render step 'score grid-invariant'; yet the step changes between regimes (bench_13c.py line 53: 0.02/0.2 vs default 0.01/0.10) and S10 line 320 reports separation falls to 0.36 when coarsened, so invariance is only approximate and the step is a third, non-spectroscopic knob.
  - *Suggested fix:* Either fix the step by a stated rule (e.g. sigma/3) or acknowledge the grid step as a third parameter and report the actual grid dependence rather than 'grid-invariant'.

### Referee 3 — MAJOR REVISION

*Focus: Whether bench.py no longer hard-codes private paths; whether the sparse + retrieval benchmarks are fully reproducible from public data; soundness of the olivetol duplicate guard and the honesty of removing it; adequacy of the "deposition in preparation" promise for the flagship dense benchmark; and any new over/under-claiming introduced by the revision.*

**Assessment.** This is a substantially improved, unusually honest revision, and most of my first-round reproducibility/fairness concerns are genuinely fixed. I verified in the actual code that bench.py no longer contains the private absolute paths I flagged — it now reads $PRL3_DIR/$OAA_DIR (or --prl3/--oaa) and errors cleanly when they are absent. I ran the public benchmarks: the sparse 1H-13C harness auto-downloads the simpleNMR peak lists and reproduces Table S3/Table 1 to the digit (LCC sep 0.7778, cosine 0.7773, bin 0.6457, tree 0.4137, NN 0.1009), bench_retrieval reproduces Table S7's CIs and the top-1=1.00 null result, and 26/26 tests pass. The olivetol removal is sound and honest: the _assert_distinct guard genuinely rejects a same-pair that renders identically, and removing the byte-identical positive correctly LOWERED the sparse separation (0.81->0.78) rather than flattering it. The un-centred cosine ablation is now a reproducible row, and "faithful" was softened to "representative." What still stands between this and acceptance are two things under my lens. First, my one CONFIRMED first-round MAJOR is only PARTIALLY resolved: the flagship dense benchmark — which the paper itself now calls "its most decisive number" — remains non-reproducible from any public artifact, and "a public deposition is in preparation" (no DOI, no date, no embargo reason) is not an adequate substitute for a demanding general journal. Second, the revision introduces a new abstract framing that promotes exactly that non-reproducible, un-bootstrappable n=1 result as "decisive" while downgrading the actually-quantified sparse result to "within noise" — statistically inverted, and the n=1 caveat that appears in the body is missing from the abstract. The remaining items are cosmetic round-one minors left untouched.

**Prior concerns:** RESOLVED: (1) hard-coded private paths in bench.py — verified gone, now $PRL3_DIR/$OAA_DIR with a clear error; (2) missing/asserted cosine ablation — now a reproducible benchmark row matching the JSON (0.777); (3) "faithful reimplementations" — softened to "representative" (L181); (4) olivetol byte-identical duplicate — removed, guarded (_assert_distinct), regression-tested, and disclosed, with the score moving in the honest direction. NEWLY ADDED and verified reproducible: bootstrap 95% CIs and the dereplication/retrieval test on public data. PARTIALLY RESOLVED: my CONFIRMED MAJOR (dense benchmark not independently reproducible) — the path half is fixed and the paper leans on the fully-public sparse+retrieval legs, but the dense flagship data are still request-only with only an "in preparation" deposition promise. REMAIN (round-1 minors, untouched): 2.6x-vs-3x multiplier basis; method-name hyphenation; grid step still changed per regime and still labelled "grid-invariant"/"set by spectroscopy" while S10 shows the score does move (0.36 at coarse grid).

**Findings:**

- **[MAJOR · PARTIALLY_RESOLVED] Flagship dense benchmark still not independently reproducible; 'deposition in preparation' is not adequate** — *verification: OVERSTATED → MINOR*
  - *Location:* LCC_angewandte.md L268-270 (Data Availability); bench.py run(); results/method_comparison.json
  - Data Availability (L268-270): dense spectra 'available from the author on reasonable request ... a public deposition of the processed spectra is in preparation so that the dense benchmark also becomes independently reproducible.' No DOI, no date, no embargo reason. The stored results/method_comparison.json (lcc_new 0.7549, bin 0.2924) cannot be regenerated by a referee. The paper itself calls this 'its most decisive number, the dense-protein margin' (L192-193).
  - *Suggested fix:* Actually deposit the processed/windowed PRL3 and OAA matrices (Zenodo/BMRB) and cite the DOI in-text before acceptance. A promise 'in preparation' does not satisfy the data-availability bar for the datum the paper labels most decisive; the fix is trivial now that bench.py reads paths from env/flags.
  - *Verifier:* Referee's factual claim is true: dense PRL3/OAA Bruker data are not public. bench.py L77-80 errors without --prl3/--oaa; the stored margin in results/method_comparison.json (lcc_new 0.7549 vs bin 0.2924) cannot be regenerated by a referee from public data. Data Availability (L263-266) confirms only "available from the author on reasonable request ... a public deposition ... is in preparation" — no DOI, date, or embargo reason in the manuscript. So the gap exists (not REFUTED).

But the MAJOR severity is overstated because the revision pre-emptively rebuts the "most decisive number" framing: (1) L192-193 explicitly de-weights it — "its most decisive number, the dense-protein margin, rests on a single negative spectrum, which is why we add the larger public sparse benchmark and the retrieval test"; the paper deliberately moves the load-bearing evidence off the irreproducible number. (2) L263-266: "the sparse 1H–13C benchmark and the dereplication/retrieval experiment are fully reproducible from public data: the harnesses download the simpleNMR peak lists ... automatically and contain no hard-coded paths" — verified in bench.py L68-69 and the parametrized --prl3/--oaa / $PRL3_DIR/$OAA_DIR interface. (3) A legitimate embargo reason exists in bench.py L67: "the dense PRL3/OAA Bruker data are not redistributable" — the withheld item is proprietary raw spectra, not analysis code (all open MIT). The author can therefore rebut MAJOR: the irreproducibility is acknowledged, justified, and no longer load-bearing. Residual weakness (soft "in preparation" pledge, no DOI/date) is minor reproducibility hygiene.
- **[MAJOR · NEW] Abstract promotes the non-reproducible n=1 dense result as 'decisive' while downgrading the quantified sparse result to 'within noise' — statistically inverted** — *verification: OVERSTATED → MINOR*
  - *Location:* LCC_angewandte.md L39-41 (abstract) vs L169, L192-193 (body)
  - Abstract L39-41: 'A bootstrap shows LCC's sparse-regime edge over the bin method is within noise at five compounds, while its dense-fingerprint advantage is decisive.' The dense advantage rests on a single negative spectrum (bench.py diff = one OAA spectrum) and therefore cannot be bootstrapped at all, whereas the sparse result is the one with an actual computed CI. The n=1 caveat is stated in the body (L169, L192-193) but omitted from the abstract, so the abstract elevates the least-reproducible, least-quantifiable number.
  - *Suggested fix:* In the abstract, attach the n=1 caveat to the word 'decisive' (or replace it): e.g. 'its dense-fingerprint margin is large but rests on a single negative spectrum.' Do not let a bootstrapped 'within noise' sit rhetorically below an un-bootstrappable 'decisive.'
  - *Verifier:* The referee's factual observations are all TRUE: (1) the dense different-class is a single spectrum — bench.py:48 `diff = read_bruker_2d(oaa / "103")`, with margin computed as `min_same - diff_s` against that one negative, so it cannot be bootstrapped; (2) the sparse result carries the actual computed CI (abstract L39-40 "bootstrap ... 95% CI 0.55–0.98"; body L160-162); (3) the abstract does not contain the n=1 caveat, while the body discloses it twice — L169 "The dense-protein margin (0.71) remains LCC's strongest result, but it rests on a single negative spectrum" and L192-193 "its most decisive number, the dense-protein margin, rests on a single negative spectrum, which is why we add the larger public sparse benchmark and the retrieval test." So the referee did NOT miss a fix: the abstract genuinely omits the caveat.

But the "statistically inverted" framing is OVERSTATED. The abstract does not claim the dense advantage was bootstrapped or is statistically significant. The two words describe two correctly-characterized, different quantities: "within noise" reports the bootstrap CI overlap for sparse (LCC 0.78 CI 0.55–0.98 vs bin 0.65 CI 0.42–0.84), and "decisive" describes the dense effect-size gap (margin 0.71 vs 0.24, per body L175 "decisively on the dense amide fingerprint (margin 0.71 vs 0.24)" and Table 1) — a real, computed ~3x separation, not a significance claim. "A bootstrap shows [sparse within noise], while [dense advantage is decisive]" reads as a contrast; the "decisive" clause is independently supported by the margin, not attributed to a bootstrap. The dense advantage being decisive in effect-size terms is TRUE on its own merits. The only genuine residue is an editorial judgment — whether the abstract should also carry the n=1 reproducibility caveat that the body states twice. That is a MINOR presentation/emphasis point (abstracts routinely defer caveats to the body), not a factual error or a MAJOR statistical inversion a competent author could not rebut.
- **[MINOR · REMAINING] Two different multipliers for the same headline comparison, neither stating its basis**
  - *Location:* LCC_angewandte.md L133; LCC_angewandte_SI.md L252
  - MS L133: 'a margin of 0.71 — three times the bin method's 0.24' (margin basis). SI L252: 'LCC separates same-protein from different-protein spectra ~2.6x better than the previous best' (separation basis, 0.75/0.29). Both correct on their own basis but a cross-referencing reader sees 3x and 2.6x for what looks like one claim. Flagged in round 1, untouched.
  - *Suggested fix:* State the basis at each mention ('3x on worst-case margin', '2.6x on separation') or unify to one multiplier.
- **[MINOR · NEW] Bootstrap treats the 40 different-compound pairs as independent when they are C(10,2) combinations of 10 spectra**
  - *Location:* bench_retrieval.py bootstrap_sep(); LCC_angewandte_SI.md S13 L392-393
  - bench_retrieval.bootstrap_sep resamples diff_vals (40 values) i.i.d., and S13 says '(independently resampling the 5 same-compound and 40 different-compound pairs)'. Each of the 10 spectra recurs in up to 8 different-pairs, so the pairs are not independent and the naive resample understates variance. This is conservative for the paper's stated 'CIs overlap / within noise' conclusion (a wider true CI only reinforces it), so it does not overturn any claim, but the non-independence is undisclosed.
  - *Suggested fix:* Add one sentence noting the different-pairs are non-independent (pair-level resample understates variance, so the reported CIs are if anything optimistically narrow), or bootstrap at the spectrum/compound level.
- **[MINOR · NEW] 'Widest separation in both regimes' is a numerical tie with LCC's own ablation in the sparse regime**
  - *Location:* LCC_angewandte.md L36-37 vs L39-40; Table 1 (L232-233)
  - Abstract L36-37: 'LCC gives the widest same/different separation in both regimes (0.75 dense, 0.78 sparse)'. In the sparse regime LCC 0.778 barely edges its own un-centred-cosine ablation 0.777 (Table 1 / comparison_13c.json) and its CI overlaps the bin method's. So the sparse 'lead' is within its own ablation and within the bin CI — mild tension with the abstract's own 'within noise' statement two sentences later.
  - *Suggested fix:* In the sparse regime, phrase LCC as matching the bin method and its own centred/un-centred forms (i.e. lead is over tree/NN, not over bin/cosine), consistent with the honest 'within noise' sentence.
- **[MINOR · REMAINING] 'Both knobs set by spectroscopy' still undercut by the per-regime grid-step change labelled 'grid-invariant'**
  - *Location:* LCC_angewandte.md L44; LCC_angewandte_SI.md Table S1 L217, S10 L320-321; bench_13c.py METHODS
  - Abstract L44 and Table S1 (L217, '$\geq$3 px per $\sigma$; score grid-invariant') claim spectroscopic, non-tuned knobs, yet the render step is changed 0.01/0.10 -> 0.02/0.2 between regimes (bench_13c METHODS) and S10 (L320-321) shows the score does move with grid (0.36 at bin resolution). Round-1 minor, untouched.
  - *Suggested fix:* Either tie the step to sigma by a stated rule (e.g. sigma/3) so it is not a free knob, or acknowledge the grid step as a third, non-spectroscopic parameter and state the actual (approximate) grid dependence.

### Referee 4 — REJECT

*Focus: Internal numerical consistency across abstract/Table 1/Figure 2/SI Tables S3-S7; clarity of the new cosine and retrieval rows; format compliance (references, captions, keywords, TOC); residual overclaims/inconsistencies.*

**Assessment.** From my presentation/format lens the revision is a clear, substantial improvement: I verified that the sparse separation is 0.78 everywhere (JSON lcc_new = 0.7778), the benchmark is 5 pairs / 40 different in every table, the bootstrap CIs in Table 1, the main text and SI Table S7 all match results/retrieval_13c.json (LCC 0.55-0.98, bin 0.42-0.84, tree 0.34-0.49, NN 0.09-0.11), the test count 26 is consistent and confirmed by pytest, olivetol appears only as an excluded duplicate with no dangling positive, and the references are 27 in correct Angewandte style. Almost my entire first-round presentation cluster is resolved (ZNCC removed, \"faithful\"->\"representative\", ablation now a consistent 0.16 with no stray 0.17, self-similarity no longer framed as distinguishing LCC, keywords=5, self-contained captions). What remains is minor polish plus two structural items outside presentation: the dense benchmark is still not independently reproducible (data on request), and — decisively for Angewandte — the paper still demonstrates no chemical result. My recommendation is therefore reject for Angewandte's significance bar, but I want the editor to note that on presentation grounds this is now a clean minor-revision-quality manuscript that would transfer directly to a specialized methods venue.

**Prior concerns:** RESOLVED from my first round: ZNCC undefined-term removed from all captions; \"faithful reimplementations\" softened to \"representative\" (MS L181, SI S10); the mean-centring ablation figure is now a consistent +0.16 (SI S6b L202 = 0.75-0.59) with no residual 0.17; self-similarity is no longer framed as LCC-distinguishing (MS L127 \"Every method self-scores exactly 1.00\", L189); \"strictly decreasing\" is now bounded by an explicit idealization (SI S5 L162-163). PARTIALLY: bare multipliers — CIs added but 3x/2.6x juxtaposition persists; lead-metric alternation — abstract now foregrounds separation but dense paragraph still leads with margin. REMAIN: method-name hyphenation title-vs-body; TOC graphic still \"suggested\" not supplied as a file (acceptable at text-submission stage, not separately flagged); grid step still not disclosed as a third parameter (carryover finding #7).

**Findings:**

- **[MAJOR · REMAINING] Reject-driving significance gate unaddressed: no chemical result** — *verification: OVERSTATED → MINOR*
  - *Location:* MS abstract; MS L191-198; SI S13
  - MS L191-193: "this is a **methodological** contribution — a robust, reproducible similarity primitive, not a new chemical discovery — and its most decisive number, the dense-protein margin, rests on a single negative spectrum". The revision honestly re-scopes itself as methodological but adds no binder identified, no CSP-mapped site, no dereplicated unknown; the new retrieval test is an explicit null result (SI S13 L418: "every method achieves top-1 = 1.00").
  - *Suggested fix:* For Angewandte, add one real chemical result where LCC changes a conclusion an incumbent gets wrong (a blind BMRB/HMDB dereplication hit, or a CSP/titration ranking where bins/tree/NN fail). Otherwise the honest home is a specialized methods journal, to which the manuscript transfers directly.
  - *Verifier:* The referee's factual descriptions are all correct — but they describe caveats the revision itself volunteers, not a defect the author cannot rebut. The manuscript pre-emptively and explicitly rebuts the exact charge. MS L191-198: "We are candid about scope: this is a **methodological** contribution — a robust, reproducible similarity primitive, not a new chemical discovery — and its most decisive number, the dense-protein margin, rests on a single negative spectrum, which is why we add the larger public sparse benchmark and the retrieval test." The abstract (L30-36) frames LCC purely as a similarity method ("Here we introduce the Lineshape Correlation Coefficient..."), never claiming a chemical discovery, so there is no overclaim to gate against. The retrieval null result is not buried but reported with full candor — SI S13 L418: "This is an honest null result... We report the null result rather than omit an experiment that did not favour the proposed method" — and its non-discriminating nature is explained ("on a small, clean library the correct partner is always the single most-similar spectrum"), with a stated path to a discriminating study (BMRB/HMDB at-scale dereplication). The referee is therefore not identifying a factual error or an unaddressed weakness; the referee is asserting that a methodological contribution lacking a novel chemical result should be rejected. That is a significance/scope judgment reserved to the editor, and it is one the revision confronts head-on and legitimately answers. A competent author can and does rebut it: the paper is honestly and consistently positioned as a method note contributing a similarity primitive, a recognized category. The "reject-driving MAJOR" characterization overstates an acknowledged, disclosed limitation as an incurable defect.
- **[MAJOR · PARTIALLY_RESOLVED] Flagship dense benchmark still not independently reproducible** — *verification: OVERSTATED → MINOR*
  - *Location:* MS Data Availability L260-270; bench.py L45-63; results/method_comparison.json
  - Hardcoded paths are fixed (bench.py L71-74 now reads --prl3/--oaa or $PRL3_DIR/$OAA_DIR) and the data-availability statement is updated. But the flagship 0.75/0.71 numbers remain non-regenerable from any public artefact: DAS L267-270 "PRL3 and OAA Bruker spectra available from the author on reasonable request ... a public deposition ... is in preparation". The committed results/method_comparison.json separation block contains only bin/bin_rot45/tree/nn/lcc_new — no negative-set distribution and no cosine row.
  - *Suggested fix:* Deposit the processed PRL3/OAA spectra (or the windowed rendered matrices) with a DOI so the dense benchmark reruns from public data; until then the paper's single most striking number cannot be verified by a referee.
  - *Verifier:* Referee's factual core is partly true, partly wrong. TRUE: bench.py L71-74 reads --prl3/--oaa or $PRL3_DIR/$OAA_DIR (no hardcoded paths); and the dense 0.75/0.71 numbers are not regenerable from any public artefact because the PRL3/OAA Bruker data are not deposited — DAS L267-270: "PRL3 and OAA Bruker spectra available from the author on reasonable request … a public deposition … is in preparation." FALSE: the referee's claim that method_comparison.json "contains only bin/bin_rot45/tree/nn/lcc_new — no negative-set distribution" is wrong; the committed json has 8 rows including `2_base vs OAA_103_diffprotein`, the negative (different-protein) score for all five methods. The dense benchmark has a single negative by design, which the manuscript states (L169 "it rests on a single negative spectrum"). The referee is correct only that the cosine_uncentred row (the 0.75→0.59 ablation) is absent from the committed json. Crucially the revision responded to this exact concern rather than leaving it: DAS is fully transparent about data-on-request + deposition-in-preparation, and the manuscript demotes the dense result and adds a public reproducible regime — L169-171: "The dense-protein margin (0.71) remains LCC's strongest result, but it rests on a single negative spectrum; the sparse regime and this retrieval test broaden the evidence with a larger, public, fully-reproducible negative set." The sparse 1H-13C benchmark (bench_13c.py) auto-downloads public simpleNMR peak lists (DAS L264-267) and is fully reproducible. So the reproducibility gap for the dense headline is real and unrebutted, but the referee mischaracterizes the artefact (negative row present) and overstates severity: transparent disclosure plus a corroborating public benchmark plus "available on request + deposition in preparation" (an accepted norm) reduce this to a MINOR deposition/documentation item, not a MAJOR blocker.
- **[MINOR · NEW] Dense un-centred-cosine bar (0.59) is a hardcoded constant, contradicting the "reproducible in both harnesses" claim**
  - *Location:* SI S10 L314; results/make_plots.py L26-33; MS Fig. 2 caption L217
  - SI S10 L314 asserts the cosine baseline is "a real, reproducible benchmark row in both harnesses (cosine_similarity in hsqc_lcc.py), not an asserted number." But for the DENSE regime the 0.59 is not in any committed data file; it is hardcoded in results/make_plots.py L33 `"cosine_uncentred": (0.59, None)` with the comment (L26-33) that method_comparison.json "is not regenerable here" and 0.59 is "the SI ablation value." So Figure 2's dense cosine bar is an asserted constant, reproducible only for the sparse harness from public data.
  - *Suggested fix:* Either regenerate and commit method_comparison.json with the cosine_uncentred row (bench.py already computes it, L36-39), or soften SI S10 to state the dense cosine value is the ablation figure and only the sparse cosine row is reproducible from public data.
- **[MINOR · REMAINING] Method-name hyphenation still inconsistent (title vs body)**
  - *Location:* MS L2 vs L30, L77, L204; SI L2
  - Title (MS L2) and SI title use "Lineshape-Correlation Coefficient" (hyphen); every in-text mention drops the hyphen — MS L30/L77/L204 "Lineshape Correlation Coefficient (LCC)". First-round finding #8 (presentation cluster) is unresolved on this point.
  - *Suggested fix:* Pick one form (Angewandte would take the unhyphenated "Lineshape Correlation Coefficient") and use it in title, running text and SI uniformly.
- **[MINOR · PARTIALLY_RESOLVED] Bare superiority multipliers still juxtapose 3x (main) and 2.6x (SI) for the same LCC-vs-bin dense comparison**
  - *Location:* MS L133-134; SI S8 L252
  - MS L133-134 "a margin of 0.71 — three times the bin method's 0.24" (margin ratio) while SI S8 L252 "LCC separates ... $\approx 2.6\times$ better than the previous best" (separation ratio 0.75/0.29). Each is correct for its own metric, and bootstrap CIs are now added, but the un-labelled 3x-vs-2.6x juxtaposition still reads as an internal inconsistency to a reader who does not track which metric each refers to.
  - *Suggested fix:* State the metric next to each multiplier ("3x on margin", "2.6x on separation") or drop the bare multipliers in favour of the now-available CI language.
- **[MINOR · REMAINING] "Two-knob, both set by spectroscopy" still omits the grid step as a de facto third parameter**
  - *Location:* MS L43-44, L176-178; SI Table S1 L217; SI S10 L320-321
  - MS L43-44 "a two-knob ... measure with both knobs set by spectroscopy rather than by tuning" and L176-178 repeat this. The render grid step (SI Table S1, `--step-f1/--step-f2`) is a non-spectroscopic choice, and grid-invariance is only approximate — SI S10 L320-321 shows separation drops to 0.36 when the grid is coarsened to the bin method's resolution. Carryover of first-round finding #7 (CONFIRMED MINOR).
  - *Suggested fix:* Either fix the grid step by a stated rule (e.g. step = sigma/3) so it is not free, or acknowledge it as a third, non-spectroscopic parameter and cite the actual grid dependence.
