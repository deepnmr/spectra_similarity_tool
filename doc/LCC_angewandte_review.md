---
title: "Simulated Peer Review — Angewandte Chemie submission: A Lineshape-Correlation Coefficient for 2D HSQC Similarity"
author: "Editorial simulation (4 referees + adversarial verification + handling editor)"
date: " "
geometry: margin=1in
colorlinks: true
linkcolor: RoyalBlue
urlcolor: RoyalBlue
fontsize: 11pt
---

# Simulated peer review report

**Manuscript:** *A Lineshape-Correlation Coefficient for Robust Similarity of Two-Dimensional
HSQC NMR Spectra* (D. Lee, KBSI). Files reviewed: `doc/LCC_angewandte.md` (manuscript),
`doc/LCC_angewandte_SI.md` (Supporting Information), and the repository code/benchmarks/tests.

**How this review was produced.** Four independent referees, each with a distinct lens (novelty &
significance; technical & statistical rigor; benchmarking fairness & reproducibility; presentation
& format compliance), read the actual manuscript, SI and code. Every FATAL/MAJOR finding was then
passed to an independent adversarial verifier instructed to *refute* it — to check whether the
defect is real and unrebuttable, or whether the manuscript already discloses/handles it. A
handling editor then weighed the four reports against the verification verdicts, discounting
refuted findings and down-ranking overstated ones. This is a simulation intended to harden the
manuscript before submission, not an actual editorial decision.

---

## Post-review corrections applied (data-integrity + reproducibility)

Three concrete defects surfaced during and after this review have since been **fixed in code**;
the manuscript, SI, `Table 1` and Figure 2 now reflect the corrected numbers:

1. **13C same-compound duplicate (real finding).** The olivetol "same-compound" pair was
   **byte-identical** (`Olivetol.json` ≡ `Olivetol_A.json`, same MD5, identical peak list) — a
   self-comparison scoring 1.00 for every method that inflated `mean_same`. It is **removed**; the
   sparse benchmark now uses **5 genuinely-distinct pairs / 40 different** (was 6/60), and
   `bench_13c.py` now **raises** on any same-pair that renders to an identical image. LCC's sparse
   separation corrects from 0.81 to **0.78**.
2. **`clip`/`abs` baseline bug.** The default `clip` silently deletes negative CH$_2$ crosspeaks in
   multiplicity-edited HSQC (sign encodes multiplicity, not absence). A **`baseline="abs"`** mode was
   added across all methods and the CLI; the peak-list loader already used $|$intensity$|$.
3. **Missing cosine ablation baseline.** The "$0.75 \to 0.59$" mean-centring ablation was asserted
   but not reproducible (no cosine method existed). An **un-centred cosine** (`cosine_similarity`,
   `center=False`) is now a first-class, benchmarked row in both harnesses. It revealed an honest,
   important nuance: mean-centring is decisive **only in the dense protein regime** (0.75 vs 0.59);
   on the sparse regime centred and un-centred coincide (0.778 vs 0.777).

The test suite grew from 19 to **23** (all passing) with one regression per fix. These corrections
address the reproducibility and internal-consistency findings below; the decisive
scope/significance issue (No. 1 in the ranked list) is unchanged.

---

## Editorial decision: **REJECT** (resubmit to a specialized methods journal, *or* return with a genuine chemical-discovery application)

| Referee | Lens | Recommendation |
| --- | --- | :---: |
| Referee 1 | Novelty, significance, journal fit | **Reject** |
| Referee 2 | Technical correctness & statistics | **Major revision** |
| Referee 3 | Benchmarking fairness & reproducibility | **Major revision** |
| Referee 4 | Presentation & Angewandte format | **Minor revision** |
| **Handling editor** | **Consolidated** | **Reject (for Angewandte's bar)** |

**Verification outcome:** of the 16 MAJOR findings the four referees raised, adversarial
verification **confirmed 2** and judged **14 OVERSTATED** (downgraded to MINOR) — because the
manuscript already discloses them, applies its render/window symmetrically to every method, and
scopes its claim to the robust *ranking* rather than the absolute values. In other words, the
paper is **not technically broken**; the reject is a **venue-fit judgment specific to Angewandte's
general-chemistry significance bar**, plus one cheaply-fixable reproducibility gap.

### Editor's rationale

Weighing the four referees against the verification verdicts, the technical case for LCC largely
holds: Referee 2 verified the core math (the $\exp(-\Delta^2/4\sigma^2)$ tolerance factor,
self-similarity $=1$, symmetry, and the mean-centring/covariance argument), and all four referees
confirm the code runs and matches the tabulated numbers. Crucially, the verification pass dissolved
almost every technical MAJOR. The four most-cited objections — sparse-benchmark circularity, the
$n=1$ dense negative and absent error bars, handicapped baselines, and incremental novelty — were
all judged OVERSTATED and revised to MINOR, because the manuscript itself discloses each one,
applies the render/window symmetrically to every method, restricts its claim to the robust ranking
rather than the absolute values, and shows via the blur sweep ($\Sigma = 0.78$–$0.84$) that the
result is not a tuned artifact. So the paper is not technically broken and would be a
minor-revision accept at a specialized methods venue.

Only two findings survived verification as real, load-bearing issues. **First**, the flagship dense
$^1$H–$^{15}$N benchmark is **not independently reproducible** — data only "on reasonable request",
hard-coded private paths (CONFIRMED MAJOR). This is genuine but trivially fixable by a
Zenodo/BMRB deposit. **Second, and decisively**, Referee 1's **scope/significance** finding was
CONFIRMED: the paper produces no chemical result — no binder identified, no CSP-mapped site, no
dereplicated unknown, no QC failure caught — only discrimination scores on two small benchmarks (7
protein comparisons, 66 mostly-synthetic small-molecule pairs). The verification correctly frames
this as a venue-fit judgment rather than a technical defect, but for Angewandte that judgment is
the gate. Even the two "major revision" referees explicitly flagged the venue fit as borderline and
named specialized journals (*J. Magn. Reson.*, *Anal. Chem.*, *J. Cheminform.*) as the natural
home. Being technically clean is not Angewandte's bar; demonstrated broad chemical significance is,
and lifting it requires new science (a real prospective dereplication or CSP-screening result), not
a revision cycle.

---

## Confirmed strengths

- **Mathematically sound core.** Verification confirmed the displaced-Gaussian tolerance factor
  $\exp(-\Delta^2/4\sigma^2)$, self-similarity $=1$, symmetry, and the covariance/mean-centring
  discrimination argument.
- **Reproducible and internally consistent (sparse regime).** The code runs; all headline numbers
  trace correctly to the results JSON and to Tables 1, S2, S3, S4; the sparse $^1$H–$^{13}$C
  benchmark auto-downloads public data and is fully reproducible.
- **Unusually honest self-caveating.** The manuscript discloses the single dense negative, the
  synthetic render, the reimplemented baselines, and restricts its claim to the robust ranking —
  which is *why* most technical objections were judged OVERSTATED.
- **Physically-motivated, non-fitted parameters** (blur $=$ NMR linewidth) plus a robustness sweep
  (Table S5) and a clean ablation (mean-centring $0.75 \to 0.59$) isolating the discriminating step.
- **Genuine cross-regime result.** LCC is the only measure tested that leads in *both* the dense
  protein $^1$H–$^{15}$N and sparse small-molecule $^1$H–$^{13}$C regimes.
- **Clear, well-written, Angewandte format-compliant:** 5 keywords, TOC entry, correctly styled
  references, self-contained captions, two genuinely distinct figures whose bar-chart values match
  Table 1.

---

## Consolidated verified issues (ranked)

| # | Severity | Verified status | Issue | Raised by |
| :-: | :--: | :--: | --- | --- |
| 1 | **MAJOR** | **CONFIRMED** | **No demonstrated chemical discovery** — only similarity/discrimination scores on two small benchmarks; no binder, CSP-mapped site, dereplicated unknown, or QC catch. Below Angewandte's breadth/significance bar; needs new science, not a revision. **This is the decisive, reject-driving issue.** | R1 (also R2, R3) |
| 2 | **MAJOR** | **CONFIRMED** | **Flagship dense benchmark not independently reproducible:** PRL3/OAA spectra only "on reasonable request"; `bench.py` hard-codes private paths. The most striking number (0.75/0.71) cannot be regenerated from a public artefact. Trivially fixable by a DOI deposit. | R3 (also R2, R4) |
| 3 | MINOR | OVERSTATED (was MAJOR) | **No uncertainty quantification:** dense benchmark $n=1$ negative + 6 non-independent positives; sparse benchmark 6 (mostly near-duplicate) positives + 60 non-independent negatives; no CIs/bootstrap/ROC, yet "2.6×/3×" ratios quoted. *Disclosed in text; second benchmark added to broaden negatives; claim scoped to ranking.* Residual valid ask: add bootstrap CIs, more decoy proteins, soften bare multipliers. | R1, R2, R3 |
| 4 | MINOR | OVERSTATED (was MAJOR) | **Baselines are author reimplementations** run in one wide common window with frozen, unswept parameters and no fidelity check vs the originals' published results; "faithful reimplementations" overstates it. *Window applied symmetrically to LCC too; caveat stated twice; saturation argued mechanistically.* Residual ask: sweep each baseline's own knobs, soften "faithful", label "Bin+45°" as this work. | R1, R2, R3 |
| 5 | MINOR | OVERSTATED (was MAJOR) | **Sparse benchmark is synthetic** (peak lists rasterized as Gaussians) and LCC's nominal blur equals the render width on the $^{13}$C axis (potential circularity), with no real measured $^{13}$C-HSQC validation. *Render is shared preprocessing for all five methods; effective widths differ after LCC's second blur; disclosed blur sweep 0.78–0.84 spans widths above and below the render width.* Residual ask: state the $\sigma$-vs-render-width relation; add noise/baseline or a measured pair. | R1, R2, R3 |
| 6 | MINOR | OVERSTATED (was MAJOR) | **Novelty is a recombination** of three prior methods plus a mean-centred correlation (Pearson/ZNCC of blurred images) — the analogue of the mass-spec contrast-angle. *The paper honestly labels every borrowed ingredient and never claims to invent the primitive; the contribution is the physically-motivated construction, backed by the ablation.* | R1 |
| 7 | MINOR | CONFIRMED (as MINOR) | **"Parameters set by spectroscopy, not tuning" is undercut:** $\sigma$'s change per regime and the non-spectroscopic grid step also changes; grid-invariance is only approximate (separation still 0.36 when coarsened). Fix by a stated rule (e.g. $\sigma/3$) or acknowledge the grid as a third parameter. | R3 |
| 8 | MINOR | CONFIRMED (as MINOR) | **Presentation/consistency cluster:** inconsistent multiplier (2.6× separation vs 3× margin) and ablation figure (0.16 vs 0.17); inconsistent hyphenation of "Lineshape(-)Correlation Coefficient"; "ZNCC" undefined outside Figure 1; TOC graphic described but not supplied as a file; lead metric alternates between separation and margin; "strictly decreasing" monotonicity overstated under centring+clamp. | R4 (also R2) |

---

## Revision plan (editor's ordered actions)

1. **Deposit the data.** Put the processed PRL3 and OAA HSQC spectra (or the windowed rendered
   matrices sufficient to rerun `bench.py`) in a public repository (Zenodo/BMRB) with a DOI, and
   de-hardcode the paths in `bench.py` so the dense benchmark reruns from public data. Closes the
   one surviving reproducibility MAJOR.
2. **Add the significance lift Angewandte requires.** Show at least one real chemical result where
   LCC changes a conclusion an incumbent gets wrong — a blind dereplication hit against BMRB/HMDB,
   or a titration/CSP screen where LCC correctly ranks binding/epitope while bins/tree/NN fail. If
   not achievable, redirect to a specialized methods journal (*J. Magn. Reson.*, *Anal. Chem.*, or
   *J. Cheminform.*), where the current evidence base is already sufficient.
3. **Add uncertainty quantification.** Bootstrap confidence intervals on separation/margin over the
   66 sparse-regime pairs; add several decoy proteins (different fold/size) to the dense benchmark;
   report same-vs-different distributions with ROC/AUC. Then soften or drop the bare "2.6×/3×"
   language until it is backed by an interval.
4. **Give the baselines the same latitude LCC receives.** Sweep each baseline's own parameters and
   window and report its best-case separation/margin (or a sensitivity band), and/or validate each
   reimplementation against a published discrimination result. Soften "faithful reimplementations"
   to "representative reimplementations", document departures from the source, and label
   "Bin + 45° rotation" as a variant of this work.
5. **Strengthen the sparse benchmark.** State explicitly the relationship between LCC's blur width
   and the render width and cite the Table S5 sweep; ideally add realistic noise/baseline/lineshape
   distortion or a genuinely measured $^1$H–$^{13}$C HSQC pair, and show the different-class score
   survives when the render width differs from LCC's blur.
6. **Tighten wording/scope claims.** Qualify "strictly decreasing" as weakly monotone under
   mean-centring plus the $\max(0,\cdot)$ clamp; drop the framing that self-similarity $=1$
   distinguishes LCC (all methods share it); either fix the grid step by a stated rule ($\sigma/3$)
   or acknowledge it as a third, non-spectroscopic parameter and report the actual grid dependence.
7. **Clean up presentation.** Unify the multiplier and the mean-centring ablation figure (0.16 vs
   0.17); use one hyphenation of the method name throughout; expand "ZNCC" in the Figure 1 caption;
   produce and submit the actual TOC graphic as an image file; and foreground a single primary
   metric (separation) consistently across the abstract and both result paragraphs.

---

## Bottom line for the author

The manuscript is technically sound and honestly written — the adversarial verification pass is the
proof: 14 of 16 "major" referee objections collapsed to minor once checked against the actual text,
precisely because the paper already discloses its limitations. Two things stand between it and
publication:

- **For any venue:** deposit the dense-benchmark data with a DOI (item 1) and add the quantitative
  rigor of items 3–7. These make it a clean minor-revision accept at a specialized methods journal.
- **For Angewandte specifically:** the paper needs a *chemical result*, not just a better metric
  (item 2). Without a demonstrated discovery, the honest path is a specialized methods journal;
  the Angewandte-format manuscript here transfers directly.

---

# Appendix: full per-referee reports

The complete referee findings, with the adversarial verifier's reasoning for each MAJOR finding,
follow verbatim.

### Referee 1 — REJECT

*Focus: Whether LCC is a genuine conceptual advance and whether a 2D-HSQC similarity metric with no demonstrated chemical discovery meets the breadth/significance bar of a general-chemistry journal.*

**Assessment.** The manuscript is clearly written and the implementation is clean, correct, and honestly caveated. However, judged against Angewandte's bar for a general-chemistry readership, it has a scope/significance mismatch that a revision cannot easily fix. The core construction is, by the authors' own description, a recombination of three existing methods plus a mean-centred correlation — i.e. Pearson correlation of Gaussian-blurred rendered images, a textbook image cross-correlation that the paper itself equates with the long-established mass-spectral contrast-angle/cosine. No chemical discovery is enabled: the paper is a metric benchmark, not a result. The validation is thin (a single different-protein negative; six small-molecule positives, mostly synthetic peak-list rasterizations, of which only two are non-trivial) and the headline "beats every method" rests on the authors' own reimplementations run outside their design envelope. This is strong, useful work for a specialized venue (J. Magn. Reson., Magn. Reson. Chem., J. Cheminform., or Anal. Chem.), but not a broad, significant advance for Angewandte.

**Findings:**

- **[MAJOR] No demonstrated chemical impact — the paper is a metric, not a discovery** — *verification: CONFIRMED → MAJOR*
  - *Location:* MS abstract and Conclusion; Table 1; SI S8/S9
  - The entire result set is discrimination scores on two benchmarks. Nothing chemical is discovered or resolved: no binder identified, no CSP-mapped binding site, no dereplicated unknown, no QC failure caught. The forward-looking language stays hypothetical — 'We expect LCC to be directly useful wherever a robust 2D-HSQC similarity is needed: automated titration tracking and CSP read-out, library and dereplication searches' (MS, penultimate paragraph). Angewandte weighs broad significance for a general-chemistry readership; a similarity number with no downstream chemical result does not clear that bar.
  - *Suggested fix:* Add at least one real application in which LCC changes a chemical conclusion an incumbent method gets wrong — e.g. an actual dereplication hit against BMRB/HMDB on a blind unknown, or a titration where LCC correctly ranks binding affinity/epitope where bins/tree/NN fail. Absent a demonstrated discovery, the natural home is a specialized methods journal.
  - *Verifier:* Every factual assertion in the finding checks out against the manuscript. All results are discrimination numbers: separation/margin on two benchmarks (Table 1, lines 196-202; Abstract lines 37-41). No chemical result is produced from either benchmark — the PRL3 titration (lines 125-133) only supplies "same-class" scores and never maps a binding site or interprets CSP; the simpleNMR set (lines 135-149) only scores same/different pairs and dereplicates no unknown. The paper's own caveat confirms the referee: "LCC measures position-and-intensity coincidence, not molecular identity" (lines 155-157). The forward-looking utility is hypothetical, verbatim as quoted: "We expect LCC to be directly useful wherever a robust 2D-HSQC similarity is needed: automated titration tracking and CSP read-out, library and dereplication searches, and quality control of repeated measurements" (lines 167-170). The referee misread nothing and there is no rebutting hedge or table value; the caveat paragraph reinforces the point. Distinction: this is a broad-significance / journal-fit judgment, not a technical or reproducibility defect. The absence of a chemical result is real and cannot be rebutted, but the author could legitimately argue a methods/metric paper is a valid contribution type — a weaker counter than the referee's factual claim, so the concern stands but as a scope/significance risk rather than an error.
- **[MAJOR] Novelty is an incremental recombination; the scoring step is textbook** — *verification: OVERSTATED → MINOR*
  - *Location:* MS 'We resolve this tension...' and 'Score.' paragraphs; SI S4, S6b, Table S6
  - The paper states LCC 'keeps the physically meaningful ingredients of all three predecessors ... and adds a mean-centred correlation as the discriminating step' (MS) and SI Table S6 tabulates three of four ingredients as 'taken from' prior work. The one new ingredient — 'the mean-centred zero-lag correlation (the discriminating step)' — is the Pearson correlation coefficient of two blurred images, which the authors themselves call 'the spectral-imaging analogue of the cosine / contrast-angle similarity long used in mass-spectral library search' (MS, Score step; refs 22-23). Mean-centring a cosine into a correlation is a standard statistical choice, not a conceptual advance; the ablation (0.75→0.59) merely re-discovers that Pearson separates better than uncentred cosine. ZNCC / template matching on blurred rendered images is established image processing. This is recombination, not a new idea.
  - *Suggested fix:* Either sharpen a genuinely new conceptual claim (what does the NMR setting demand that generic ZNCC/Pearson does not already provide?), or reposition honestly as an engineering synthesis for a methods venue rather than a conceptual advance for a general journal.
  - *Verifier:* The referee's factual reading is correct — the manuscript itself openly frames LCC as a recombination. MS line 74-78: "a measure that keeps the physically meaningful ingredients of all three predecessors ... and adds a mean-centred correlation as the discriminating step." MS line 96-98: "This is the spectral-imaging analogue of the cosine / contrast-angle similarity long used in mass-spectral library search." SI Table S6 lists three ingredients as "taken from" Bodis/Castillo/Pierens and only "mean-centred zero-lag correlation (the discriminating step) | this work." The ablation is exactly as quoted (SI S10: "drops the protein separation from 0.75 to 0.59"). So the referee misread nothing; there is no factual defect to confirm.

But the "MAJOR" severity is overstated for three reasons the manuscript itself supplies. (1) The paper never overclaims the scoring primitive as invented — it explicitly calls it the analogue of an established mass-spec similarity and honestly labels each borrowed ingredient, so there is no misrepresentation to correct. (2) The stated contribution is not "we invented mean-centred correlation" but the whole physically-motivated construction that reconciles two requirements incumbents satisfy only one at a time: MS line 163-165 "a single, physically-motivated construction — render onto a shared grid, blur by the NMR linewidth, and take the mean-centred zero-lag correlation — reconciles shift tolerance with discrimination, the two requirements that the incumbent measures satisfy only one at a time," with knobs that "are spectroscopic quantities rather than fitting parameters." (3) This claim is backed empirically across two regimes: LCC gives the widest same/different margin on both dense protein 1H-15N (margin 0.71 vs bin 0.24) and sparse small-molecule 1H-13C (separation 0.81), the only measure to lead in both (Table 1). A method that measurably outperforms all three standard NMR incumbents in both data regimes, using non-tuned physical parameters, is a legitimate applied-methods contribution for a chemistry venue — "recombination" describes its form, not a deficiency a competent author cannot rebut. The referee's objection is a subjective bar-height/taste concern, appropriately MINOR, not a MAJOR defect.
- **[MAJOR] Dense-regime 'different class' is a single spectrum (n=1 negative)** — *verification: OVERSTATED → MINOR*
  - *Location:* bench.py run(); MS 'Dense protein' paragraph; SI Table S2 (2 vs OAA row)
  - bench.py sets diff = read_bruker_2d(oaa / '103') — one OAA spectrum. Every dense-regime number in Table 1 and Table S2 (separation 0.75, margin 0.71, the '3x the bin method' headline) is defined against that single point: 'separation = mean_same - diff_s', 'margin = min_same - diff_s'. A one-point negative class admits no distribution, no false-positive rate, no significance test; the impressive 0.18 for the different protein could be idiosyncratic to OAA. The claim 'LCC ... pushing the different protein (0.18) below every same-protein score' (SI S8) generalizes from N=1.
  - *Suggested fix:* Expand the negative set to many unrelated proteins (BMRB HSQCs are plentiful) and report a full same-vs-different score distribution with an ROC/AUC, so the discrimination claim rests on a population rather than one comparison.
  - *Verifier:* FACTUAL CORE CONFIRMED: bench.py run() sets `diff = read_bruker_2d(oaa / "103")` — a single OAA spectrum — and computes `"separation": mean_same - diff_s, "margin": min_same - diff_s`, so every dense-regime number in Table 1 (sep 0.75, margin 0.71) and Table S2 (2 vs OAA = 0.1815) is indeed defined against one negative point. A single negative admits no distribution or false-positive rate. The referee reads the code correctly.

BUT THE MANUSCRIPT EXPLICITLY ADDRESSES IT, which the finding ignores. Line 135-137: "To probe the opposite regime — and to replace the single negative spectrum above with a large negative set — we used public HSQC peak lists... and asked each method to score the 6 same-compound pairs above the 60 different-compound pairs." The paper names "the single negative spectrum" as a known limitation and adds Benchmark 2 (60 negatives, Table S3) for exactly this reason. The headline claim is regime-spanning ("the only measure tested that keeps same-class similarity high while pushing different-class similarity down in both regimes", lines 151-152), so it does not rest on the N=1 dense negative alone. SI S10 caveat (iii) even notes "the same-compound side is only 6 pairs," showing the authors track their sample sizes.

The referee's charge that SI S8 "generalizes from N=1" overstates: the sentence "pushing the different protein (0.18) below every same-protein score" is a literal description of the seven observed values in Table S2, not a statistical inference. The "3x the bin method" headline is a ratio of two margins measured on the same single OAA comparison for both methods, so it is an apples-to-apples statement about that one spectrum, not a claimed population effect.

Net: the single-negative dense benchmark is a real limitation, but it is disclosed in the main text and deliberately compensated by a second benchmark with a large negative set. This is not an unaddressed MAJOR flaw that a competent author could not rebut — the author already rebutted it in the submitted text.
- **[MAJOR] Sparse-regime positives are six pairs, mostly synthetic and near-trivial** — *verification: OVERSTATED → MINOR*
  - *Location:* bench_13c.py (SIG_H/SIG_C render vs LCC sigma; PAIRS); SI S9 Table S4; S10 caveat iii
  - bench_13c.py yields n_same=6 (confirmed). Spectra are not measured 2D matrices but 'one Gaussian per peak' rasterizations of peak lists (SI S9, caveat i), and LCC's own blur (sigma_C=0.5) equals the render blur used to build them — LCC scores shapes close to those it would itself generate. Of the six same-compound pairs, the SI concedes 'the olivetol pair is near-identical data (a trivial positive)' (S10 caveat iii) and Table S4 shows olivetol/santonin/chartreusin/indanone are all near-duplicate peak lists; only menthol (two solvents) and rotenone are 'genuinely shifted'. So the meaningful positive set is effectively n=2. Separation 0.81 on n=6 (n≈2 non-trivial) is not a robust significance claim.
  - *Suggested fix:* Benchmark on real measured 1H-13C HSQC matrices (not peak-list rasterizations) with render blur independent of LCC's blur, and with enough genuinely re-recorded same-compound pairs to support a distributional claim.
  - *Verifier:* The referee's raw facts are all accurate and verifiable in the code/SI:

- n_same=6: `bench_13c.py` PAIRS has 6 compounds, `same = [(a[0], b[0]) for _, a, b in PAIRS]` → 6 pairs. Confirmed.
- Rasterization: SI S9 "Peak lists are rasterized to a synthetic 2D spectrum (one Gaussian per peak)"; caveat (i) "spectra are rendered from peak lists, not measured 2D matrices." Confirmed.
- Blur equality: render `SIG_H, SIG_C = 0.03, 0.5` (line 41); LCC in METHODS uses `sigma_f2=0.05, sigma_f1=0.5` (line 49). The 13C axis matches exactly (0.5=0.5); the 1H axis does NOT (0.03 render vs 0.05 LCC). So "sigma_C=0.5 equals render blur" is true only on the 13C axis.
- olivetol trivial: caveat (iii) verbatim "The olivetol pair is near-identical data (a trivial positive)."
- Only 2 genuinely shifted: Table S4 caption verbatim "Two pairs are genuinely shifted (menthol in two solvents; a rotenone re-measurement)." Table S4 LCC scores confirm the other four are near-duplicate (santonin 1.00, olivetol 1.00, indanone 0.99, chartreusin 0.94).

But the referee's two load-bearing INFERENCES are wrong or already handled:

1. The rigging mechanism ("LCC scores shapes close to those it would itself generate") is a misread of the benchmark. The rendered, already-blurred image is the SHARED input to all five methods (bin, bin+45, tree, NN, LCC — `specs` dict, line 96, passed to every `fn`); LCC then applies its OWN blur on top. The render blur is common preprocessing, not an LCC privilege. More decisively, the render blur is applied identically to same-compound AND different-compound pairs, so it cannot manufacture the same-vs-different separation — that separation (LCC mean-diff = 0.003) comes from mean-centring (§S6b, S10 ablation), not from blur matching.

2. The paper never makes the "significance claim" the referee attacks. It reports separation descriptively and explicitly hedges: main text lines 159-160 "the robust ranking, not the exact numerical values, is the result"; SI caveat (i) same; and caveat (iii) already discloses the exact limitation the referee raises verbatim: "the between-compound contrast is strong (60 pairs) but the same-compound side is only 6 pairs." No statistical significance test is run or claimed.

3. The benchmark's stated purpose (main text line 135-136) is "to replace the single negative spectrum above with a large negative set" — the discrimination story rests on the well-powered 60-pair negative side and on beating tree/NN's saturation (a 60-pair result), not on the thin positive set. And this is the secondary benchmark; the headline quantitative claim is the protein benchmark (0.75).

So the residual legitimate point — the positive set is small and partly near-duplicate — is real but already disclosed by the authors, and the "blur rigging" and "not a robust significance claim" framing that earns the MAJOR label is a misread of what the benchmark tests and what the paper claims.
- **[MAJOR] The 'widest separation of every method tested' claim rests on handicapped reimplementations** — *verification: OVERSTATED → MINOR*
  - *Location:* MS abstract and Table 1; SI S10 caveat ii; bench.py/bench_13c.py METHODS (per-regime LCC sigma vs fixed baseline params)
  - The three baselines are the authors' own reimplementations run in 'one common wide window' that the authors admit is outside the competitors' design: 'the original papers tune per-spectrum windows and parameters, so this is not a verdict on those methods as their authors deployed them' (SI S10 caveat ii; echoed MS final caveat). Meanwhile LCC's two knobs are tuned per regime (sigma 0.03/0.30 for 15N vs 0.05/0.50 for 13C). Beating tree/NN when they are run as 'drop-in global similarity scores' they were not designed for does not establish LCC's superiority over the published methods, yet the abstract asserts LCC 'delivers the widest same/different separation of every method tested'.
  - *Suggested fix:* Either compare against each method under its authors' own recommended windowing/parameters, or soften the headline to a same-window controlled comparison and drop the implication of a general ranking. Give the baselines the same per-regime tuning latitude LCC receives.
  - *Verifier:* The referee's factual premises all check out, but the manuscript already discloses and bounds every one of them, so the "overclaim" charge is overstated.

CONFIRMED facts: (1) LCC sigma is changed per regime — code shows bench.py uses defaults sigma_f2=0.03/sigma_f1=0.30 (via WINDOW only) while bench_13c.py line 49 passes sigma_f2=0.05, sigma_f1=0.5. (2) The three baselines are the authors' reimplementations run in one common wide window — MS lines 158-161 and SI S10 caveat (ii) both say so. The referee even acknowledges these are "echoed" in the paper.

Why OVERSTATED, not a MAJOR defect:

1. The abstract is already scoped to what the referee complains is untested-as-deployed: "LCC delivers the widest same/different separation of every method tested" (line 37). It does NOT claim superiority over the published methods as their authors deployed them.

2. The paper states the exact caveat the referee demands, twice. MS lines 158-161: "the three reference methods are our own faithful reimplementations compared inside one common window... accordingly the robust ranking, not the exact numerical values, is the result, and this is not a verdict on the original methods as their authors deployed them with per-spectrum windows." SI S10(ii) repeats it. The referee is asking the paper to add a disclaimer it already contains.

3. The "LCC tuned vs baselines frozen" asymmetry is weaker than claimed. The bin baseline IS re-parameterized per regime too — bench_13c.py line 45 passes min_bin_width_f2=0.1, min_bin_width_f1=1.0 for the 13C regime, differing from bench.py's defaults. So baselines are not uniformly frozen while only LCC is adjusted.

4. The LCC per-regime change is defended as physical (sigma = sqrt(linewidth^2 + drift^2), SI S3/S7) rather than fit-to-win, and SI S10 supplies a robustness sweep showing the advantage is not a tuned artifact: "separation is Sigma = 0.71-0.77 at sigma = 0.02-0.04 / 0.20-0.40 ppm, and still Sigma = 0.36 (vs the bin method's 0.29) when the grid is coarsened all the way to the bin method's own resolution." Table S5 shows the same for 13C. So calling the two knobs "tuned per regime" while the paper shows the ranking survives a wide sweep is a mischaracterization.

5. The tree/NN saturation is argued mechanistically ("in any crowded window every peak finds a coincidental near neighbour," MS lines 67-69; SI S6, S9), i.e. intrinsic to drop-in global scoring, not merely a parameter-choice artifact.

Residual kernel of truth (why not fully REFUTED): tree and NN get no per-regime parameter adjustment while LCC's blur is, and the abstract's confident "widest separation of every method tested" could read as stronger than the caveated body. That is a legitimate framing nit, but it is MINOR, not MAJOR — a competent author rebuts it by pointing to the "every method tested" scoping, the two explicit caveats, the regime-parameterized bin baseline, and the robustness sweep.
- **[MINOR] Classification framing without any classification statistic**
  - *Location:* MS 'We summarize discrimination by...'; Table 1; SI Tables S2, S3
  - The paper motivates 'margin (worst same-class - best different-class score, which is what a single classification threshold actually sees)' (MS) and repeatedly frames a same/different decision, yet reports no ROC curve, AUC, sensitivity/specificity, or confidence interval anywhere in MS/SI. 'Margin' from n=1 (protein) or n=6 (small molecule) is not a substitute for a classification-performance measure.
  - *Suggested fix:* Report AUC with a bootstrap CI on adequately sized same/different populations; this is also the natural way to make the significance claim quantitative for reviewers.

### Referee 2 — MAJOR REVISION

*Focus: Correctness of the ZNCC/Pearson math and tolerance derivation; statistical validity of the separation/margin claims given tiny, non-independent positive sets and a single negative; circularity of the synthetic sparse benchmark; fairness of the baselines; reproducibility; and scope/significance for a general-chemistry journal.*

**Assessment.** The core construction is mathematically sound: I verified the exp(-Δ²/4σ²) tolerance factor is the correct normalized overlap of two displaced Gaussians of width σ, self-similarity=1 and symmetry follow trivially from the Pearson form in _zncc, and the covariance/mean-centring argument in S6/S10 is correct. The claimed numbers match results/method_comparison.json and results/comparison_13c.json. The problem is not the math, it is the evidence base. The headline "2.6×/3×" superiority rests on a dense benchmark with exactly ONE negative spectrum and six non-independent positives (a single PRL3 titration series), with no confidence intervals or significance tests anywhere; the sparse benchmark is entirely synthetic (peak lists rasterized as Gaussians) and, worse, LCC's blur width σ_13C=0.5 is set equal to the very render width SIG_C=0.5 that generated the data, a home-field match the paper does not disclose. The three baselines are the authors' own reimplementations run in a single wide window the authors admit is "not how their authors deployed them," so "beats prior methods" is partly a strawman. For a general high-impact journal these gaps are serious; the method may be better suited to a specialized methods journal unless the statistical and reproducibility case is substantially strengthened.

**Findings:**

- **[MAJOR] Headline superiority rests on n=1 negative and six non-independent positives; no error bars or significance tests** — *verification: OVERSTATED → MINOR*
  - *Location:* SI Table S2 (L238), bench.py SAME=[4,6,8,10,12,14]; MS Table 1, MS L129, SI L242
  - The dense benchmark's entire 'different' class is a single spectrum: Table S2 lists one negative, '2 vs OAA (diff. protein)' = 0.1815, and the six positives are one PRL3 ligand titration series (SAME = [4,6,8,10,12,14] in bench.py), i.e. the same protein tracked across a titration — highly correlated, not independent replicates. Separation 0.75 and margin 0.71 (Table 1) and the '≈2.6×'/'three times' claims (MS L129; SI L242) are therefore ratios of point estimates with effective sample size ~1 vs ~1. No confidence interval, bootstrap, or significance test appears anywhere in MS or SI. A ratio of two single-sample means cannot support a quantitative superiority claim in a demanding general journal.
  - *Suggested fix:* Add multiple negative proteins (several distinct apo spectra) and, ideally, independent same-protein replicates rather than one titration series; report bootstrap CIs on separation/margin and a rank-based significance statistic (e.g. AUC with CI, or Mann–Whitney) for each method. Downgrade the '2.6×/3×' language until it is backed by an interval.
  - *Verifier:* The referee's factual claims are almost all verified. (1) Dense-benchmark negative class = a single spectrum: Table S2 (SI L238) lists exactly one different-class row, "2 vs OAA (diff. protein) = 0.1815". (2) The six positives are one PRL3 ligand titration series: bench.py L22 `SAME = [4, 6, 8, 10, 12, 14]  # PRL3 + ligand titration`, REF=2, all the same protein tracked across increasing ligand — correlated, not independent replicates. (3) The headline dense numbers are point-estimate ratios: MS L129 "a margin of 0.71 — three times the bin method's 0.24"; SI L242 "≈2.6× better". (4) No error bar, bootstrap, CI, std, or significance test appears anywhere — a grep for bootstrap|confidence|significan|p-value|std|error bar|standard error|t-test|resampl across MS, SI, bench.py, bench_13c.py returns nothing. So "no error bars or significance tests," "n=1 negative," and "six non-independent positives" are all literally true.

BUT the referee's framing — "Headline superiority rests on n=1 negative" — is overstated, and the manuscript already pre-empts the core concern. MS L135–136 introduces the second benchmark explicitly: "To probe the opposite regime — and to replace the single negative spectrum above with a large negative set — we used public HSQC peak lists ... six compounds each recorded twice ... score the 6 same-compound pairs above the 60 different-compound pairs." That is the author openly acknowledging the dense benchmark's single negative AND remedying it with a 66-pair benchmark where LCC also leads (0.81 sep, 0.37 margin; Table 1, Table S3). The abstract headline is dual-regime ("0.75 and 0.81"), so overall superiority does not rest on n=1 alone. The author further hedges to ranking rather than values: MS L159–161 "the robust ranking, not the exact numerical values, is the result," and SI S10 caveat (iii) L313–314 "the between-compound contrast is strong (60 pairs) but the same-compound side is only 6 pairs." A competent author can therefore rebut the "rests on n=1" thesis by pointing to the second regime and the explicit ranking framing. What the author cannot rebut is the total absence of any uncertainty quantification and the bare "2.6×"/"three times" ratio language, which does overclaim precision from a single negative and should be softened; a bootstrap over the 66 sparse-regime pairs would be feasible and would strengthen the paper.
- **[MAJOR] Sparse benchmark is circular: LCC blur width equals the Gaussian width used to synthesize the data** — *verification: OVERSTATED → MINOR*
  - *Location:* bench_13c.py L41 vs L49; SI S9/S10 caveat (i) L308
  - bench_13c.py rasterizes peak lists with SIG_H, SIG_C = 0.03, 0.5 (L41) and then scores LCC with sigma_f2=0.05, sigma_f1=0.5 (L49). The 13C blur σ_C=0.5 is identical to the render width SIG_C=0.5 that generated the synthetic spectra, and σ_H is within a hair of SIG_H. LCC's model (correlate Gaussian lineshapes) exactly matches the data-generating process, an advantage the position-insensitive/peak-based baselines do not get. The SI discloses that spectra are 'rendered from peak lists' (caveat i, SI L308) but never discloses that LCC's tuning parameter coincides with the generator's width. This directly inflates the headline sparse separation (0.81).
  - *Suggested fix:* Either use real measured 2D 13C-HSQC matrices, or generate the synthetic spectra with a lineshape/width different from (and unknown to) LCC's blur, and show the separation survives a mismatch. State the σ-vs-render-width relationship explicitly.
  - *Verifier:* The referee's factual observation is partly true but the "circularity inflates the headline" conclusion is refuted by the manuscript's own code and SI.

WHAT'S TRUE: bench_13c.py L41 sets render widths `SIG_H, SIG_C = 0.03, 0.5` and L49 sets `lcc_similarity(... sigma_f2=0.05, sigma_f1=0.5 ...)`. The nominal σ_13C=0.5 equals SIG_C=0.5. The SI does not explicitly flag this numeric coincidence.

WHY THE "CIRCULARITY / INFLATION" CLAIM FAILS:

(1) The render is SHARED preprocessing for ALL methods, not an LCC-only advantage. load_peaklist (L80) blurs the image once — `img = _blur(_blur(img, SIG_C/C_STEP, 0), SIG_H/H_STEP, 1)` — and returns ONE Spectrum2D. run() (L104-106) passes that same `specs[a]` to every method (bin, bin_rot45, tree, nn, lcc). The Gaussian lineshape is baked into the images all five methods consume; it is not something LCC uniquely matches while baselines don't see it. The peak-based baselines pick peaks back out of the same rendered image, and peak centers are recoverable from a Gaussian at any width.

(2) The numeric "match" is not even real at the operative level. hsqc_lcc.render_image ALWAYS re-blurs, so LCC applies its σ on top of the already-blurred image: effective 13C width = sqrt(0.5^2 + 0.5^2) = 0.707 ppm, and 1H = sqrt(0.03^2 + 0.05^2) = 0.058 (referee's "within a hair" is off — nominal 0.05 vs 0.03 is 67% larger). LCC's effective lineshape does NOT equal the generator width.

(3) The manuscript directly tests sensitivity to the coincidence and discloses it. SI Table S5 (L297-306) sweeps ONLY LCC's blur (render SIG fixed in bench): separation = 0.78 at σ=0.03/0.3, 0.81 at 0.05/0.5, 0.84 at 0.12/1.2. At σ_13C=0.3 — well BELOW the render width 0.5 — separation is 0.78, only 0.03 under the headline 0.81; at 1.2 (far above) it rises to 0.84. So de-matching the blur from the render width does not collapse the result; the headline is not an artifact of the coincidence. Text L293-296 makes the robustness claim explicitly.

(4) Caveat (i) (SI L308-310) discloses the underlying dependence: "spectra are rendered from peak lists, not measured 2D matrices, so absolute values depend on the render width; the ranking is the result, not the absolute number." The claim rests on ranking/separation ordering, not the absolute 0.81.

Net: the referee spotted a genuine parameter coincidence and a small disclosure gap (one sentence noting σ nominally equals the render width, and pointing to Table S5), but the asserted consequence — a circular benchmark that "directly inflates" the headline — is contradicted by the shared-render design, the convolution making effective widths unequal, and the disclosed blur sweep showing 0.78–0.84 across widths spanning above and below the render width.
- **[MAJOR] Neither benchmark is simultaneously real, public, and reproducible** — *verification: OVERSTATED → MINOR*
  - *Location:* MS L223–224 Data Availability; bench.py L65–66
  - The only real measured 2D data (dense PRL3/OAA) is not public — 'available from the author on reasonable request' (MS L223–224) — and the hard-coded paths in bench.py (L65–66: /Users/donghanlee/PRL3_mark_5nov24, /Users/donghanlee/Downloads/OAA...) confirm it cannot be reproduced by a referee. The public benchmark (sparse) is synthetic (peak lists, not measured matrices). So there is no reproducible test of LCC on real, measured 2D-HSQC matrices with public data. For a general journal that weighs reproducibility heavily, the primary quantitative result is unverifiable.
  - *Suggested fix:* Deposit at least one public set of real measured 2D-HSQC matrices (e.g. a BMRB/Zenodo deposition of the PRL3 titration and a negative protein) so the dense benchmark is independently reproducible.
  - *Verifier:* All factual claims verify. MS L223–224: "The dense 1H–15N benchmark uses PRL3 and OAA Bruker spectra available from the author on reasonable request" — not public. bench.py L65–66: hard-coded defaults Path("/Users/donghanlee/PRL3_mark_5nov24") and Path("/Users/donghanlee/Downloads/OAA_CEST_277K_04may15") — confirmed (argparse defaults, but moot without data). Sparse benchmark is synthetic: bench_13c.py L10–12 "peak lists are rasterized to a synthetic 2D spectrum (one Gaussian per peak)"; MS L159 "the 1H–13C spectra are rendered from peak lists." So the structural gap is real: no benchmark is simultaneously real-measured + public + reproducible. However, the referee's impact claim ("the primary quantitative result is unverifiable") is overstated. (1) The sparse 1H–13C regime IS fully reproducible — bench_13c.py L26 auto-downloads public simpleNMR peak lists and runs with public code, so a referee can independently reproduce the 0.81 separation and the entire cross-regime ranking in one regime. (2) The manuscript explicitly frames the result as the ranking not the values (L159–161: "the robust ranking, not the exact numerical values, is the result") and discloses the synthetic-render caveat honestly (L155–161, "Two caveats bound the claim honestly"). (3) "Available on reasonable request" is a standard data-availability mode. The narrow true point — the real-measured 2D result cannot be independently reproduced — stands and is a legitimate, cheaply-fixable reproducibility weakness (deposit the Bruker data, de-hardcode paths), but it does not render the primary result "unverifiable."
- **[MAJOR] Baselines are the authors' reimplementations run in a regime the authors admit is not how they were deployed** — *verification: OVERSTATED → MINOR*
  - *Location:* MS L157–161; SI S10 caveat (ii) L310–312
  - The tree and NN 'saturation' that motivates the whole paper is produced by the repository's own reimplementations 'run inside one common wide window,' whereas 'the original papers tune per-spectrum windows and parameters' (SI L310–312; MS L158–161). The reimplementation fidelity is not independently validated against the original authors' published results. The central narrative — that incumbent methods 'leave the two classes overlapping' (MS L38) — may therefore be an artifact of the chosen single-window drop-in setting rather than a property of the methods. Combined with the 'ranking, not exact values' hedge, the comparison establishes less than the abstract claims.
  - *Suggested fix:* Validate each reimplementation by reproducing a published discrimination result from the original papers, or benchmark against the authors' original code. Soften abstract/intro claims that these methods 'fail' to a claim about their behavior as global drop-in scores.
  - *Verifier:* The referee's underlying facts are all true and, crucially, all disclosed by the manuscript itself. MS L158-161: "the three reference methods are our own faithful reimplementations compared inside one common window, and the $^1$H–$^{13}$C spectra are rendered from peak lists; accordingly the robust **ranking**, not the exact numerical values, is the result, and this is not a verdict on the original methods as their authors deployed them with per-spectrum windows." SI S10 caveat (ii) L310-312 says the same: "the original papers tune per-spectrum windows and parameters, so this is not a verdict on those methods as their authors deployed them — only on how they behave as drop-in global similarity scores." So the referee is essentially quoting the manuscript's own caveats back at it and re-labeling them a MAJOR flaw.

The referee's strongest sub-point — "reimplementation fidelity is not independently validated against the original authors' published results" — is factually correct: the repo contains the reimplementations (hsqc_methods.py, hsqc_similarity.py) but no test compares them to published numbers. That is a genuine, undisputed gap. However, it does not carry the weight the referee assigns. The paper's stated goal is a parameter-free, single-window, drop-in measure; running every method in one common window WITHOUT per-spectrum tuning is the fair, apples-to-apples comparison for that goal, and the manuscript is explicit that this is the scope. The abstract claim "the tree and nearest-neighbour measures leave the two classes overlapping (margin ≤ 0.04)" is a factually accurate report of the benchmark (Table 1 / Table S2: margins −0.01, 0.03, 0.00, 0.04) about the drop-in regime, not a claim about the original methods as deployed.

The referee's assertion that the saturation "may be an artifact of the single-window drop-in setting rather than a property of the methods" is also independently rebutted by the mechanistic argument in S6b/S9 (NN's one-way nearest-neighbour distance and the tree's mass-centre overlap inherently saturate in any crowded/wide window because every peak finds some near neighbour). That argument holds regardless of reimplementation fidelity.

So: real underlying limitation, fully disclosed, with the claim explicitly scoped and hedged ("ranking, not exact values"; "not a verdict on the original methods"). This is not a claim-invalidating MAJOR that a competent author cannot rebut — the author can point to two explicit caveats and the parameter-free framing. It survives as at most a minor "add an independent fidelity check against published numbers" request.
- **[MAJOR] Scope/significance for a general chemistry journal** — *verification: OVERSTATED → MINOR*
  - *Location:* MS abstract and L98; overall scope
  - The contribution is an image-processing similarity metric (render, Gaussian blur, mean-centred zero-lag Pearson) — the paper itself calls it 'the spectral-imaging analogue of the cosine / contrast-angle similarity long used in mass-spectral library search' (MS L98). Validated on 7 protein comparisons and 66 synthetic small-molecule pairs, it is a competent but incremental methods advance. For Angewandte (≈30% acceptance, breadth expectation) the demonstrated impact is narrow; a specialized venue (J. Magn. Reson., J. Cheminform., Anal. Chem.) would be a more natural home unless the significance is broadened (e.g. a real prospective screening/dereplication application at scale).
  - *Suggested fix:* Either add a substantially larger, real-world application demonstrating broad impact, or target a specialized methods journal.
  - *Verifier:* The referee's factual predicates are all accurate. L98 does read verbatim: "This is the spectral-imaging analogue of the cosine / contrast-angle similarity long used in mass-spectral library search." The benchmark scale is as stated: the protein regime is one apo reference vs "six spectra of the same protein along a ligand titration ... and with a **different** protein" (7 comparisons), and the small-molecule regime is "6 same-compound pairs above the 60 different-compound pairs" (66 pairs). So the referee did not misread the numbers or the self-description.

But three things make the MAJOR framing overstated:

1. It is an editorial venue-fit opinion, not a defect grounded in a manuscript error. The referee catches no overclaim: the paper never asserts large-scale or prospective validation. It explicitly frames the two datasets as "two deliberately **opposite** regimes, so that any measure claiming general utility must succeed on both" (L119) — a breadth-by-design argument, not a scale claim — and bounds impact honestly with two caveats (L155–162): LCC "measures position-and-intensity coincidence, not molecular identity," and "the robust **ranking**, not the exact numerical values, is the result." There is nothing here the author is caught oversically claiming.

2. The "incremental / just cosine" characterization understates the manuscript's explicitly stated and ablated contribution. The paper positions the mean-centring (Pearson, not cosine) as the novel discriminating operation: "Mean-centring is what turns overlap into discrimination ... intensity that is not co-located is actively penalized rather than merely ignored, as it is by the one-way nearest-neighbour distance" (L98–103), and backs this with an ablation — "replacing it with an un-centred cosine overlap on the *same* images drops the separation from 0.75 to 0.59" (L131–132), also headlined in the abstract. The referee's own L98 quote is the sentence the paper immediately qualifies as differing from cosine by the centring step; calling it "the spectral-imaging analogue of cosine" therefore undersells what the paper demonstrates.

3. The unrebutted core — no prospective, at-scale screening/dereplication application, and modest demonstrated impact — is factually true and grounded in the text, so the concern is not empty. That part a competent author cannot make disappear; they can only argue relevance. This keeps it from being REFUTED.

Net: the observation is real and grounded, but (a) it is a subjective journal-fit judgment rather than a scientific flaw, and (b) it partly rests on downplaying the manuscript's demonstrated novelty. That is a legitimate MINOR editorial point an editor may weigh, not a MAJOR defect the authors could not rebut.
- **[MINOR] Self-similarity=1 is promoted as a distinguishing property but is shared by every method**
  - *Location:* MS L109, L123, L166–167; SI S5
  - The MS repeatedly sells 'self-similarity is exactly one' (MS L109, L166–167; abstract), yet the paper also states 'Every method self-scores exactly 1.00' (MS L123; Table S2/S3). A property held by all compared methods carries no discriminating information and should not be framed as an LCC advantage.
  - *Suggested fix:* Keep self-similarity=1 as a sanity property but drop the framing that it distinguishes LCC.
- **[MINOR] Monotonicity is stated as 'strictly decreasing' but only holds in the un-centred idealization, and the numerical check is thin**
  - *Location:* SI S5 L153–155, Eq. S9; hsqc_lcc.py L102; tests/test_hsqc_lcc.py
  - SI S5 calls S 'a strictly decreasing function of the displacement magnitude (Eq. S9)', but Eq. S9 is the UN-centred Gaussian overlap; under mean-centring plus the max(0,·) clamp in _zncc (hsqc_lcc.py L102) the centred correlation can flatten to 0 for large Δ, so strict monotonicity fails there. The test (test_lcc_monotonic_in_drift) checks only 5 coupled drift points with a >= tolerance, not strict decrease.
  - *Suggested fix:* State that monotonicity holds for the un-centred/large-grid limit and is weakly monotone under centring+clamp; qualify 'strictly'.
- **[MINOR] Inconsistent value for the mean-centring ablation gain (0.16 vs 0.17)**
  - *Location:* SI L193; MS abstract/L131; hsqc_lcc.py L17
  - The ablation gain is quoted as '≈ +0.16' in SI S6b (L193) and as '0.75 to 0.59' (=0.16) in the abstract/MS, but the code comment states '~+0.17 separation step' (hsqc_lcc.py L17). Minor but should be internally consistent, and the 0.59 un-centred value is asserted without appearing in the results JSON I could locate.
  - *Suggested fix:* Reconcile to a single figure and add the un-centred-cosine ablation output to the results directory so it is reproducible.

### Referee 3 — MAJOR REVISION

*Focus: Whether the head-to-head comparison against the three literature baselines is fair, whether the two benchmarks are independently reproducible, and whether the "ranking, not values" hedge covers the construction and statistical weaknesses.*

**Assessment.** The LCC construction is clean, physically motivated and the code runs and matches the tabulated numbers. My concern is not the method but the evidence offered for its superiority. The head-to-head win rests on an evaluation design that structurally favours LCC: the baselines are the authors' own reimplementations run in one deliberately wide common window with frozen, arbitrary parameters and no sensitivity analysis, while LCC's knobs are both swept (S10) and re-tuned per regime; the sparse benchmark data are generated by one-Gaussian-per-peak rendering, which is LCC's own model; the dense-protein headline separation is (mean of six) minus a single negative spectrum; and the flagship protein raw data are not available (only "on reasonable request"), so half the paper — the more impressive half — cannot be reproduced. The authors' "robust ranking, not the exact numerical values, is the result" hedge does not cover these issues, because the ranking itself is what the design biases. These are addressable in revision, but as it stands the central claim that "LCC is the only measure tested that keeps same-class high while pushing different-class down in both regimes" is not fairly supported. Scope/significance for a general-chemistry readership is also borderline.

**Findings:**

- **[MAJOR] The baselines are handicapped, and the "ranking not values" hedge does not cover it** — *verification: OVERSTATED → MINOR*
  - *Location:* MS main text ("LCC is thus the only measure..."; "the robust ranking, not the exact numerical values, is the result"); bench.py L21,L26-33; bench_13c.py L44-50; hsqc_methods.py L201-215, L302-313; SI S10 caveat (ii)
  - The paper's headline is a ranking claim: "LCC is thus the only measure tested that keeps same-class similarity high while pushing different-class similarity down in both regimes." Yet the baselines are run in a single deliberately wide window (bench.py L21 range_f2=(6.5,10), range_f1=(105,130); bench_13c.py WIN 0-10 x 0-165) with hard-coded arbitrary parameters (quadtree_similarity alpha=0.1, beta=0.33, gamma=3.0, threshold_frac=0.005; nn_peak_similarity threshold_frac=0.05) that are never justified and never varied. Meanwhile LCC receives a full robustness sweep (Table S5, S10) AND is re-tuned between regimes (S7: sigma 0.03/0.30 -> 0.05/0.5, and bench_13c L49 also changes step 0.01/0.10 -> 0.02/0.2). The tree/NN saturation the paper attributes to the methods is a direct consequence of the wide-window choice the authors admit the originals avoid: "the original papers tune per-spectrum windows and parameters" (S10 caveat ii). The hedge that "the robust ranking, not the exact numerical values, is the result" (main text) is self-defeating: the ranking is exactly the quantity the asymmetric evaluation biases. This is a fairness problem, not a values problem.
  - *Suggested fix:* Give each baseline the same courtesy LCC gets: sweep each baseline's own parameters and window and report its best-case separation/margin, or at minimum a sensitivity band. If LCC still wins against each baseline at that baseline's best setting, the ranking claim becomes defensible. Otherwise soften the claim from "only measure that works" to "works out-of-the-box with physically-set knobs, unlike methods that require per-spectrum tuning."
  - *Verifier:* The referee's factual claims are all accurate: bench.py L21 uses one wide window (6.5-10 x 105-130); bench_13c.py WIN is 0-10 x 0-165; quadtree defaults are alpha=0.1/beta=0.33/gamma=3.0/threshold_frac=0.005 (hsqc_methods.py L207-210); nn_peak threshold_frac=0.05 (L308); LCC is re-parameterized between regimes (defaults 0.03/0.30 in bench.py vs sigma_f2=0.05/sigma_f1=0.5/step_f2=0.02/step_f1=0.2 in bench_13c.py L49); and LCC alone gets robustness sweeps (Table S5, S10). But the MAJOR fairness/severity framing is overstated for three text-grounded reasons. (1) The window is applied SYMMETRICALLY: bench.py _methods() passes **WINDOW to lcc_new too (L29-36), and bench_13c METHODS passes **WIN to all five including lcc_new (L45-49). The referee's central lever -- 'tree/NN saturation is a direct consequence of the wide-window choice' -- backfires: LCC runs in the same wide window and discriminates, which is the paper's demonstrated point, not a bias applied only to baselines. (2) The claim is explicitly and repeatedly bounded: main text 'the robust ranking, not the exact numerical values, is the result, and this is not a verdict on the original methods as their authors deployed them with per-spectrum windows'; SI S10 caveat (ii) 'the original papers tune per-spectrum windows and parameters, so this is not a verdict on those methods ... only on how they behave as drop-in global similarity scores.' The ranking is scoped to drop-in global scores in one common window with default parameters -- fair for the claim actually made. (3) LCC's re-parameterization is a physical regime change (15N->13C, S7 'blur widths are scaled to that regime') and is shown tuning-INSENSITIVE by the very sweep the referee calls unfair: S5/S10 separation is flat (0.71-0.77; 0.78-0.84), and the step change is grid resolution shown grid-invariant (S10 separation 0.36 even at bin resolution) -- not a discrimination knob. Corroborating: the same fixed baseline parameters give bin 0.29 (protein) vs 0.69 (sparse), so baselines are NOT uniformly handicapped; tree/NN saturate in both regimes with fixed parameters, supporting the paper's structural (parameter-independent) saturation mechanism in S6. Residual valid point: the paper never sweeps the baselines' own knobs (gamma, threshold_frac) to prove the saturation is parameter-independent as asserted -- a real completeness gap, but MINOR robustness-reporting, not a MAJOR flaw that invalidates the bounded ranking claim.
- **[MAJOR] Sparse-regime data are rendered by LCC's own generative model (circularity)** — *verification: OVERSTATED → MINOR*
  - *Location:* bench_13c.py L65-89 (load_peaklist, _blur) and L49; SI S9 ("Peak lists are rasterized ... one Gaussian per peak"); Table S3 LCC mean diff 0.003
  - In bench_13c.py load_peaklist (L65-81) every peak list is rasterized to a synthetic image as one Gaussian per peak (img[...] += inten then _blur with SIG_H=0.03, SIG_C=0.5), and ALL methods score this same synthetic image. LCC (L49, sigma_f1=0.5) is a render-Gaussians-and-correlate method, so the data-generating process and the method under test share the same model. Two different clean sums of narrow non-overlapping Gaussians are near-orthogonal after mean-centring, which is precisely why LCC reports mean different = 0.003 (Table S3) — an artefact of noiseless, baseline-free, artefact-free synthetic data, not empirical performance. The SI acknowledges rendering (S10 caveat i) but frames it only as a values caveat; the deeper issue is that the construction favours the method being validated.
  - *Suggested fix:* Validate the sparse regime on real measured 2D 13C-HSQC matrices (not peak lists), or at least add realistic noise/baseline/lineshape distortion to the render and show the 0.003 different-class score survives. Report how LCC's lead changes when the render width differs from LCC's blur width (the current bench uses matched widths).
  - *Verifier:* The factual substrate is correct: bench_13c.py L65-81 rasterizes each peak list to one Gaussian per peak (img[...] += inten, then _blur with SIG_H=0.03, SIG_C=0.5, L41/L80), and LCC (L49, sigma_f1=0.5) is a render-Gaussians-and-correlate method. But the "circularity favours the method being validated" framing does not survive the code and the manuscript's own hedges.

(1) The render is method-AGNOSTIC, not LCC-private. All five methods consume the identical rasterized Spectrum2D (bench_13c.py L96, L102-104). The bin method (hsqc_similarity) is an equally position-sensitive rasterize-and-compare method and it also wins decisively on this data (Table S3: separation 0.692, margin 0.195) — the manuscript openly states "The bin method is far stronger here (0.69) than on the dense fingerprint (0.29)". If Gaussian rendering created circular advantage for LCC specifically, the co-competitor bin method would not be similarly favored. The tree/NN methods fail (margins 0.004, 0.043) because they are shift-BLIND (hsqc_methods.py: quadtree exp(-gamma*dist), NN 1/(1+d)), not because of the render. The benchmark tests position-sensitivity vs shift-tolerance, which is the actual claim.

(2) The near-orthogonal "diff = 0.003" is the discrimination signal, not an artefact. Same-compound pairs are rendered by the exact same Gaussian generator yet score mean 0.817 (Table S3) because their peaks co-locate. If summing narrow Gaussians trivially forced orthogonality, same-compound pairs would also collapse to ~0 — they do not. The referee's "two different Gaussian sums are near-orthogonal" is precisely the intended mechanism working.

(3) LCC does not privilege-read the generative model. The generative blur (SIG_H=0.03/SIG_C=0.5 on a 0.01/0.10 grid) is followed by LCC's OWN independent regrid (render_image step_f2=0.02/step_f1=0.2) and OWN second blur (sigma_f2=0.05/sigma_f1=0.5, hsqc_lcc.py L48-89), so LCC re-derives its image from the pixels like every other method — it is not fed the peak list directly.

(4) The primary validation (Benchmark 1, dense protein) uses REAL measured Bruker 2D matrices, not rendered peak lists (manuscript L224 "OAA Bruker spectra available from the author"; SI S8 PRL3/OAA experiments), and LCC wins there too (separation 0.75, margin 0.71, Table S2). The central claim rests on real data; the 13C set is an explicit complementary cross-regime probe.

(5) The manuscript already discloses the render limitation and restricts the claim to the ranking, not the values: main text "the 1H-13C spectra are rendered from peak lists; accordingly the robust ranking, not the exact numerical values, is the result"; SI S10 caveat (i) "absolute values depend on the render width; the ranking is the result, not the exact numbers".

The referee's genuine kernel — that noiseless, baseline-free synthetic images make the exact 0.003 unrealistically clean — is real but is exactly the value the authors decline to claim; it warrants at most a one-line acknowledgement that real noise/baseline/t1-ridges would raise the different-class floor. The stronger "construction favours the validated method" charge is refuted by (1)-(4), so this is not a MAJOR defect.
- **[MAJOR] Dense-protein headline rests on a single negative spectrum** — *verification: OVERSTATED → MINOR*
  - *Location:* bench.py L23-24, L44, L57; MS Table 1 (H-15N sep 0.75 / margin 0.71); SI Table S2 row "2 vs OAA"; bench_13c.py L98-99
  - The flagship result — separation 0.75, margin 0.71 (Table 1, Table S2) — is (mean of six same-protein scores) minus exactly ONE different-protein score: OAA exp 103 (bench.py L23 diff = read_bruker_2d(oaa/"103"); L57 diff_s single value). Margin = min_same - diff is likewise governed by one number. There are no error bars, no multiple decoy proteins, no negative controls of varying similarity. A metric's discrimination cannot be established against n=1 negative. The sparse benchmark's 60 "different" pairs are also non-independent: they are all pairwise combinations of just 12 spectra from 6 compounds (bench_13c.py L99), so each spectrum recurs in 10 pairs — the effective sample size is far below 60.
  - *Suggested fix:* Add several decoy proteins (different fold, different size) as negatives and report distributions (mean +/- SD, ROC/AUC) rather than a single difference. State the effective number of independent comparisons for the 13C set and prefer a per-compound or leave-one-out summary.
  - *Verifier:* The referee's mechanical facts check out. bench.py L44 `diff = read_bruker_2d(oaa / "103")` is a single spectrum; L52 `diff_s = fn(ref, diff)` yields one number (Table S2: LCC 2-vs-OAA = 0.1815), and both `separation = mean_same - diff_s` and `margin = min_same - diff_s` are governed by that one value. In bench_13c.py the 60 "different" pairs are indeed `itertools.combinations(names, 2)` filtered by compound (verified: C(12,2)=66, minus 6 same = 60), so each of the 12 spectra recurs in 10 different-pairs — the pairs are not independent. There are no error bars anywhere. So the statistical substance of the finding is real.

But the referee frames this as an undisclosed weakness, and that is where it is OVERSTATED. The manuscript itself states it plainly. Main text L135-137: "To probe the opposite regime — and to replace the single negative spectrum above with a large negative set — we used public HSQC peak lists from the simpleNMR example set for six compounds..." The author explicitly calls the dense benchmark's negative "the single negative spectrum" and introduces the sparse benchmark expressly to broaden the negative set. The whole design rationale (L118-119) is "two deliberately opposite regimes, so that any measure claiming general utility must succeed on both" — the headline claim is a cross-regime ranking (abstract: "0.75 and 0.81"), not a statistically-powered estimate from the dense benchmark alone. SI caveat (iii) (S10) further concedes "the between-compound contrast is strong (60 pairs) but the same-compound side is only 6 pairs," showing awareness of sample-size limits, and caveats (i)/(ii) restrict the claim to "the robust ranking, not the exact numerical values." The robustness sweep (S10, Σ=0.71–0.77) varies the blur, not the decoy, so it does NOT cure the n=1 decoy — the referee is right that the dense number rests on one spectrum. However, the sparse regime supplies 6 distinct decoy compounds, so "n=1 negative" is not true of the paper as a whole.

Net: the residual concern — one decoy protein, correlated sparse pairs, no error bars/CIs on a headline abstract number — is a legitimate limitation, but it is disclosed, partially mitigated by a second 6-compound benchmark, and framed as a ranking rather than a statistical estimate. It is a fair request for more decoys/error bars (MINOR), not a MAJOR flaw that the author cannot rebut, because the paper does not rest its discrimination claim on the single negative in the way the finding asserts.
- **[MAJOR] Flagship protein benchmark is not independently reproducible (Angewandte data policy)** — *verification: CONFIRMED → MAJOR*
  - *Location:* MS Data Availability Statement; bench.py L65-66; SI S12 Reproduction
  - The Data Availability Statement says the dense benchmark "uses PRL3 and OAA Bruker spectra available from the author on reasonable request", and bench.py hard-codes private local paths (L65-66 /Users/donghanlee/PRL3_mark_5nov24, /Users/donghanlee/Downloads/OAA_CEST_277K_04may15). The entire dense-protein regime — which is the paper's primary drug-discovery motivation and its most striking result (0.75 vs the next best 0.29) — therefore cannot be reproduced by a reader or referee. Only the 13C benchmark (auto-downloaded peak lists) is reproducible. Angewandte's data-availability policy expects the data underlying the central claims to be deposited, not gated behind an author request.
  - *Suggested fix:* Deposit the processed PRL3 and OAA spectra (or the windowed rendered matrices sufficient to rerun bench.py) in a public repository with a DOI. If the raw titration data are restricted, at minimum release the derived 2D matrices over the benchmark window so Table S2 is reproducible.
  - *Verifier:* Every factual claim in the finding checks out against the sources. MS Data Availability Statement (LCC_angewandte.md L223-224): "The dense $^1$H–$^{15}$N benchmark uses PRL3 and OAA Bruker spectra available from the author on reasonable request." bench.py L65-66 hard-codes private absolute paths: default=Path("/Users/donghanlee/PRL3_mark_5nov24") and default=Path("/Users/donghanlee/Downloads/OAA_CEST_277K_04may15"). The dense 1H-15N benchmark is indeed the flagship: the comparison table (L198, L202) shows LCC separation 0.75 vs the next method (Bin/Bodis) 0.29 on the dense protein column, and the abstract/intro (L37-48) frame the dense protein regime as the drug-discovery (fragment-/structure-based) motivation. Only the sparse 1H-13C benchmark auto-downloads its data (simpleNMR peak lists, L225-226 and SI S12), so it is the sole reproducible leg. SI S12 "Reproduction" lists bench.py but supplies no deposited spectra for it; the SI caveats hedge only the sparse benchmark's render-width dependence, never the dense benchmark's data availability. So the referee did not misread and the manuscript contains no hedge that neutralizes the concern. The core defect — the paper's single most striking quantitative result cannot be independently regenerated without an author request — is real and unrebuttable. Skeptical counterweight: the software, harnesses, tabulated scores, test suite, and the entire sparse benchmark ARE public, and the qualitative cross-regime claim ("LCC lowers different-class similarity in both regimes") survives on the open sparse leg alone; "available on reasonable request" is a statement Angewandte tolerates rather than an outright policy violation. These mitigations mean the referee's framing ("entire dense-protein regime cannot be reproduced," implying the paper's central claim is unverifiable) is somewhat stronger than warranted, and the remedy (deposit processed spectra on Zenodo) is trivial. Still a legitimate, standard reproducibility request touching the headline number.
- **[MAJOR] Scope and significance may be below the bar for a general-chemistry journal** — *verification: OVERSTATED → MINOR*
  - *Location:* MS abstract and throughout; MS "the spectral-imaging analogue of the cosine / contrast-angle similarity long used in mass-spectral library search"
  - The contribution is a similarity score for one specific data type (2D 1H-X HSQC), validated on one protein titration against one decoy and six rendered small molecules. The novelty is essentially a mean-centred normalized cross-correlation (Pearson/ZNCC) of Gaussian-blurred images — a standard image-similarity operation (the paper itself calls it "the spectral-imaging analogue of the cosine / contrast-angle similarity"). This is solid, useful methods work but reads as a specialized chemoinformatics/NMR-methods contribution rather than a broad advance for a general Angewandte readership; the demonstrated impact (a better number on two small benchmarks) is narrow.
  - *Suggested fix:* Either broaden the demonstration to a scale that shows general impact (e.g. a real dereplication search over a public library such as BMRB/NMRShiftDB with retrieval statistics, and a real multi-protein CSP screen), or redirect to a specialized methods journal (J. Cheminform., Anal. Chem., Magn. Reson. Chem.).
  - *Verifier:* The referee's factual basis is accurate but not misread: the paper is a single-data-type method (2D 1H–X HSQC; title/abstract), Benchmark 1 uses one protein titration vs. one decoy protein (OAA, line 128; Table S2), Benchmark 2 uses six peak-list-rendered compounds (S9), and the paper does call the score "the spectral-imaging analogue of the cosine / contrast-angle similarity long used in mass-spectral library search" (line 98). However, framing this as a MAJOR defect overstates it. (1) It is a rebuttable editorial fit judgment, not a scientific flaw, and the manuscript supplies the breadth material the referee omits: LCC is claimed "directly useful wherever a robust 2D-HSQC similarity is needed: automated titration tracking and CSP read-out, library and dereplication searches, and quality control of repeated measurements" (lines 168–170), motivated across drug discovery, metabolomics, natural-product dereplication and BMRB/NMRShiftDB/HMDB library search (lines 47–56). (2) The "novelty is just standard ZNCC" framing undersells the actual claim: the paper never claims to invent Pearson/ZNCC (line 98 explicitly names it a known analogue); S11 Table S6 credits the three predecessors for the physical ingredients and lists "mean-centred zero-lag correlation (the discriminating step)" as the new contribution, and the S10 ablation demonstrates it (removing mean-centring collapses separation 0.75→0.59). (3) Validation is two deliberately opposite regimes plus robustness/blur sweeps (Table S5) and a 19-test property suite, with LCC the only method winning both (Table 1). The genuine empirical thinness (single decoy protein; only 6 same-compound pairs) is real but already disclosed as caveats (S10 caveat iii; lines 155–161) and the second benchmark was added precisely "to replace the single negative spectrum above with a large negative set" (lines 135–136). The observation is real and text-grounded but is a subjective, partly pre-empted editorial concern, not a MAJOR defect the author cannot rebut.
- **[MINOR] "Faithful reimplementations" overstates fidelity to the original methods**
  - *Location:* MS ("our own faithful reimplementations"); hsqc_methods.py L212-215, L294-327; MS/SI Table 1 row "Bin + 45 rotation"
  - The reimplemented baselines deviate from the originals in ways that affect discrimination. nn_peak_similarity uses a plain symmetric mean nearest-neighbour distance mapped by 1/(1+d) (hsqc_methods.py L294-327), dropping Pierens et al.'s genetic-algorithm optimal assignment entirely. quadtree_similarity adds a self-normalization s(xy)/sqrt(s(xx)s(yy)) that the code comment admits is not in the original ("The raw score is not 1 for identical spectra, so it is normalized...", L212-215). "Bin + 45 rotation" (Table 1) is an author-introduced variant, not a literature method. Calling these "our own faithful reimplementations" (MS) is too strong.
  - *Suggested fix:* State explicitly where each reimplementation departs from the source (NN assignment, tree normalization) and label the 45-degree bin as a variant of this work. Soften "faithful" to "representative reimplementations."
- **[MINOR] "Set by spectroscopy, not tuning" is undercut by the per-regime and grid changes**
  - *Location:* MS abstract/conclusion; SI Table S1 ("score grid-invariant"); SI S10; bench_13c.py L49
  - The paper repeatedly claims both knobs are "set by spectroscopy rather than by tuning" (abstract; MS conclusion; Table S1). But between the two benchmarks not only the physical sigmas change (0.03/0.30 -> 0.05/0.5) but also the render step sizes (bench_13c.py L49 step_f2=0.02, step_f1=0.2 vs the 0.01/0.10 defaults), and the step is not a spectroscopic quantity. The score is asserted "grid-invariant" (Table S1) yet S10 reports it does change with grid ("still 0.36 ... when the grid is coarsened"), so the invariance is approximate. This weakens the "no tuning" framing.
  - *Suggested fix:* Report the score's actual grid dependence and either fix the step relative to sigma by a stated rule (e.g. always sigma/3) or acknowledge the grid as a third, non-spectroscopic parameter.

### Referee 4 — MINOR REVISION

*Focus: Manuscript presentation and Angewandte Chem. Int. Ed. format compliance: abstract/keywords/TOC norms, reference style, figure existence and distinctness, caption self-containment, in-text claim traceability, terminology consistency, and wording-level over-claim.*

**Assessment.** On the presentation axis this manuscript is in good shape and largely Angewandte-compliant: it carries 5 keywords, a TOC entry, a correctly formatted reference list ([n], journal italic, year bold, volume italic, page range), a 221-word abstract, and self-contained figure captions. I specifically checked the flagged concern and can confirm Figures 1 and 2 both exist and are genuinely distinct images (../results/lcc_schematic.png is the render/blur/score schematic; ../results/comparison_all.png is the five-method separation bar chart), and the bar-chart values match Table 1 exactly. All headline numbers in the main text trace correctly to Tables 1, S2, S3 and S4. The remaining issues are compliance and consistency items: the headline dense-protein benchmark rests on data only available "on reasonable request" rather than deposited (a real reproducibility/FAIR-compliance gap), the graphical abstract is described rather than supplied as a file, and there are minor terminology inconsistencies (hyphenation of the method's own name; an undefined "ZNCC" that appears only inside Figure 1; two different multipliers quoted for the same headline comparison). None invalidate a result, so I recommend minor revision.

**Findings:**

- **[MAJOR] Headline benchmark data not deposited — only 'on reasonable request'** — *verification: OVERSTATED → MINOR*
  - *Location:* MS 'Data Availability Statement' (L223-224); headline values Table 1 (L202) and Table S2 (L238-240)
  - The Data Availability Statement states: 'The dense 1H–15N benchmark uses PRL3 and OAA Bruker spectra available from the author on reasonable request.' The paper's flagship result — the 0.75 separation / 0.71 margin that anchors the abstract and Table 1 — derives entirely from these PRL3/OAA spectra. Although the numbers are internally traceable to Table S2, an independent reader cannot reproduce the central protein claim from any public artefact; only the sparse 1H–13C (simpleNMR) benchmark is fully reproducible. This falls short of Angewandte/Wiley FAIR data expectations for the datum that carries the paper's principal claim.
  - *Suggested fix:* Deposit the processed PRL3 and OAA HSQC spectra (or a derived, redistributable rendered-grid form) in a public repository (e.g. Zenodo/BMRB) with a DOI, or add a small synthetic dense-protein surrogate that reproduces the qualitative ranking from public data. At minimum, state the licensing/consent barrier explicitly.
  - *Verifier:* Every factual claim in the finding is verifiable and the manuscript does NOT rebut it — this is not a misread. Data Availability Statement (L223-224) reads verbatim: "The dense $^1$H–$^{15}$N benchmark uses PRL3 and OAA Bruker spectra available from the author on reasonable request." The headline 0.75 / 0.71 in Table 1 (L202: "**LCC (this work)** | **0.75** | **0.71** ...") and the abstract (L37: "0.75 and 0.81"; L39 "collapses the protein separation from 0.75") both anchor on this dense benchmark. Table S2 (SI L228-240) makes the numbers internally traceable (2 vs OAA = 0.1815; sep 0.75; margin 0.71). The dense harness bench.py confirms the raw data is not in the repo: it reads Bruker spectra from external absolute paths (`default=Path("/Users/donghanlee/PRL3_mark_5nov24")` and `.../OAA_CEST_277K_04may15`), with the docstring stating "Data lives outside the repo." The sparse 1H-13C benchmark, by contrast, "uses public HSQC peak lists from the simpleNMR example set ... downloaded automatically by the benchmark script" — genuinely reproducible. So the referee is factually correct: the flagship dense-protein datum cannot be reproduced from any public artefact.

However, the MAJOR severity is overstated. (1) "Available on reasonable request" is an explicitly recognized/permitted tier under Wiley/Angewandte data policy for raw experimental data such as Bruker spectra — it is a weaker option, not a policy violation, so "falls short of FAIR expectations" is stronger than warranted. (2) All derived scores are public and fully traceable in Table S2, and the scoring harness (bench.py) plus all four methods are public under MIT — an independent reader can audit exactly how the numbers were produced. (3) The core scientific claim (LCC separates same/different better than prior methods) is independently reproducible on the fully-public simpleNMR sparse benchmark, which also shows LCC winning (0.81 sep, headline alongside 0.75). The specific dense raw spectra are request-only, but the paper's methodological conclusion is not hostage to them.
- **[MINOR] Graphical abstract (TOC graphic) described but not supplied**
  - *Location:* MS 'Table of Contents Entry' (L313-314)
  - The TOC section ends with '*Suggested TOC graphic: the separation bar chart (Figure 2) with the three-step LCC schematic inset.*' Angewandte requires an actual graphical-abstract image submitted with the TOC text, not a description of one. As it stands there is no TOC graphic file.
  - *Suggested fix:* Produce and submit the composite TOC graphic (bar chart + schematic inset) as an image asset and remove the italic 'suggested' placeholder line.
- **[MINOR] Inconsistent hyphenation of the method's own name**
  - *Location:* MS title (L2) vs abstract (L30), L74; SI title (L2)
  - The manuscript title and SI title read 'Lineshape-Correlation Coefficient' (hyphenated), but the body uses the unhyphenated 'Lineshape Correlation Coefficient' three times (e.g. abstract L30 '...introduce the Lineshape Correlation Coefficient (LCC)', and L74). A method's proper name should be spelled identically throughout.
  - *Suggested fix:* Pick one form (recommend unhyphenated 'Lineshape Correlation Coefficient' to match the LCC expansion) and apply it consistently in both title and body of MS and SI.
- **[MINOR] 'ZNCC' acronym appears only inside Figure 1 and is never defined in text**
  - *Location:* Figure 1 image (results/lcc_schematic.png), panel 3 label; MS caption L176-183 does not expand ZNCC
  - Figure 1a panel 3 is labelled 'Score (mean-centred ZNCC)'. The acronym ZNCC (zero-lag normalized cross-correlation) does not occur anywhere in the manuscript or SI text — the running text always writes it out. A reader meets an undefined acronym in the figure. Figure captions/labels must be self-contained.
  - *Suggested fix:* Expand 'ZNCC' at first use in the Figure 1 caption, e.g. 'mean-centred zero-lag normalized cross-correlation (ZNCC)', so the label is self-explanatory.
- **[MINOR] Two different multipliers quoted for the same headline comparison**
  - *Location:* MS L129-130 vs SI S8 L242
  - The main text states LCC's protein advantage as 'a margin of 0.71 — three times the bin method's 0.24' (L129-130, margin basis, ~3x), while SI S8 states 'LCC separates same-protein from different-protein spectra ≈ 2.6× better than the previous best' (L242, separation basis, 0.75/0.29). Both are arithmetically correct on their respective bases, but a reader cross-referencing MS and SI sees '3×' and '2.6×' for what looks like the same claim without the basis being flagged.
  - *Suggested fix:* State the basis explicitly at each mention ('3x on worst-case margin', '2.6x on separation') or unify to a single quoted multiplier.
- **[MINOR] Headline metric alternates between 'separation' and 'margin' without a single lead figure**
  - *Location:* MS abstract L37-38; results L128-130, L139; definitions L120-122
  - The abstract leads with separation ('0.75 and 0.81 on a [0,1] scale', L37) and margin ('margin <= 0.04', L38), while the protein results paragraph leads with '0.94 ... margin of 0.71' (L128-130) and the small-molecule paragraph leads with 'separation 0.81' (L139). Both metrics are defined (L120-122), but the shifting lead metric adds friction. This is presentation only; the numbers are all correct and traceable to Table 1 / Tables S2-S3.
  - *Suggested fix:* Foreground one primary metric (separation, since it is what the abstract and Figure 2/Table 1 emphasize) consistently across abstract and both result paragraphs, quoting margin as the secondary/threshold metric.
