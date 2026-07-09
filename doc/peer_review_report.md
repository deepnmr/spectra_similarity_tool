# Peer-review report — LCC Application Note (simulated)

Simulated Oxford Bioinformatics review: three independent referees (NMR spectroscopy;
methods/reproducibility; applications/novelty) plus an editor decision. Each referee read
the manuscript and **inspected the repository** (ran the tests, checked the result JSONs,
grepped the code), so most findings are code/data-verified rather than opinion.

## Decision: **MAJOR REVISION** (unanimous)

| Referee | Focus | Recommendation |
| --- | --- | --- |
| R1 | Biomolecular NMR spectroscopy | major revision |
| R2 | Methods / reproducibility | major revision |
| R3 | Applications / novelty | major revision |
| Editor | — | **major revision** |

**Not in question:** the method is physically sound (the editor and R1 numerically reproduced
the shift-tolerance law `exp(-d²/4σ²)`, 0.778 vs 0.779), every number in Table 1 regenerates
from the committed JSONs, the 19-test suite passes, and the software is small, MIT-licensed and
dependency-light — a good fit for the Application Note format.

**What blocks acceptance:** the central quantitative claims are not backed by deposited, runnable,
honestly-described computation, and two are outright defects in the current artifacts.

## Consolidated concerns (by severity)

Status: **V** = verified against code/data by a referee · **A** = analytical/expert argument.

| # | Sev | Theme | Raised by | Status | Issue |
| --- | --- | --- | --- | --- | --- |
| 1 | MAJOR | Missing ablation code | R2, R3, Ed | V | The headline mean-centring ablation (0.75→0.59) is computed by **no code**; `0.59` appears only in the prose. Grep of the repo confirms it. |
| 2 | MAJOR | 13C positive class mislabeled | R2, Ed | **V (bug)** | The olivetol pair is **byte-identical** (LCC=1.000 by construction, not measurement); rotenone/santonin/indanone are `-bckup` files, not independent recordings. Only **menthol & chartreusin** are genuine repeats — and LCC is weakest there (0.58). "Six compounds each recorded twice" overstates the data. |
| 3 | MAJOR | Edited-HSQC negatives (clip vs abs) | R1 | **V (bug)** | The shipped default `render_image` uses `baseline='clip'` (zeroes negatives) → deletes **CH2 crosspeaks** in multiplicity-edited HSQC. But `bench_13c.py` uses `abs(intensity)`. The validated path ≠ the shipped path; on a real edited 2rr, ~⅓ of crosspeaks vanish. |
| 4 | MAJOR | 1H-15N single negative + no stats | R1, R2, R3, Ed | V | Discrimination (0.75, margin 0.71, "3×") rests on **one** different-protein spectrum (OAA). No variance, no error bar, no significance test anywhere in the code. |
| 5 | MAJOR | 1H-15N not reader-reproducible | R1, R2, R3, Ed | V | `bench.py` hardcodes non-public absolute paths; data "on request". Table 1's 15N column and Fig 1's blue bars cannot be regenerated. |
| 6 | MAJOR | Baseline fairness + missing baseline | R2, R3, Ed | V/A | Baselines are the authors' reimplementations; NN's `1/(1+d)` over a 0–165 ppm window may manufacture the ~0.9–1.0 "saturation". No comparison to **plain cosine/Pearson of the matrices** or a **no-blur (σ→0)** variant, so the blur's contribution and novelty over standard 2D correlation are unproven. |
| 7 | MAJOR | Global shift vs local alignment | R1 | A | "No alignment" conflates a rigid **global** referencing offset (an artifact, correctable, cannot cause a false match) with local per-peak matching. Real 15N referencing offsets 0.1–0.3 ppm erode scores unrelated to identity. |
| 8 | MAJOR | Intensity reproducibility | R1 | A | LCC is intensity-weighted, but HSQC intensities are far less reproducible than positions (relaxation, concentration, transfer efficiency). Undermines the cross-condition library-matching use case; visible as min_same 0.39/0.58. |
| 9 | MINOR | Non-independent 13C pairs | R2, R3 | V | 60 "different" pairs come from 12 spectra (each in ~10 pairs) → effective N far below 60; needs compound-level (clustered) resampling, not per-pair means. |
| 10 | MINOR | Render circularity | R3 | A | 13C test spectra are Gaussians blurred by a Gaussian — LCC's own forward model. Add a Lorentzian/noise robustness check or a real 2rr pair. |
| 11 | MINOR | Timing not reproduced | R3 | V | "~6 ms/pair" measured ~30 ms on a 250×350 grid; 13C grids slower. State input resolution and cost vs grid size. |
| 12 | MINOR | Clamp hides anti-correlation | R1, R3 | V | `[0,1]` clamp discards genuinely negative correlations; report raw different-class values or justify. |
| 13 | MINOR | Lineshape terminology | R1 | A | Kernel is called "physical linewidth" but σ is a variance not FWHM and σ_15N=0.30 ppm is drift-dominated, not a real linewidth. |
| 14 | MINOR | Dependency pinning / CI / coverage | R3 | V | No requirements.txt/pyproject; only 3 of 19 tests exercise LCC; add version pin + a benchmark-regression test. |
| 15 | MINOR | Ablation only on dense regime | R1 | A | Mean-centring is "the discriminating step" only in the dense image; in the ~99%-zero sparse image Pearson≈cosine. Scope the claim or measure it for 13C. |
| 16 | MINOR | 13C window / px-per-σ | R1 | V | 0–165 ppm 13C truncates aldehyde carbons; benchmark runs at 2.5 px/σ, not the stated ≥3. |

## Prioritized revision checklist (editor)

1. **Ablation harness** — a `--no-center` switch (plus a `σ→0` no-blur variant and a plain cosine-of-matrices baseline) that regenerates 0.59 and writes it to a results JSON, for **both** regimes; cite it.
2. **More 1H-15N negatives** — several unrelated proteins of varying size/dispersion; report per-class spread (bootstrap CI) and a rank test (Mann–Whitney). If unavailable, drop the "widest separation"/"3×" wording for this regime and present it as one illustrative case.
3. **Deposit the dense data** — PRL3/OAA (or a derivative sufficient to reproduce the scores) with a DOI; make `bench.py` fetch it like `bench_13c.py`, or demote the regime from a benchmark.
4. **Relabel the 13C pairs truthfully** — drop or flag the byte-identical olivetol pair; report same-class stats on the genuine-repeat subset (menthol, chartreusin) separately from the backup pairs.
5. **Validate/justify baselines** — reproduce at least one published value from Bodis/Castillo/Pierens, state the `1/(1+d)` map and normalization, and soften "saturate/overlapping" to "in our reimplementation, with these parameters".
6. **Statistics respecting non-independence** — compound-level bootstrap or paired per-compound comparison for the 13C set; show per-pair distributions, not only means.
7. **Robustness to the render model** — Lorentzian render or peak-height noise, or a real experimental 2rr pair, to show the ranking survives.
8. **Report raw (unclamped) different-class correlations** or justify the `[0,1]` clamp.
9. **Fix timing** — native input resolution behind "~6 ms" and cost vs grid size.
10. **Pin dependencies + minimal CI** — requirements/pyproject with tested Python/NumPy, a symmetry assertion, and a benchmark-regression test.
11. *(optional, raises impact)* a small library-retrieval demo on the simpleNMR set.

## Genuine defects to fix now (independent of the paper's fate)

- **#2** olivetol is byte-identical and the "-bckup" pairs are not independent recordings — the `1H-13C` positive class and its description are partly wrong.
- **#3** clip-vs-abs: the shipped default silently deletes CH2 peaks on edited HSQC.
- **#1** the mean-centring ablation number has no committed code.
- **#11** the "~6 ms" timing does not reproduce (~30 ms).

The rest are strengthening/validity requests, not errors.
