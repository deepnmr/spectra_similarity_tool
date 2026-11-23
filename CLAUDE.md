# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Interpreter

Use `uv run` for everything; `pyproject.toml` + `uv.lock` pin the environment (numpy; pytest and matplotlib in the `dev` group, matplotlib only for `--plot` and the figure scripts). Bare `python3`/`python3.11` on this machine have no numpy. `pyproject.toml` is environment-only (`package = false`): the layout stays flat modules, nothing is installed.

## Commands

```bash
uv run pytest                                         # all tests (run from repo root; tests import the top-level modules directly)
uv run pytest tests/test_hsqc_lcc.py -k monotonic   # one file / one test
uv run python hsqc_lcc.py EXP1 EXP2 --f2-min 6.5 --f2-max 10 --f1-min 105 --f1-max 130 --json
uv run python hsqc_methods.py EXP1 EXP2 --method nn      # quadtree | nn
uv run python hsqc_similarity.py EXP1 EXP2               # 2D bin method
uv run python spectrum_similarity.py EXP1 EXP2           # 1D bin method
```

Benchmarks (each writes JSON under `results/`, which the figure scripts and tests read):

```bash
uv run python bench_13c.py          # sparse 1H-13C, public data auto-downloaded into data_13c/ (gitignored); ~15 s
uv run python bench_retrieval.py    # retrieval/bootstrap stats on the same set -> results/retrieval_13c.json
uv run python bench_nhsqc.py        # dense 1H-15N, 23 titration points + 2 decoys, needs ./Nhsqc (gitignored, not redistributed)
uv run python bench.py --prl3 DIR --oaa DIR   # legacy single-decoy dense benchmark; or set $PRL3_DIR/$OAA_DIR
for f in make_plots make_fig3 make_si_figs; do uv run python results/$f.py; done   # regenerate figures from the JSON
```

The bin method dominates benchmark time (one `histogram2d` per grid size, up to 100 per pair); `hsqc_similarity` drops zero-weight points before that loop, which is bit-exact and is what keeps `bench_13c.py` at seconds instead of ~10 min. The dense-benchmark test in `tests/test_bench_nhsqc.py` runs the full `bench_nhsqc.run` when `Nhsqc/` is present and skips otherwise.

## Architecture

Flat module layout, no package. Dependency order:

```
spectrum_similarity.py   1D bin method + JCAMP `procs` parser (parse_jcamp, _number, upper_envelope)
   └─ hsqc_similarity.py 2D bin method + Bruker `2rr` reader -> Spectrum2D(ppm_f2, ppm_f1, intensity, source)
         ├─ hsqc_methods.py  Castillo quadtree, Pierens NN; owns the shared _window()/_overlap_ranges() preprocessing
         └─ hsqc_lcc.py      STCC (CLI name `lcc`), un-centred cosine ablation, ECHOS baseline (`echos`: box bin + Pearson, no blur), experimental local-contrast
bench*.py                harnesses; bench_nhsqc imports bench._methods, bench_retrieval imports bench_13c.METHODS
```

- `Spectrum2D` is the only data contract. `intensity[f1_index, f2_index]`, ppm axes descend. Bench scripts synthesize it directly from peak lists, so any method that accepts `Spectrum2D` works on both Bruker data and rasterized sticks.
- `_baselined_window` in `hsqc_similarity.py` is the single baseline + ppm-window step; `prepare_spectrum_2d`, `hsqc_methods._window` and the contour plot all go through it.
- Every method returns a dict with at least `method`, `similarity`, `range_f2`, `range_f1`, `source_x`, `source_y`; benches only consume `["similarity"]`. Self-similarity must be exactly 1.0 (`bench_nhsqc.run` asserts it).
- Ranges: `None` means "common ppm overlap of the two spectra", resolved once by `_overlap_range` in `hsqc_similarity.py` (reused by `hsqc_methods._overlap_ranges`). CLIs share `_paired_range` for `--f2-min/--f2-max` pairs.
- STCC pipeline in `hsqc_lcc.py`: `_window` (baseline + area-weight + normalize) -> `render_image` (histogram2d onto a fixed ppm grid, separable Gaussian blur, `RenderedImage` carries step and sigma) -> `_zncc` (mean-centred, zero-lag, clamped). `_render_pair`/`_report` are the shared scaffolding; a new render-based method needs only a feature transform and a `_report` call. Zero lag and no alignment are deliberate design choices documented in the module docstring; do not add a shift search.
- Local-contrast validates its own inputs (sigma/step must be finite and > 0) because its DoG background needs a real kernel. LCC allows sigma 0 as "no blur"; keep that contract.
- Method registries are hand-maintained dicts: `bench._methods()` (dense) and `bench_13c.METHODS` (sparse, different sigma/step because the stick raster is coarser). Adding a method means adding it to both; `tests/test_bench_*.py` assert `local_contrast` is registered.

## Data and results

- Raw Bruker data (`Nhsqc/`, PRL3/OAA dirs) and `data_13c/` are gitignored. The derived JSON in `results/` is tracked and is the source of every number in `README.md`, `results/*.md`, and `doc/LCC_angewandte*.md`. If a method or benchmark changes, regenerate the JSON, the figures, and then the prose; a refactor must leave `results/comparison_13c.json` value-for-value reproducible (the 13C bench was re-run after the last refactor and matched).
- `doc/` holds the manuscript (`LCC_angewandte.md` + `_SI.md`, rendered docx/pdf). Two reference PDFs there are copyrighted and gitignored; only the Pierens open-access PDF is committed.
- Naming: the paper calls the method STCC; code and CLI still say `lcc`/`LCC` (Lineshape Correlation Coefficient). `lcc` stays the default CLI method; `local-contrast` is experimental until it improves held-out metrics.
