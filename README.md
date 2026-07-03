# Spectrum similarity for Bruker NMR

This tool implements the bin-based spectrum similarity measure described in:

Lorant Bodis, Alfred Ross, Erno Pretsch, "A novel spectra similarity measure",
Chemometrics and Intelligent Laboratory Systems 85 (2007) 1-8.

`spectrum_similarity.py` handles 1D spectra; `hsqc_similarity.py` extends the
same method to 2D (HSQC and other processed 2D experiments).

## 1D

The input is a processed Bruker 1D spectrum. Pass either an experiment directory,
a `pdata` directory, or a processed `1r` file. The program reads `procs` and `1r`.

## Usage

```bash
python3 spectrum_similarity.py /path/to/exp1 /path/to/exp2
```

Useful options:

```bash
python3 spectrum_similarity.py exp1 exp2 --procno 1 --ppm-min 0 --ppm-max 10
python3 spectrum_similarity.py exp1 exp2 --min-bin-width 0.4 --json
python3 spectrum_similarity.py exp1 exp2 --norm-x 19 --norm-y 19
python3 spectrum_similarity.py exp1 exp2 --plot result.png
python3 spectrum_similarity.py exp1 exp2 --plot --show
```

Defaults:

- `--min-bin-width 0.4`: the value recommended in the paper for 1H NMR.
- `--baseline clip`: negative intensities are set to zero before integration.
- If no ppm range is supplied, the common ppm overlap of the two spectra is used.
- Each spectrum is normalized to total integral 1 unless `--norm-x` and
  `--norm-y` are supplied. For predicted or assigned proton spectra, set these
  to the proton counts when desired.
- `--plot` saves a graph with the two normalized spectra and the `SI_n` /
  `SI*_n` similarity curves. Without a path it writes `spectrum_similarity.png`.
- `--show` displays the same graph interactively when a GUI backend is available.

## 2D HSQC

`hsqc_similarity.py` applies the same bin-method concept to processed Bruker 2D
spectra (HSQC, HMBC, COSY, …). Each spectrum is subdivided into an `n` by `n`
grid of bins whose widths shrink as `n` grows, and the per-grid `SI_n` index and
its `SI*_n` envelope are averaged into one score in `[0, 1]`.

```bash
python3 hsqc_similarity.py /path/to/exp1 /path/to/exp2
python3 hsqc_similarity.py exp1 exp2 --f2-min 0 --f2-max 10 --f1-min 0 --f1-max 160
python3 hsqc_similarity.py exp1 exp2 --min-bin-width-f2 0.4 --min-bin-width-f1 4.0
python3 hsqc_similarity.py exp1 exp2 --plot result.png
python3 hsqc_similarity.py exp1 exp2 --json
```

The reader targets processed 2D Bruker data and reconstructs the submatrix
(tile) layout:

```text
experiment/
  pdata/
    1/
      2rr
      procs    # direct dimension F2 (e.g. 1H)
      proc2s   # indirect dimension F1 (e.g. 13C or 15N)
```

It reads `SI`, `OFFSET`, `SW_p`, `SF`, and `XDIM` from `procs`/`proc2s`, plus
`BYTORDP`, `DTYPP`, and `NC_proc` from `procs`, to decode `2rr` and build both
ppm axes.

Defaults:

- `--min-bin-width-f2 0.4`: the 1H value recommended in the paper.
- `--min-bin-width-f1 4.0`: a wider bin for the 13C/15N axis, whose ppm range is
  roughly ten times larger.
- If no ppm range is supplied for a dimension, the common overlap of the two
  spectra is used in that dimension.

## Bruker assumptions

The reader targets processed 1D Bruker data:

```text
experiment/
  pdata/
    1/
      1r
      procs
```

It uses `SI`, `OFFSET`, `SW_p`, `SF`, `BYTORDP`, `DTYPP`, and `NC_proc` from
`procs` to decode the data and build the ppm axis.

## Test

```bash
python3 -m pytest
```
