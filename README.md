# Spectrum similarity for Bruker 1D NMR

This tool implements the bin-based spectrum similarity measure described in:

Lorant Bodis, Alfred Ross, Erno Pretsch, "A novel spectra similarity measure",
Chemometrics and Intelligent Laboratory Systems 85 (2007) 1-8.

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
