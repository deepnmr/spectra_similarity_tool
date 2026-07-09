#!/usr/bin/env python3
"""1H-13C HSQC cross-regime benchmark on SPARSE small-molecule spectra.

Complements bench.py (dense protein 1H-15N) with the sparse small-molecule 1H-13C
regime the tree/NN methods were designed for. Data: public HSQC peak lists from the
simpleNMR example set (EricHughesABC/simpleNMR) — six compounds each recorded twice
(a variant pair, e.g. menthol in two solvents), giving 12 spectra. A good measure
scores the 6 same-compound pairs HIGH and the 60 different-compound pairs LOW.

The peak lists are rasterized to a synthetic 2D spectrum (one Gaussian per peak) so
the existing Spectrum2D methods run unchanged. Run with python3.11 (numpy/matplotlib).
The peak-list data are downloaded on first run (not redistributed here).
"""
from __future__ import annotations
import itertools
import json
import urllib.request
from pathlib import Path

import numpy as np

from hsqc_similarity import Spectrum2D, hsqc_similarity
from hsqc_methods import nn_peak_similarity, quadtree_similarity
from hsqc_lcc import lcc_similarity, cosine_similarity

RAW = "https://raw.githubusercontent.com/EricHughesABC/simpleNMR/main/exampleProblems"
# local_name -> (folder, remote_filename); two backups store the json under the compound name.
# NOTE: the olivetol pair (Olivetol/Olivetol_A) was dropped -- its two files are BYTE-IDENTICAL
# (same MD5, identical peak lists), so it was a self-comparison scoring 1.00 for every method,
# an information-free positive that inflated mean_same. run() also guards against any such
# duplicate at load time (see _assert_distinct). The remaining 5 pairs are genuinely distinct
# recordings/re-picks (verified: different peak coordinates).
PAIRS = [
    ("menthol",     ("menthol_CDCl3", "menthol_CDCl3.json"), ("menthol_eeh", "menthol_eeh.json")),
    ("rotenone",    ("Rotenone", "Rotenone.json"), ("Rotenone-bckup", "Rotenone.json")),
    ("santonin",    ("santonin", "santonin.json"), ("santonin_bckup", "santonin.json")),
    ("chartreusin", ("Chartreusin", "Chartreusin.json"), ("Chartreusin_eh", "Chartreusin_eh.json")),
    ("indanone",    ("2-ethyl-1-indanone", "2-ethyl-1-indanone.json"),
                    ("2-ethyl-1-indanone-bckup", "2-ethyl-1-indanone-bckup.json")),
]

# Common grid over small-molecule 1H-13C HSQC. 1H direct (F2), 13C indirect (F1).
H_LO, H_HI, H_STEP = 0.0, 10.0, 0.01
C_LO, C_HI, C_STEP = 0.0, 165.0, 0.10
SIG_H, SIG_C = 0.03, 0.5  # render linewidth in ppm
WIN = dict(range_f2=(H_LO, H_HI), range_f1=(C_LO, C_HI))

METHODS = {
    "bin_Bodis09":  lambda x, y: hsqc_similarity(x, y, min_bin_width_f2=0.1, min_bin_width_f1=1.0, **WIN)["similarity"],
    "bin_rot45":    lambda x, y: hsqc_similarity(x, y, min_bin_width_f2=0.1, min_bin_width_f1=1.0, rotate_deg=45.0, **WIN)["similarity"],
    "tree_Castillo13": lambda x, y: quadtree_similarity(x, y, **WIN)["similarity"],
    "nn_Pierens12": lambda x, y: nn_peak_similarity(x, y, **WIN)["similarity"],
    "lcc_new":      lambda x, y: lcc_similarity(x, y, sigma_f2=0.05, sigma_f1=0.5, step_f2=0.02, step_f1=0.2, **WIN)["similarity"],
    # Ablation baseline: same render/blur as LCC but WITHOUT mean-centring (un-centred cosine /
    # contrast angle). Isolates what mean-centring buys -- see the ablation in the SI.
    "cosine_uncentred": lambda x, y: cosine_similarity(x, y, sigma_f2=0.05, sigma_f1=0.5, step_f2=0.02, step_f1=0.2, **WIN)["similarity"],
}


def download(data_dir: Path) -> None:
    data_dir.mkdir(parents=True, exist_ok=True)
    for _, va, vb in PAIRS:
        for local, (folder, remote) in ((va[0], va), (vb[0], vb)):
            dest = data_dir / f"{local}.json"
            if dest.exists():
                continue
            url = f"{RAW}/{folder}/{remote}"
            urllib.request.urlretrieve(url, dest)
            print(f"downloaded {local}.json")


def load_peaklist(path: Path) -> Spectrum2D:
    d = json.load(open(path))
    hs = d["hsqc"]
    keys = list(hs["f2_ppm"].keys())
    stype = hs.get("signaltype", {})
    ppm_h = np.arange(H_LO, H_HI, H_STEP)
    ppm_c = np.arange(C_LO, C_HI, C_STEP)
    img = np.zeros((ppm_c.size, ppm_h.size))
    for k in keys:
        if stype.get(k, "Compound") != "Compound":
            continue
        h, c = float(hs["f2_ppm"][k]), float(hs["f1_ppm"][k])
        inten = abs(float(hs["intensity"][k]))  # sign encodes CH2 multiplicity, not absence
        if H_LO <= h <= H_HI and C_LO <= c <= C_HI:
            img[int(round((c - C_LO) / C_STEP)), int(round((h - H_LO) / H_STEP))] += inten
    img = _blur(_blur(img, SIG_C / C_STEP, 0), SIG_H / H_STEP, 1)
    return Spectrum2D(ppm_f2=ppm_h, ppm_f1=ppm_c, intensity=img, source=path)


def _blur(m, sig_px, axis):
    r = max(1, int(round(3 * sig_px)))
    x = np.arange(-r, r + 1)
    k = np.exp(-0.5 * (x / sig_px) ** 2)
    k /= k.sum()
    return np.apply_along_axis(lambda row: np.convolve(row, k, mode="same"), axis, m)


def _assert_distinct(specs: dict, same: list) -> None:
    """A same-compound pair must be two DISTINCT recordings, not the same file twice.
    Byte-identical inputs render to identical images and score 1.00 for every method, a
    fake positive that inflates mean_same (this caught the olivetol duplicate)."""
    for a, b in same:
        if np.array_equal(specs[a].intensity, specs[b].intensity):
            raise ValueError(
                f"same-compound pair ({a}, {b}) renders to an identical image -- "
                f"the two source files are duplicates, not independent recordings"
            )


def run(data_dir: Path) -> dict:
    download(data_dir)
    names = [v for _, a, b in PAIRS for v in (a[0], b[0])]
    comp = {v: c for c, a, b in PAIRS for v in (a[0], b[0])}
    specs = {n: load_peaklist(data_dir / f"{n}.json") for n in names}

    same = [(a[0], b[0]) for _, a, b in PAIRS]
    _assert_distinct(specs, same)
    diff = [(x, y) for x, y in itertools.combinations(names, 2) if comp[x] != comp[y]]

    rows, per_pair = {}, {}
    for m, fn in METHODS.items():
        s = [fn(specs[a], specs[b]) for a, b in same]
        dvals = [fn(specs[x], specs[y]) for x, y in diff]
        ms, md = float(np.mean(s)), float(np.mean(dvals))
        rows[m] = dict(self=float(fn(specs[names[0]], specs[names[0]])),
                       mean_same=ms, min_same=float(min(s)),
                       mean_diff=md, max_diff=float(max(dvals)),
                       separation=ms - md, margin=float(min(s) - max(dvals)))
        per_pair[m] = {c: float(v) for (c, _, _), v in zip(PAIRS, s)}
    return {"rows": rows, "per_pair": per_pair, "n_same": len(same), "n_diff": len(diff)}


def main() -> int:
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data-dir", type=Path, default=Path("data_13c"))
    p.add_argument("--out", type=Path, default=Path("results/comparison_13c.json"))
    args = p.parse_args()

    result = run(args.data_dir)
    rows = result["rows"]
    cols = ["self", "mean_same", "min_same", "mean_diff", "max_diff", "separation", "margin"]
    print(f"{'method':18} " + " ".join(f"{c:>10}" for c in cols))
    for m in sorted(rows, key=lambda n: -rows[n]["separation"]):
        r = rows[m]
        print(f"{m:18} " + " ".join(f"{r[c]:10.4f}" for c in cols))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    json.dump(result, open(args.out, "w"), indent=2)
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
