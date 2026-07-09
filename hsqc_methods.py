#!/usr/bin/env python3
"""Alternative HSQC similarity methods from the literature, for comparison with the
bin method in ``hsqc_similarity.py``.

Implemented here:

* ``quadtree_similarity`` -- the tree representation of Castillo et al.,
  Chemometrics and Intelligent Laboratory Systems 127 (2013) 1-6. A 2D spectrum
  is recursively split at its centre of mass into a quad-tree; two trees are
  compared node-by-node with a similarity that blends an intensity ratio and a
  shift-tolerant term.
* ``nn_peak_similarity`` -- the nearest-neighbour peak matching of Pierens et
  al., Journal of Cheminformatics 4 (2012) 25. Peaks are picked from each
  spectrum and matched to their nearest neighbour in the other, and the average
  peak-to-peak distance is turned into a similarity.

Both are self-contained apart from the Bruker reader in ``hsqc_similarity.py``.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import numpy as np

from hsqc_similarity import Spectrum2D, _axis_spacing, read_bruker_2d


# --------------------------------------------------------------------------- #
# Shared preprocessing: window, clip, area-weight, normalize.
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Windowed2D:
    ppm_f2: np.ndarray          # length nx, direct dimension
    ppm_f1: np.ndarray          # length ny, indirect dimension
    weight: np.ndarray          # shape (ny, nx), area-weighted intensity, sum == norm
    source: Path


def _window(
    spectrum: Spectrum2D,
    range_f2: tuple[float, float],
    range_f1: tuple[float, float],
    baseline: str = "clip",
    norm: float = 1.0,
) -> Windowed2D:
    intensity = spectrum.intensity.astype(np.float64, copy=True)
    if baseline == "clip":
        intensity[intensity < 0] = 0.0
    elif baseline == "abs":  # sign encodes multiplicity (edited HSQC CH2), not absence
        np.abs(intensity, out=intensity)
    elif baseline == "shift":
        intensity -= intensity.min()
    elif baseline != "none":
        raise ValueError(f"Unknown baseline mode: {baseline}")

    lo2, hi2 = sorted(range_f2)
    lo1, hi1 = sorted(range_f1)
    mask_f2 = (spectrum.ppm_f2 >= lo2) & (spectrum.ppm_f2 <= hi2)
    mask_f1 = (spectrum.ppm_f1 >= lo1) & (spectrum.ppm_f1 <= hi1)
    if not np.any(mask_f2) or not np.any(mask_f1):
        raise ValueError(f"{spectrum.source} has no points in requested ppm window")

    ppm_f2 = spectrum.ppm_f2[mask_f2]
    ppm_f1 = spectrum.ppm_f1[mask_f1]
    sub = intensity[np.ix_(mask_f1, mask_f2)]
    weight = sub * np.outer(_axis_spacing(ppm_f1), _axis_spacing(ppm_f2))
    total = float(weight.sum())
    if total <= 0:
        raise ValueError(f"{spectrum.source} has no positive integrated intensity in the window")
    weight = weight * (norm / total)
    return Windowed2D(ppm_f2=ppm_f2, ppm_f1=ppm_f1, weight=weight, source=spectrum.source)


def _overlap_ranges(
    x: Spectrum2D,
    y: Spectrum2D,
    range_f2: tuple[float, float] | None,
    range_f1: tuple[float, float] | None,
) -> tuple[tuple[float, float], tuple[float, float]]:
    def common(ax, ay, supplied, label):
        if supplied is not None:
            lo, hi = sorted(supplied)
        else:
            lo = max(float(np.min(ax)), float(np.min(ay)))
            hi = min(float(np.max(ax)), float(np.max(ay)))
        if hi <= lo:
            raise ValueError(f"The two spectra do not overlap in {label}")
        return lo, hi

    return (
        common(x.ppm_f2, y.ppm_f2, range_f2, "F2"),
        common(x.ppm_f1, y.ppm_f1, range_f1, "F1"),
    )


# --------------------------------------------------------------------------- #
# Castillo 2013 -- quad-tree similarity.
# --------------------------------------------------------------------------- #
@dataclass
class QuadNode:
    # Centre of mass in unit-square coordinates (both axes scaled to [0, 1]).
    u: float
    v: float
    intensity: float
    children: list["QuadNode | None"] = field(default_factory=lambda: [None, None, None, None])


def _build_quadtree(
    weight: np.ndarray,
    u_axis: np.ndarray,
    v_axis: np.ndarray,
    r0: int,
    r1: int,
    c0: int,
    c1: int,
    threshold: float,
    min_span: float,
    depth: int,
    max_depth: int,
) -> QuadNode | None:
    if r1 <= r0 or c1 <= c0:
        return None
    block = weight[r0:r1, c0:c1]
    total = float(block.sum())
    if total <= threshold:
        return None

    row_mass = block.sum(axis=1)
    col_mass = block.sum(axis=0)
    rows = np.arange(r0, r1)
    cols = np.arange(c0, c1)
    u = float(np.dot(u_axis[cols], col_mass) / total)
    v = float(np.dot(v_axis[rows], row_mass) / total)
    node = QuadNode(u=u, v=v, intensity=total)

    span_v = abs(float(v_axis[r1 - 1] - v_axis[r0]))
    span_u = abs(float(u_axis[c1 - 1] - u_axis[c0]))
    if depth >= max_depth or (span_v <= min_span and span_u <= min_span):
        return node

    # Split at the centre-of-mass index; force progress so the recursion terminates.
    rmid = int(round(float(np.dot(rows, row_mass) / total)))
    cmid = int(round(float(np.dot(cols, col_mass) / total)))
    rmid = min(max(rmid, r0 + 1), r1 - 1) if r1 - r0 >= 2 else r0
    cmid = min(max(cmid, c0 + 1), c1 - 1) if c1 - c0 >= 2 else c0

    quadrants = [
        (r0, rmid, c0, cmid),  # canonical order: UL, UR, LL, LR
        (r0, rmid, cmid, c1),
        (rmid, r1, c0, cmid),
        (rmid, r1, cmid, c1),
    ]
    for k, (a, b, c, d) in enumerate(quadrants):
        node.children[k] = _build_quadtree(
            weight, u_axis, v_axis, a, b, c, d, threshold, min_span, depth + 1, max_depth
        )
    return node


def build_quadtree(
    win: Windowed2D,
    threshold_frac: float = 0.005,
    min_span: float = 0.02,
    max_depth: int = 14,
) -> QuadNode | None:
    total = float(win.weight.sum())
    # Scale both axes to [0, 1] so the shift term treats F1 and F2 on equal footing.
    lo2, hi2 = float(win.ppm_f2.min()), float(win.ppm_f2.max())
    lo1, hi1 = float(win.ppm_f1.min()), float(win.ppm_f1.max())
    u_axis = (win.ppm_f2 - lo2) / (hi2 - lo2) if hi2 > lo2 else np.zeros_like(win.ppm_f2)
    v_axis = (win.ppm_f1 - lo1) / (hi1 - lo1) if hi1 > lo1 else np.zeros_like(win.ppm_f1)
    return _build_quadtree(
        win.weight, u_axis, v_axis, 0, win.weight.shape[0], 0, win.weight.shape[1],
        threshold=threshold_frac * total, min_span=min_span, depth=0, max_depth=max_depth,
    )


def _node_C(a: QuadNode | None, b: QuadNode | None, alpha: float, gamma: float) -> float:
    if a is None or b is None:
        return 0.0
    hi = max(a.intensity, b.intensity)
    ratio = min(a.intensity, b.intensity) / hi if hi > 0 else 0.0
    dist = math.hypot(a.u - b.u, a.v - b.v)
    return alpha * ratio + (1.0 - alpha) * math.exp(-gamma * dist)


def _tree_sim(
    a: QuadNode | None, b: QuadNode | None, alpha: float, beta: float, gamma: float
) -> float:
    if a is None or b is None:
        return 0.0
    c = _node_C(a, b, alpha, gamma)
    child = sum(_tree_sim(a.children[k], b.children[k], alpha, beta, gamma) for k in range(4)) / 4.0
    return beta * c + (1.0 - beta) * child


def quadtree_similarity(
    spectrum_x: Spectrum2D,
    spectrum_y: Spectrum2D,
    range_f2: tuple[float, float] | None = None,
    range_f1: tuple[float, float] | None = None,
    baseline: str = "clip",
    alpha: float = 0.1,
    beta: float = 0.33,
    gamma: float = 3.0,
    threshold_frac: float = 0.005,
) -> dict[str, object]:
    """Castillo 2013 quad-tree similarity. The raw score is not 1 for identical
    spectra, so it is normalized as s(x, y) / sqrt(s(x, x) s(y, y)); the normalized
    value is 1 for identical spectra and stays in [0, 1] for the intensity/shift
    blend used here (alpha=0.1, beta=0.33; gamma acts on unit-square distances)."""
    range_f2, range_f1 = _overlap_ranges(spectrum_x, spectrum_y, range_f2, range_f1)
    win_x = _window(spectrum_x, range_f2, range_f1, baseline=baseline)
    win_y = _window(spectrum_y, range_f2, range_f1, baseline=baseline)
    tree_x = build_quadtree(win_x, threshold_frac=threshold_frac)
    tree_y = build_quadtree(win_y, threshold_frac=threshold_frac)

    sxy = _tree_sim(tree_x, tree_y, alpha, beta, gamma)
    sxx = _tree_sim(tree_x, tree_x, alpha, beta, gamma)
    syy = _tree_sim(tree_y, tree_y, alpha, beta, gamma)
    denom = math.sqrt(sxx * syy)
    normalized = sxy / denom if denom > 0 else 0.0

    return {
        "method": "quadtree",
        "similarity": float(min(1.0, max(0.0, normalized))),
        "raw": float(sxy),
        "self_x": float(sxx),
        "self_y": float(syy),
        "n_nodes_x": _count_nodes(tree_x),
        "n_nodes_y": _count_nodes(tree_y),
        "range_f2": [float(range_f2[0]), float(range_f2[1])],
        "range_f1": [float(range_f1[0]), float(range_f1[1])],
        "source_x": str(spectrum_x.source),
        "source_y": str(spectrum_y.source),
    }


def _count_nodes(node: QuadNode | None) -> int:
    if node is None:
        return 0
    return 1 + sum(_count_nodes(child) for child in node.children)


# --------------------------------------------------------------------------- #
# Pierens 2012 -- nearest-neighbour peak matching.
# --------------------------------------------------------------------------- #
def pick_peaks(
    win: Windowed2D,
    threshold_frac: float = 0.05,
    max_peaks: int = 400,
) -> np.ndarray:
    """Return picked peaks as an array of shape (n, 3): unit-F2, unit-F1, intensity.
    A peak is a local maximum in its 3x3 neighbourhood above ``threshold_frac`` of
    the global maximum."""
    w = win.weight
    peak = float(w.max())
    if peak <= 0:
        return np.zeros((0, 3))
    interior = np.ones_like(w, dtype=bool)
    greater = np.ones_like(w, dtype=bool)
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            if dr == 0 and dc == 0:
                continue
            shifted = np.full_like(w, -np.inf)
            r_src = slice(max(0, -dr), w.shape[0] - max(0, dr))
            r_dst = slice(max(0, dr), w.shape[0] - max(0, -dr))
            c_src = slice(max(0, -dc), w.shape[1] - max(0, dc))
            c_dst = slice(max(0, dc), w.shape[1] - max(0, -dc))
            shifted[r_dst, c_dst] = w[r_src, c_src]
            greater &= w >= shifted
    mask = greater & interior & (w >= threshold_frac * peak)
    rows, cols = np.nonzero(mask)
    if rows.size == 0:
        return np.zeros((0, 3))

    lo2, hi2 = float(win.ppm_f2.min()), float(win.ppm_f2.max())
    lo1, hi1 = float(win.ppm_f1.min()), float(win.ppm_f1.max())
    u = (win.ppm_f2[cols] - lo2) / (hi2 - lo2) if hi2 > lo2 else np.zeros(cols.shape)
    v = (win.ppm_f1[rows] - lo1) / (hi1 - lo1) if hi1 > lo1 else np.zeros(rows.shape)
    inten = w[rows, cols]
    peaks = np.column_stack([u, v, inten])
    if peaks.shape[0] > max_peaks:
        keep = np.argsort(peaks[:, 2])[::-1][:max_peaks]
        peaks = peaks[keep]
    return peaks


def _mean_nn_distance(a: np.ndarray, b: np.ndarray) -> float:
    # Average nearest-neighbour distance from every peak in a to the closest peak in b.
    diff_u = a[:, 0][:, None] - b[:, 0][None, :]
    diff_v = a[:, 1][:, None] - b[:, 1][None, :]
    dist = np.hypot(diff_u, diff_v)
    return float(dist.min(axis=1).mean())


def nn_peak_similarity(
    spectrum_x: Spectrum2D,
    spectrum_y: Spectrum2D,
    range_f2: tuple[float, float] | None = None,
    range_f1: tuple[float, float] | None = None,
    baseline: str = "clip",
    threshold_frac: float = 0.05,
) -> dict[str, object]:
    """Pierens 2012 nearest-neighbour matching. Peaks are picked from both spectra,
    matched to their nearest neighbour in both directions, and the symmetric average
    peak-to-peak distance ``d`` (in unit-square coordinates) is mapped to a similarity
    ``1 / (1 + d)``, which is 1 for identical peak lists."""
    range_f2, range_f1 = _overlap_ranges(spectrum_x, spectrum_y, range_f2, range_f1)
    win_x = _window(spectrum_x, range_f2, range_f1, baseline=baseline)
    win_y = _window(spectrum_y, range_f2, range_f1, baseline=baseline)
    peaks_x = pick_peaks(win_x, threshold_frac=threshold_frac)
    peaks_y = pick_peaks(win_y, threshold_frac=threshold_frac)
    if peaks_x.shape[0] == 0 or peaks_y.shape[0] == 0:
        raise ValueError("No peaks were picked from one of the spectra; lower --threshold-frac")

    d_xy = _mean_nn_distance(peaks_x, peaks_y)
    d_yx = _mean_nn_distance(peaks_y, peaks_x)
    d = 0.5 * (d_xy + d_yx)
    return {
        "method": "nn_peak",
        "similarity": float(1.0 / (1.0 + d)),
        "mean_distance": float(d),
        "n_peaks_x": int(peaks_x.shape[0]),
        "n_peaks_y": int(peaks_y.shape[0]),
        "range_f2": [float(range_f2[0]), float(range_f2[1])],
        "range_f1": [float(range_f1[0]), float(range_f1[1])],
        "source_x": str(spectrum_x.source),
        "source_y": str(spectrum_y.source),
    }


# --------------------------------------------------------------------------- #
# CLI.
# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Alternative HSQC similarity methods (Castillo tree, Pierens nearest-neighbour)."
    )
    parser.add_argument("spectrum_x", help="Bruker experiment directory, pdata directory, or 2rr file")
    parser.add_argument("spectrum_y", help="Bruker experiment directory, pdata directory, or 2rr file")
    parser.add_argument("--method", choices=["quadtree", "nn"], default="quadtree")
    parser.add_argument("--procno", type=int, default=None)
    parser.add_argument("--f2-min", type=float, default=None)
    parser.add_argument("--f2-max", type=float, default=None)
    parser.add_argument("--f1-min", type=float, default=None)
    parser.add_argument("--f1-max", type=float, default=None)
    parser.add_argument("--baseline", choices=["clip", "abs", "shift", "none"], default="clip")
    parser.add_argument("--threshold-frac", type=float, default=None, help="Peak/leaf threshold fraction")
    parser.add_argument("--json", action="store_true")
    return parser


def _paired(lo, hi, name):
    if lo is None and hi is None:
        return None
    if lo is None or hi is None:
        raise SystemExit(f"--{name}-min and --{name}-max must be supplied together")
    return (lo, hi)


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    range_f2 = _paired(args.f2_min, args.f2_max, "f2")
    range_f1 = _paired(args.f1_min, args.f1_max, "f1")
    x = read_bruker_2d(args.spectrum_x, procno=args.procno)
    y = read_bruker_2d(args.spectrum_y, procno=args.procno)

    if args.method == "quadtree":
        kw = {} if args.threshold_frac is None else {"threshold_frac": args.threshold_frac}
        result = quadtree_similarity(x, y, range_f2=range_f2, range_f1=range_f1, baseline=args.baseline, **kw)
    else:
        kw = {} if args.threshold_frac is None else {"threshold_frac": args.threshold_frac}
        result = nn_peak_similarity(x, y, range_f2=range_f2, range_f1=range_f1, baseline=args.baseline, **kw)

    if args.json:
        print(json.dumps(result, indent=2))
    else:
        print(f"method: {result['method']}")
        print(f"similarity: {result['similarity']:.6f}")
        for key in ("raw", "mean_distance", "n_nodes_x", "n_nodes_y", "n_peaks_x", "n_peaks_y"):
            if key in result:
                print(f"{key}: {result[key]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
