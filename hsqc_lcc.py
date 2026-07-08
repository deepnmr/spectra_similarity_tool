#!/usr/bin/env python3
"""Lineshape Correlation Coefficient (LCC) similarity for Bruker 2D HSQC spectra.

A synthesis of the three literature methods (Bodis bins, Castillo tree, Pierens NN),
designed to keep the bin method's discrimination on dense protein fingerprints while
removing its shift-brittleness and without the tree/NN saturation.

Idea: render each spectrum to ONE shared grid, blur by the physical NMR linewidth
(a Gaussian per axis), then score with the mean-centred normalized cross-correlation
(Pearson / ZNCC) of the two images at zero lag -- no shift search.

Why it works where bins/tree/NN fail:
* No hard bin edges -- the Gaussian render is continuous, so a small titration drift
  lowers the correlation smoothly instead of dropping mass across a bin boundary.
* Mean-centring rewards *co-located* intensity and penalises intensity where the
  other spectrum has none, so a differently-scattered protein genuinely decorrelates
  (this is the ~+0.17 separation step over an un-centred overlap).
* Zero-lag (no alignment): aligning would let a different protein slide into
  registration and saturate, exactly the tree/NN failure. We deliberately do not align.

Self-similarity is exactly 1 (Pearson of an image with itself).

References: Bodis 2009 (bins / area weighting), Castillo 2013 (shift tolerance via a
physical blur width), Pierens 2012 (peak/lineshape view). The blur width is
sigma = sqrt(linewidth^2 + expected_drift^2) per axis (Williamson amide-CSP tolerance).
"""
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from hsqc_similarity import _gaussian_kernel, _smooth_axis, read_bruker_2d, Spectrum2D
from hsqc_methods import _overlap_ranges, _window


@dataclass(frozen=True)
class RenderedImage:
    image: np.ndarray  # blurred, shape (n1, n2) with F1 rows, F2 cols
    step_f2: float     # ppm per pixel along F2
    step_f1: float     # ppm per pixel along F1


def render_image(
    spectrum: Spectrum2D,
    range_f2: tuple[float, float],
    range_f1: tuple[float, float],
    sigma_f2: float = 0.03,
    sigma_f1: float = 0.30,
    step_f2: float = 0.01,
    step_f1: float = 0.10,
    baseline: str = "clip",
    intensity_p: float = 1.0,
) -> RenderedImage:
    """Rebin a spectrum onto a shared fixed grid and blur by the physical linewidth."""
    win = _window(spectrum, range_f2, range_f1, baseline=baseline, norm=1.0)

    lo2, hi2 = sorted(range_f2)
    lo1, hi1 = sorted(range_f1)
    nb2 = max(1, int(round((hi2 - lo2) / step_f2)))
    nb1 = max(1, int(round((hi1 - lo1) / step_f1)))
    edges2 = np.linspace(lo2, hi2, nb2 + 1)
    edges1 = np.linspace(lo1, hi1, nb1 + 1)
    dx2 = (hi2 - lo2) / nb2
    dx1 = (hi1 - lo1) / nb1

    # Render both spectra onto the SAME grid so images are directly comparable
    # regardless of native digital resolution.
    grid_f2, grid_f1 = np.meshgrid(win.ppm_f2, win.ppm_f1)
    image, _, _ = np.histogram2d(
        grid_f1.ravel(), grid_f2.ravel(), bins=(edges1, edges2), weights=win.weight.ravel()
    )

    if intensity_p != 1.0:
        image = np.power(image, intensity_p)

    # Blur = the lineshape. Separable Gaussian, sigma given in ppm -> pixels.
    kernel_f2 = _gaussian_kernel(sigma_f2 / dx2)
    if kernel_f2 is not None:
        image = _smooth_axis(image, kernel_f2, axis=1)
    kernel_f1 = _gaussian_kernel(sigma_f1 / dx1)
    if kernel_f1 is not None:
        image = _smooth_axis(image, kernel_f1, axis=0)

    return RenderedImage(image=image, step_f2=dx2, step_f1=dx1)


def _zncc(a: np.ndarray, b: np.ndarray) -> float:
    """Mean-centred normalized cross-correlation (Pearson) at zero lag, clamped to [0,1].

    Self-value is exactly 1: with a == b the ratio is ||ah||^2 / ||ah||^2, and the
    min(1, .) clamp absorbs the last-ULP fp wobble so identical images return 1.0."""
    ah = a - a.mean()
    bh = b - b.mean()
    denom = float(np.sqrt(np.sum(ah * ah) * np.sum(bh * bh)))
    if denom <= 0:  # a flat image has no lineshape to correlate
        return 0.0
    return min(1.0, max(0.0, float(np.sum(ah * bh) / denom)))


def lcc_similarity(
    spectrum_x: Spectrum2D,
    spectrum_y: Spectrum2D,
    range_f2: tuple[float, float] | None = None,
    range_f1: tuple[float, float] | None = None,
    sigma_f2: float = 0.03,
    sigma_f1: float = 0.30,
    step_f2: float = 0.01,
    step_f1: float = 0.10,
    baseline: str = "clip",
    intensity_p: float = 1.0,
) -> dict[str, object]:
    """Lineshape Correlation Coefficient similarity in [0, 1]; 1 for identical spectra."""
    range_f2, range_f1 = _overlap_ranges(spectrum_x, spectrum_y, range_f2, range_f1)
    render = lambda s: render_image(
        s, range_f2, range_f1, sigma_f2=sigma_f2, sigma_f1=sigma_f1,
        step_f2=step_f2, step_f1=step_f1, baseline=baseline, intensity_p=intensity_p,
    )
    img_x = render(spectrum_x)
    img_y = render(spectrum_y)
    similarity = _zncc(img_x.image, img_y.image)

    return {
        "method": "lcc",
        "similarity": float(similarity),
        "sigma_f2": float(sigma_f2),
        "sigma_f1": float(sigma_f1),
        "step_f2": float(img_x.step_f2),
        "step_f1": float(img_x.step_f1),
        "intensity_p": float(intensity_p),
        "grid": [int(img_x.image.shape[0]), int(img_x.image.shape[1])],
        "range_f2": [float(range_f2[0]), float(range_f2[1])],
        "range_f1": [float(range_f1[0]), float(range_f1[1])],
        "source_x": str(spectrum_x.source),
        "source_y": str(spectrum_y.source),
    }


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Lineshape Correlation Coefficient (LCC) similarity of two Bruker 2D HSQC spectra."
    )
    p.add_argument("spectrum_x", help="Bruker experiment directory, pdata directory, or 2rr file")
    p.add_argument("spectrum_y", help="Bruker experiment directory, pdata directory, or 2rr file")
    p.add_argument("--procno", type=int, default=None)
    p.add_argument("--f2-min", type=float, default=None)
    p.add_argument("--f2-max", type=float, default=None)
    p.add_argument("--f1-min", type=float, default=None)
    p.add_argument("--f1-max", type=float, default=None)
    p.add_argument("--sigma-f2", type=float, default=0.03, help="1H lineshape+drift blur sigma in ppm")
    p.add_argument("--sigma-f1", type=float, default=0.30, help="15N/13C lineshape+drift blur sigma in ppm")
    p.add_argument("--step-f2", type=float, default=0.01, help="F2 render pixel width in ppm")
    p.add_argument("--step-f1", type=float, default=0.10, help="F1 render pixel width in ppm")
    p.add_argument("--p", type=float, default=1.0, help="Intensity compression exponent (1.0 = off)")
    p.add_argument("--baseline", choices=["clip", "shift", "none"], default="clip")
    p.add_argument("--json", action="store_true")
    return p


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
    result = lcc_similarity(
        x, y, range_f2=range_f2, range_f1=range_f1,
        sigma_f2=args.sigma_f2, sigma_f1=args.sigma_f1,
        step_f2=args.step_f2, step_f1=args.step_f1,
        baseline=args.baseline, intensity_p=args.p,
    )
    if args.json:
        print(json.dumps(result, indent=2))
    else:
        print(f"method: lcc")
        print(f"similarity: {result['similarity']:.6f}")
        print(f"grid: {result['grid'][0]}x{result['grid'][1]}  "
              f"sigma: {result['sigma_f2']}/{result['sigma_f1']} ppm")
        print(f"F2 range: {result['range_f2'][0]:.6g}..{result['range_f2'][1]:.6g}  "
              f"F1 range: {result['range_f1'][0]:.6g}..{result['range_f1'][1]:.6g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
