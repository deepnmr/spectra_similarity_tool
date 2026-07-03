#!/usr/bin/env python3
"""Calculate Bodis/Ross/Pretsch spectrum similarity for Bruker 2D HSQC spectra.

This extends the 1D bin-method similarity to two dimensions. Each spectrum is
subdivided into an ``n`` by ``n`` grid of bins whose widths shrink as ``n``
grows, exactly as the 1D method subdivides a chemical-shift axis into ``n``
bins. The per-grid similarity index ``SI_n`` and its upper envelope ``SI*_n``
are averaged into a single similarity score in ``[0, 1]``.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from spectrum_similarity import _number, parse_jcamp, upper_envelope


@dataclass(frozen=True)
class Spectrum2D:
    ppm_f2: np.ndarray  # direct dimension axis (typically 1H), length nx
    ppm_f1: np.ndarray  # indirect dimension axis (typically 13C/15N), length ny
    intensity: np.ndarray  # shape (ny, nx), intensity[f1_index, f2_index]
    source: Path


def find_processed_2rr(path: Path, procno: int | None = None) -> Path:
    if path.is_file() and path.name == "2rr":
        return path
    candidates: list[Path] = []
    if procno is not None:
        candidates.extend([path / "pdata" / str(procno) / "2rr", path / str(procno) / "2rr"])
    candidates.extend(sorted(path.glob("pdata/*/2rr")))
    candidates.extend(sorted(path.glob("*/2rr")))
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(f"Could not find processed Bruker 2D file '2rr' under {path}")


def _axis_ppm(procs: dict[str, object], size: int) -> np.ndarray:
    offset = _number(procs, "OFFSET")
    sw_hz = _number(procs, "SW_p")
    sf_mhz = _number(procs, "SF")
    sw_ppm = sw_hz / sf_mhz
    return offset - np.arange(size, dtype=np.float64) * (sw_ppm / size)


def _deinterleave_submatrix(
    raw: np.ndarray, si_f1: int, si_f2: int, xdim_f1: int, xdim_f2: int
) -> np.ndarray:
    """Reassemble a Bruker processed 2D matrix from its submatrix (tile) layout."""
    if xdim_f1 <= 0 or xdim_f1 > si_f1:
        xdim_f1 = si_f1
    if xdim_f2 <= 0 or xdim_f2 > si_f2:
        xdim_f2 = si_f2
    tiles_f1 = si_f1 // xdim_f1
    tiles_f2 = si_f2 // xdim_f2
    if tiles_f1 * xdim_f1 != si_f1 or tiles_f2 * xdim_f2 != si_f2:
        # Fall back to a plain row-major matrix when tiling does not divide evenly.
        return raw.reshape(si_f1, si_f2)

    matrix = np.empty((si_f1, si_f2), dtype=raw.dtype)
    tile_size = xdim_f1 * xdim_f2
    index = 0
    for tf1 in range(tiles_f1):
        row0 = tf1 * xdim_f1
        for tf2 in range(tiles_f2):
            col0 = tf2 * xdim_f2
            tile = raw[index : index + tile_size].reshape(xdim_f1, xdim_f2)
            matrix[row0 : row0 + xdim_f1, col0 : col0 + xdim_f2] = tile
            index += tile_size
    return matrix


def read_bruker_2d(path: str | Path, procno: int | None = None) -> Spectrum2D:
    data_path = find_processed_2rr(Path(path).expanduser().resolve(), procno=procno)
    proc_dir = data_path.parent
    procs_path = proc_dir / "procs"      # direct dimension F2
    proc2s_path = proc_dir / "proc2s"    # indirect dimension F1
    if not procs_path.is_file():
        raise FileNotFoundError(f"Missing Bruker procs file next to {data_path}")
    if not proc2s_path.is_file():
        raise FileNotFoundError(f"Missing Bruker proc2s file next to {data_path}")

    procs = parse_jcamp(procs_path)
    proc2s = parse_jcamp(proc2s_path)

    si_f2 = int(_number(procs, "SI"))
    si_f1 = int(_number(proc2s, "SI"))
    xdim_f2 = int(_number(procs, "XDIM", si_f2))
    xdim_f1 = int(_number(proc2s, "XDIM", si_f1))

    endian = ">" if int(_number(procs, "BYTORDP", 0)) else "<"
    dtyp = int(_number(procs, "DTYPP", 0))
    dtype = np.dtype(f"{endian}{'f8' if dtyp == 2 else 'i4'}")

    count = si_f1 * si_f2
    raw = np.fromfile(data_path, dtype=dtype, count=count).astype(np.float64)
    if raw.size != count:
        raise ValueError(f"{data_path} contains {raw.size} points, expected {count}")

    matrix = _deinterleave_submatrix(raw, si_f1, si_f2, xdim_f1, xdim_f2)

    nc_proc = int(_number(procs, "NC_proc", 0))
    if nc_proc:
        matrix = matrix * math.pow(2.0, nc_proc)

    ppm_f2 = _axis_ppm(procs, si_f2)
    ppm_f1 = _axis_ppm(proc2s, si_f1)
    return Spectrum2D(ppm_f2=ppm_f2, ppm_f1=ppm_f1, intensity=matrix, source=data_path)


def _axis_spacing(ppm: np.ndarray) -> np.ndarray:
    spacing = np.empty_like(ppm)
    if ppm.size > 1:
        spacing[1:-1] = (ppm[2:] - ppm[:-2]) / 2.0
        spacing[0] = ppm[1] - ppm[0]
        spacing[-1] = ppm[-1] - ppm[-2]
    else:
        spacing[:] = 1.0
    return np.abs(spacing)


@dataclass(frozen=True)
class PreparedSpectrum2D:
    """Flattened points inside a ppm window, weighted by local grid area."""

    f2: np.ndarray  # direct-dim ppm of each point
    f1: np.ndarray  # indirect-dim ppm of each point
    weight: np.ndarray  # intensity * area of each point


def prepare_spectrum_2d(
    spectrum: Spectrum2D,
    range_f2: tuple[float, float],
    range_f1: tuple[float, float],
    baseline: str = "clip",
) -> PreparedSpectrum2D:
    intensity = spectrum.intensity.astype(np.float64, copy=True)
    if baseline == "clip":
        intensity[intensity < 0] = 0.0
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

    area = np.outer(_axis_spacing(ppm_f1), _axis_spacing(ppm_f2))
    weight = sub * area

    grid_f2, grid_f1 = np.meshgrid(ppm_f2, ppm_f1)
    return PreparedSpectrum2D(
        f2=grid_f2.ravel(),
        f1=grid_f1.ravel(),
        weight=weight.ravel(),
    )


def bin_integrals_2d(
    prepared: PreparedSpectrum2D,
    edges_f2: np.ndarray,
    edges_f1: np.ndarray,
) -> np.ndarray:
    hist, _, _ = np.histogram2d(
        prepared.f1, prepared.f2, bins=(edges_f1, edges_f2), weights=prepared.weight
    )
    return hist


def hsqc_similarity(
    spectrum_x: Spectrum2D,
    spectrum_y: Spectrum2D,
    min_bin_width_f2: float = 0.4,
    min_bin_width_f1: float = 4.0,
    range_f2: tuple[float, float] | None = None,
    range_f1: tuple[float, float] | None = None,
    norm_x: float = 1.0,
    norm_y: float = 1.0,
    baseline: str = "clip",
) -> dict[str, object]:
    range_f2 = _overlap_range(spectrum_x.ppm_f2, spectrum_y.ppm_f2, range_f2, "F2")
    range_f1 = _overlap_range(spectrum_x.ppm_f1, spectrum_y.ppm_f1, range_f1, "F1")

    prep_x = prepare_spectrum_2d(spectrum_x, range_f2, range_f1, baseline=baseline)
    prep_y = prepare_spectrum_2d(spectrum_y, range_f2, range_f1, baseline=baseline)
    total_x = float(prep_x.weight.sum())
    total_y = float(prep_y.weight.sum())
    if total_x <= 0 or total_y <= 0:
        raise ValueError("Both spectra must have positive integrated intensity after preprocessing")
    prep_x = PreparedSpectrum2D(prep_x.f2, prep_x.f1, prep_x.weight * (norm_x / total_x))
    prep_y = PreparedSpectrum2D(prep_y.f2, prep_y.f1, prep_y.weight * (norm_y / total_y))

    lo2, hi2 = sorted(range_f2)
    lo1, hi1 = sorted(range_f1)
    width_f2 = hi2 - lo2
    width_f1 = hi1 - lo1
    max_n = min(
        max(1, int(math.floor(width_f2 / min_bin_width_f2))),
        max(1, int(math.floor(width_f1 / min_bin_width_f1))),
    )

    si = np.empty(max_n, dtype=np.float64)
    for n in range(1, max_n + 1):
        edges_f2 = np.linspace(lo2, hi2, n + 1)
        edges_f1 = np.linspace(lo1, hi1, n + 1)
        bx = bin_integrals_2d(prep_x, edges_f2, edges_f1)
        by = bin_integrals_2d(prep_y, edges_f2, edges_f1)
        overlap = float(np.minimum(bx, by).sum())
        denominator = norm_x + norm_y - overlap
        value = overlap / denominator if denominator > 0 else 0.0
        si[n - 1] = min(1.0, max(0.0, value))
    si_star = np.clip(upper_envelope(si), 0.0, 1.0)

    return {
        "similarity": float(np.mean(si_star)),
        "n_grids": int(max_n),
        "range_f2": [float(lo2), float(hi2)],
        "range_f1": [float(lo1), float(hi1)],
        "min_bin_width_f2": float(min_bin_width_f2),
        "min_bin_width_f1": float(min_bin_width_f1),
        "si": si.tolist(),
        "si_star": si_star.tolist(),
        "source_x": str(spectrum_x.source),
        "source_y": str(spectrum_y.source),
    }


def _overlap_range(
    axis_x: np.ndarray,
    axis_y: np.ndarray,
    supplied: tuple[float, float] | None,
    label: str,
) -> tuple[float, float]:
    if supplied is not None:
        lo, hi = sorted(supplied)
    else:
        lo = max(float(np.min(axis_x)), float(np.min(axis_y)))
        hi = min(float(np.max(axis_x)), float(np.max(axis_y)))
    if hi <= lo:
        raise ValueError(f"The two spectra do not have an overlapping ppm range in {label}")
    return lo, hi


def save_similarity_plot_2d(
    result: dict[str, object],
    spectrum_x: Spectrum2D,
    spectrum_y: Spectrum2D,
    output: str | Path | None = None,
    baseline: str = "clip",
    contour_levels: int = 12,
    show: bool = False,
) -> Path | None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError(
            "Plotting requires matplotlib. Install it with: python3 -m pip install matplotlib"
        ) from exc

    range_f2 = (float(result["range_f2"][0]), float(result["range_f2"][1]))
    range_f1 = (float(result["range_f1"][0]), float(result["range_f1"][1]))
    si = np.asarray(result["si"], dtype=np.float64)
    si_star = np.asarray(result["si_star"], dtype=np.float64)
    grids = np.arange(1, si.size + 1)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
    for ax, spectrum, title in (
        (axes[0], spectrum_x, "spectrum x"),
        (axes[1], spectrum_y, "spectrum y"),
    ):
        _contour_window(ax, spectrum, range_f2, range_f1, baseline, contour_levels)
        ax.set_title(title)
        ax.set_xlabel("F2 ppm")
        ax.set_ylabel("F1 ppm")
        ax.invert_xaxis()
        ax.invert_yaxis()

    axes[2].plot(grids, si, label="SI_n", linewidth=1.0, alpha=0.7)
    axes[2].plot(grids, si_star, label="SI*_n envelope", linewidth=2.0)
    axes[2].set_xlabel("bins per dimension")
    axes[2].set_ylabel("similarity index")
    axes[2].set_ylim(-0.02, 1.02)
    axes[2].legend(loc="best")
    axes[2].grid(True, linewidth=0.4, alpha=0.4)
    axes[2].set_title(f"similarity = {float(result['similarity']):.4f}")

    saved_path = None
    if output is not None:
        saved_path = Path(output).expanduser().resolve()
        fig.savefig(saved_path, dpi=200)
    if show:
        plt.show()
    plt.close(fig)
    return saved_path


def _contour_window(ax, spectrum, range_f2, range_f1, baseline, contour_levels) -> None:
    intensity = spectrum.intensity.astype(np.float64, copy=True)
    if baseline == "clip":
        intensity[intensity < 0] = 0.0
    elif baseline == "shift":
        intensity -= intensity.min()

    lo2, hi2 = sorted(range_f2)
    lo1, hi1 = sorted(range_f1)
    mask_f2 = (spectrum.ppm_f2 >= lo2) & (spectrum.ppm_f2 <= hi2)
    mask_f1 = (spectrum.ppm_f1 >= lo1) & (spectrum.ppm_f1 <= hi1)
    ppm_f2 = spectrum.ppm_f2[mask_f2]
    ppm_f1 = spectrum.ppm_f1[mask_f1]
    sub = intensity[np.ix_(mask_f1, mask_f2)]
    peak = float(sub.max())
    if peak <= 0:
        return
    levels = np.linspace(peak * 0.05, peak, contour_levels)
    ax.contour(ppm_f2, ppm_f1, sub, levels=levels, linewidths=0.6)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Calculate the Bodis/Ross/Pretsch bin-method similarity of two Bruker 2D HSQC spectra."
    )
    parser.add_argument("spectrum_x", help="Bruker experiment directory, pdata directory, or 2rr file")
    parser.add_argument("spectrum_y", help="Bruker experiment directory, pdata directory, or 2rr file")
    parser.add_argument("--procno", type=int, default=None, help="Processed spectrum number under pdata")
    parser.add_argument("--min-bin-width-f2", type=float, default=0.4, help="Smallest F2 bin width in ppm")
    parser.add_argument("--min-bin-width-f1", type=float, default=4.0, help="Smallest F1 bin width in ppm")
    parser.add_argument("--f2-min", type=float, default=None, help="Lower F2 ppm limit")
    parser.add_argument("--f2-max", type=float, default=None, help="Upper F2 ppm limit")
    parser.add_argument("--f1-min", type=float, default=None, help="Lower F1 ppm limit")
    parser.add_argument("--f1-max", type=float, default=None, help="Upper F1 ppm limit")
    parser.add_argument("--norm-x", type=float, default=1.0, help="Target integral for first spectrum")
    parser.add_argument("--norm-y", type=float, default=1.0, help="Target integral for second spectrum")
    parser.add_argument(
        "--baseline",
        choices=["clip", "shift", "none"],
        default="clip",
        help="Negative intensity handling before integration",
    )
    parser.add_argument(
        "--plot",
        nargs="?",
        const="hsqc_similarity.png",
        default=None,
        help="Save a PNG graph. Optional path defaults to hsqc_similarity.png",
    )
    parser.add_argument("--show", action="store_true", help="Display the graph interactively")
    parser.add_argument("--json", action="store_true", help="Print full result as JSON")
    return parser


def _paired_range(
    lo: float | None, hi: float | None, name: str
) -> tuple[float, float] | None:
    if lo is None and hi is None:
        return None
    if lo is None or hi is None:
        raise SystemExit(f"--{name}-min and --{name}-max must be supplied together")
    return (lo, hi)


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    range_f2 = _paired_range(args.f2_min, args.f2_max, "f2")
    range_f1 = _paired_range(args.f1_min, args.f1_max, "f1")

    spectrum_x = read_bruker_2d(args.spectrum_x, procno=args.procno)
    spectrum_y = read_bruker_2d(args.spectrum_y, procno=args.procno)
    result = hsqc_similarity(
        spectrum_x,
        spectrum_y,
        min_bin_width_f2=args.min_bin_width_f2,
        min_bin_width_f1=args.min_bin_width_f1,
        range_f2=range_f2,
        range_f1=range_f1,
        norm_x=args.norm_x,
        norm_y=args.norm_y,
        baseline=args.baseline,
    )
    if args.json:
        print(json.dumps(result, indent=2))
    else:
        print(f"similarity: {result['similarity']:.6f}")
        print(f"grids: {result['n_grids']}")
        print(
            f"F2 range: {result['range_f2'][0]:.6g}..{result['range_f2'][1]:.6g}  "
            f"F1 range: {result['range_f1'][0]:.6g}..{result['range_f1'][1]:.6g}"
        )
        print(f"x: {result['source_x']}")
        print(f"y: {result['source_y']}")
    if args.plot is not None or args.show:
        saved_path = save_similarity_plot_2d(
            result,
            spectrum_x,
            spectrum_y,
            output=args.plot,
            baseline=args.baseline,
            show=args.show,
        )
        if saved_path is not None:
            print(f"plot: {saved_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
