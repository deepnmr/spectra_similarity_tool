#!/usr/bin/env python3
"""Calculate Bodis/Ross/Pretsch spectrum similarity for Bruker 1D spectra."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class Spectrum:
    ppm: np.ndarray
    intensity: np.ndarray
    source: Path


def parse_jcamp(path: Path) -> dict[str, object]:
    params: dict[str, object] = {}
    with path.open("r", encoding="latin-1", errors="replace") as handle:
        pending_key: str | None = None
        pending_value: list[str] = []
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("##$"):
                if pending_key is not None:
                    params[pending_key] = _parse_value(" ".join(pending_value))
                key, _, value = line[3:].partition("=")
                pending_key = key.strip()
                pending_value = [value.strip()]
            elif pending_key is not None and not line.startswith("##"):
                pending_value.append(line)
        if pending_key is not None:
            params[pending_key] = _parse_value(" ".join(pending_value))
    return params


def _parse_value(value: str) -> object:
    value = value.strip()
    if value.startswith("<") and value.endswith(">"):
        return value[1:-1]
    if value.startswith("("):
        _, _, tail = value.partition(")")
        value = tail.strip()
    parts = value.split()
    if len(parts) > 1:
        parsed = [_parse_scalar(part) for part in parts]
        return parsed
    return _parse_scalar(value)


def _parse_scalar(value: str) -> object:
    try:
        if any(ch in value for ch in ".eE"):
            return float(value)
        return int(value)
    except ValueError:
        return value


def _number(params: dict[str, object], key: str, default: float | None = None) -> float:
    value = params.get(key, default)
    if value is None:
        raise ValueError(f"Missing Bruker parameter {key}")
    if isinstance(value, list):
        value = value[0]
    return float(value)


def find_processed_1r(path: Path, procno: int | None = None) -> Path:
    if path.is_file() and path.name == "1r":
        return path
    candidates: list[Path] = []
    if procno is not None:
        candidates.extend([path / "pdata" / str(procno) / "1r", path / str(procno) / "1r"])
    candidates.extend(sorted(path.glob("pdata/*/1r")))
    candidates.extend(sorted(path.glob("*/1r")))
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(f"Could not find processed Bruker 1D file '1r' under {path}")


def read_bruker_1d(path: str | Path, procno: int | None = None) -> Spectrum:
    data_path = find_processed_1r(Path(path).expanduser().resolve(), procno=procno)
    proc_dir = data_path.parent
    procs_path = proc_dir / "procs"
    if not procs_path.is_file():
        raise FileNotFoundError(f"Missing Bruker procs file next to {data_path}")

    procs = parse_jcamp(procs_path)
    size = int(_number(procs, "SI"))
    offset = _number(procs, "OFFSET")
    sw_hz = _number(procs, "SW_p")
    sf_mhz = _number(procs, "SF")
    endian = ">" if int(_number(procs, "BYTORDP", 0)) else "<"
    dtyp = int(_number(procs, "DTYPP", 0))
    dtype = np.dtype(f"{endian}{'f8' if dtyp == 2 else 'i4'}")

    intensity = np.fromfile(data_path, dtype=dtype, count=size).astype(np.float64)
    if intensity.size != size:
        raise ValueError(f"{data_path} contains {intensity.size} points, expected {size}")

    nc_proc = int(_number(procs, "NC_proc", 0))
    if nc_proc:
        intensity *= math.pow(2.0, nc_proc)

    sw_ppm = sw_hz / sf_mhz
    # Bruker processed 1D spectra are stored left-to-right from high to low ppm.
    ppm = offset - np.arange(size, dtype=np.float64) * (sw_ppm / size)
    return Spectrum(ppm=ppm, intensity=intensity, source=data_path)


def prepare_spectrum(
    spectrum: Spectrum,
    ppm_min: float,
    ppm_max: float,
    baseline: str = "clip",
) -> tuple[np.ndarray, np.ndarray]:
    ppm = spectrum.ppm
    intensity = spectrum.intensity.copy()
    if baseline == "clip":
        intensity[intensity < 0] = 0.0
    elif baseline == "shift":
        intensity -= intensity.min()
    elif baseline != "none":
        raise ValueError(f"Unknown baseline mode: {baseline}")

    lo, hi = sorted((ppm_min, ppm_max))
    mask = (ppm >= lo) & (ppm <= hi)
    if not np.any(mask):
        raise ValueError(f"{spectrum.source} has no points in requested ppm range {lo:g}..{hi:g}")
    ppm = ppm[mask]
    intensity = intensity[mask]
    order = np.argsort(ppm)
    return ppm[order], intensity[order]


def bin_integrals(ppm: np.ndarray, intensity: np.ndarray, edges: np.ndarray) -> np.ndarray:
    indices = np.searchsorted(edges, ppm, side="right") - 1
    indices[ppm == edges[-1]] = len(edges) - 2
    valid = (indices >= 0) & (indices < len(edges) - 1)
    if not np.any(valid):
        return np.zeros(len(edges) - 1, dtype=np.float64)
    weights = point_weights(ppm, intensity)
    return np.bincount(indices[valid], weights=weights[valid], minlength=len(edges) - 1)


def point_weights(ppm: np.ndarray, intensity: np.ndarray) -> np.ndarray:
    # Weight by local point spacing so spectra with different digital resolution compare fairly.
    spacing = np.empty_like(ppm)
    spacing[1:-1] = (ppm[2:] - ppm[:-2]) / 2.0
    spacing[0] = ppm[1] - ppm[0] if ppm.size > 1 else 1.0
    spacing[-1] = ppm[-1] - ppm[-2] if ppm.size > 1 else 1.0
    return intensity * np.abs(spacing)


def upper_envelope(values: np.ndarray) -> np.ndarray:
    """Piecewise-linear upper envelope through remaining maxima and the point (1, 1)."""
    x = np.arange(1, len(values) + 1, dtype=np.float64)
    y = values.astype(np.float64).copy()
    y[0] = 1.0

    points: list[int] = [0]
    pos = 1
    while pos < len(y):
        suffix = y[pos:]
        max_value = np.max(suffix)
        # Choose the farthest equal maximum to avoid artificial short oscillations.
        relative = np.flatnonzero(np.isclose(suffix, max_value, rtol=0.0, atol=1e-12))[-1]
        idx = pos + int(relative)
        points.append(idx)
        pos = idx + 1
    if points[-1] != len(y) - 1:
        points.append(len(y) - 1)

    envelope = y.copy()
    for left, right in zip(points, points[1:]):
        if right == left:
            continue
        envelope[left : right + 1] = np.interp(x[left : right + 1], [x[left], x[right]], [y[left], y[right]])
    return np.maximum(envelope, y)


def spectrum_similarity(
    spectrum_x: Spectrum,
    spectrum_y: Spectrum,
    min_bin_width: float = 0.4,
    ppm_range: tuple[float, float] | None = None,
    norm_x: float = 1.0,
    norm_y: float = 1.0,
    baseline: str = "clip",
) -> dict[str, object]:
    if ppm_range is None:
        lo = max(float(np.min(spectrum_x.ppm)), float(np.min(spectrum_y.ppm)))
        hi = min(float(np.max(spectrum_x.ppm)), float(np.max(spectrum_y.ppm)))
    else:
        lo, hi = sorted(ppm_range)
    if hi <= lo:
        raise ValueError("The two spectra do not have an overlapping ppm range")

    ppm_x, int_x = prepare_spectrum(spectrum_x, lo, hi, baseline=baseline)
    ppm_y, int_y = prepare_spectrum(spectrum_y, lo, hi, baseline=baseline)
    total_x = float(np.sum(point_weights(ppm_x, int_x)))
    total_y = float(np.sum(point_weights(ppm_y, int_y)))
    if total_x <= 0 or total_y <= 0:
        raise ValueError("Both spectra must have positive integrated intensity after preprocessing")
    int_x = int_x * (norm_x / total_x)
    int_y = int_y * (norm_y / total_y)

    width = hi - lo
    max_bins = max(1, int(math.floor(width / min_bin_width)))
    si = np.empty(max_bins, dtype=np.float64)
    for n_bins in range(1, max_bins + 1):
        edges = np.linspace(lo, hi, n_bins + 1)
        bx = bin_integrals(ppm_x, int_x, edges)
        by = bin_integrals(ppm_y, int_y, edges)
        overlap = float(np.minimum(bx, by).sum())
        denominator = norm_x + norm_y - overlap
        value = overlap / denominator if denominator > 0 else 0.0
        si[n_bins - 1] = min(1.0, max(0.0, value))
    si_star = np.clip(upper_envelope(si), 0.0, 1.0)
    return {
        "similarity": float(np.mean(si_star)),
        "n_bins": int(max_bins),
        "ppm_range": [float(lo), float(hi)],
        "min_bin_width": float(min_bin_width),
        "si": si.tolist(),
        "si_star": si_star.tolist(),
        "source_x": str(spectrum_x.source),
        "source_y": str(spectrum_y.source),
    }


def normalized_window(
    spectrum: Spectrum,
    ppm_range: tuple[float, float],
    target_integral: float = 1.0,
    baseline: str = "clip",
) -> tuple[np.ndarray, np.ndarray]:
    lo, hi = sorted(ppm_range)
    ppm, intensity = prepare_spectrum(spectrum, lo, hi, baseline=baseline)
    total = float(np.sum(point_weights(ppm, intensity)))
    if total <= 0:
        raise ValueError(f"{spectrum.source} has no positive integrated intensity in the plot range")
    return ppm, intensity * (target_integral / total)


def save_similarity_plot(
    result: dict[str, object],
    spectrum_x: Spectrum,
    spectrum_y: Spectrum,
    output: str | Path | None = None,
    norm_x: float = 1.0,
    norm_y: float = 1.0,
    baseline: str = "clip",
    show: bool = False,
) -> Path | None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("Plotting requires matplotlib. Install it with: python3 -m pip install matplotlib") from exc

    ppm_range_values = result["ppm_range"]
    ppm_range = (float(ppm_range_values[0]), float(ppm_range_values[1]))
    ppm_x, int_x = normalized_window(spectrum_x, ppm_range, target_integral=norm_x, baseline=baseline)
    ppm_y, int_y = normalized_window(spectrum_y, ppm_range, target_integral=norm_y, baseline=baseline)

    si = np.asarray(result["si"], dtype=np.float64)
    si_star = np.asarray(result["si_star"], dtype=np.float64)
    bins = np.arange(1, si.size + 1)

    fig, axes = plt.subplots(2, 1, figsize=(10, 7), constrained_layout=True)
    axes[0].plot(ppm_x, int_x, label="spectrum x", linewidth=1.0)
    axes[0].plot(ppm_y, int_y, label="spectrum y", linewidth=1.0, alpha=0.8)
    axes[0].invert_xaxis()
    axes[0].set_xlabel("ppm")
    axes[0].set_ylabel("normalized intensity")
    axes[0].legend(loc="upper right")
    axes[0].set_title(f"Spectrum comparison, similarity = {float(result['similarity']):.4f}")

    axes[1].plot(bins, si, label="SI_n", linewidth=1.0, alpha=0.7)
    axes[1].plot(bins, si_star, label="SI*_n envelope", linewidth=2.0)
    axes[1].set_xlabel("number of bins")
    axes[1].set_ylabel("similarity index")
    axes[1].set_ylim(-0.02, 1.02)
    axes[1].legend(loc="best")
    axes[1].grid(True, linewidth=0.4, alpha=0.4)

    saved_path = None
    if output is not None:
        saved_path = Path(output).expanduser().resolve()
        fig.savefig(saved_path, dpi=200)
    if show:
        plt.show()
    plt.close(fig)
    return saved_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Calculate the Bodis/Ross/Pretsch bin-method similarity of two Bruker 1D spectra."
    )
    parser.add_argument("spectrum_x", help="Bruker experiment directory, pdata directory, or 1r file")
    parser.add_argument("spectrum_y", help="Bruker experiment directory, pdata directory, or 1r file")
    parser.add_argument("--procno", type=int, default=None, help="Processed spectrum number under pdata")
    parser.add_argument("--min-bin-width", type=float, default=0.4, help="Smallest bin width in ppm")
    parser.add_argument("--ppm-min", type=float, default=None, help="Lower ppm limit")
    parser.add_argument("--ppm-max", type=float, default=None, help="Upper ppm limit")
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
        const="spectrum_similarity.png",
        default=None,
        help="Save a PNG graph. Optional path defaults to spectrum_similarity.png",
    )
    parser.add_argument("--show", action="store_true", help="Display the graph interactively")
    parser.add_argument("--json", action="store_true", help="Print full result as JSON")
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    ppm_range = None
    if args.ppm_min is not None or args.ppm_max is not None:
        if args.ppm_min is None or args.ppm_max is None:
            raise SystemExit("--ppm-min and --ppm-max must be supplied together")
        ppm_range = (args.ppm_min, args.ppm_max)

    spectrum_x = read_bruker_1d(args.spectrum_x, procno=args.procno)
    spectrum_y = read_bruker_1d(args.spectrum_y, procno=args.procno)
    result = spectrum_similarity(
        spectrum_x,
        spectrum_y,
        min_bin_width=args.min_bin_width,
        ppm_range=ppm_range,
        norm_x=args.norm_x,
        norm_y=args.norm_y,
        baseline=args.baseline,
    )
    if args.json:
        print(json.dumps(result, indent=2))
    else:
        print(f"similarity: {result['similarity']:.6f}")
        print(f"bins: {result['n_bins']}  ppm_range: {result['ppm_range'][0]:.6g}..{result['ppm_range'][1]:.6g}")
        print(f"x: {result['source_x']}")
        print(f"y: {result['source_y']}")
    if args.plot is not None or args.show:
        saved_path = save_similarity_plot(
            result,
            spectrum_x,
            spectrum_y,
            output=args.plot,
            norm_x=args.norm_x,
            norm_y=args.norm_y,
            baseline=args.baseline,
            show=args.show,
        )
        if saved_path is not None:
            print(f"plot: {saved_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
