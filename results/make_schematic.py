#!/usr/bin/env python3.11
"""Figure 1 for the Angewandte manuscript: LCC three-step schematic (a) + shift-tolerance curve (b).

ponytail: pure-matplotlib schematic, no repo imports; regenerate with `python3.11 make_schematic.py`.
"""
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

OUT = Path(__file__).with_name("lcc_schematic.png")

rng = np.random.default_rng(0)


def _peaks(centers, grid, s2, s1):
    """Sum of 2D Gaussians on a (F1, F2) grid."""
    g2, g1 = np.meshgrid(grid[1], grid[0])
    img = np.zeros_like(g2)
    for (c2, c1, amp) in centers:
        img += amp * np.exp(-((g2 - c2) ** 2) / (2 * s2 ** 2) - ((g1 - c1) ** 2) / (2 * s1 ** 2))
    return img


def main() -> None:
    fig = plt.figure(figsize=(9.2, 3.4))
    gs = fig.add_gridspec(1, 5, width_ratios=[1.0, 1.0, 1.0, 0.28, 1.25], wspace=0.35)

    ax_render = fig.add_subplot(gs[0, 0])
    ax_blur = fig.add_subplot(gs[0, 1])
    ax_score = fig.add_subplot(gs[0, 2])
    ax_curve = fig.add_subplot(gs[0, 4])

    f2 = np.linspace(6.0, 9.5, 160)
    f1 = np.linspace(105, 128, 160)
    grid = (f1, f2)
    centers = [(7.1, 112, 1.0), (8.3, 120, 0.85), (7.8, 124, 0.7),
               (6.6, 108, 0.6), (8.9, 116, 0.9), (7.4, 118, 0.75)]

    # (1) Render: sharp points on shared grid
    sharp = _peaks(centers, grid, 0.02, 0.2)
    ax_render.imshow(sharp, origin="upper", aspect="auto", cmap="viridis",
                     extent=[f2[0], f2[-1], f1[-1], f1[0]])
    ax_render.set_title("1. Render\n(shared grid, area-weighted)", fontsize=9)

    # (2) Blur: physical linewidth
    blurred = _peaks(centers, grid, 0.09, 0.9)
    ax_blur.imshow(blurred, origin="upper", aspect="auto", cmap="viridis",
                   extent=[f2[0], f2[-1], f1[-1], f1[0]])
    ax_blur.set_title(r"2. Blur" "\n" r"($\sigma=\sqrt{\ell^2+d^2}$)", fontsize=9)

    # (3) Score: mean-centred overlay of two slightly-drifted copies
    drifted = _peaks([(c2 + 0.05, c1 + 0.5, a) for c2, c1, a in centers], grid, 0.09, 0.9)
    overlay = (blurred - blurred.mean()) * (drifted - drifted.mean())
    vmax = np.abs(overlay).max()
    ax_score.imshow(overlay, origin="upper", aspect="auto", cmap="RdBu_r",
                    vmin=-vmax, vmax=vmax, extent=[f2[0], f2[-1], f1[-1], f1[0]])
    ax_score.set_title("3. Score\n(mean-centred ZNCC)", fontsize=9)

    for ax in (ax_render, ax_blur, ax_score):
        ax.set_xlabel(r"$^1$H (ppm)", fontsize=8)
        ax.set_xticks([6, 8]); ax.set_yticks([110, 120])
        ax.tick_params(labelsize=7)
    ax_render.set_ylabel(r"$^{15}$N (ppm)", fontsize=8)

    # arrows between panels
    for a, b in ((ax_render, ax_blur), (ax_blur, ax_score)):
        arr = FancyArrowPatch((1.04, 0.5), (1.20, 0.5), transform=a.transAxes,
                              mutation_scale=13, color="0.35", clip_on=False,
                              arrowstyle="-|>", lw=1.4)
        fig.add_artist(arr)

    fig.text(0.02, 0.96, "(a)", fontsize=12, fontweight="bold")
    fig.text(0.70, 0.96, "(b)", fontsize=12, fontweight="bold")

    # (b) shift-tolerance decay exp(-Δ²/4σ²)
    d = np.linspace(0, 3, 200)
    ax_curve.plot(d, np.exp(-d ** 2 / 4), color="#2c6fbb", lw=2.2)
    ax_curve.axvspan(0, 0.6, color="#7fc97f", alpha=0.25)
    ax_curve.text(0.28, 0.12, "physical\ndrift", fontsize=8, ha="center", color="#2a7d2a")
    ax_curve.text(2.0, 0.55, "random\nrelocation", fontsize=8, ha="center", color="#b03030")
    ax_curve.annotate("", xy=(2.6, 0.06), xytext=(2.6, 0.42),
                      arrowprops=dict(arrowstyle="->", color="#b03030"))
    ax_curve.set_xlabel(r"displacement $\Delta/\sigma$", fontsize=9)
    ax_curve.set_ylabel(r"correlation $\exp(-\Delta^2/4\sigma^2)$", fontsize=9)
    ax_curve.set_title("Graded shift tolerance", fontsize=9)
    ax_curve.set_xlim(0, 3); ax_curve.set_ylim(0, 1.02)
    ax_curve.tick_params(labelsize=8)
    ax_curve.grid(alpha=0.25)

    fig.savefig(OUT, dpi=200, bbox_inches="tight")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
