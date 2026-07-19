#!/usr/bin/env python3.11
"""Method-explainer figures for the Supporting Information (real STCC math on toy spectra).

Ten figures, numbered by their position in the SI (S1..S10):

- si_baseline.png      (S1, §S1)  : clip vs abs on an edited HSQC -- clip deletes CH2, abs keeps it.
- si_sharedgrid.png    (S2, §S2)  : two native resolutions rendered onto one shared grid.
- si_blur.png          (S3, §S3)  : blur width sigma = linewidth (+) drift in quadrature.
- si_pipeline.png      (S4, §S4)  : the four STCC steps end-to-end (render->blur->centre->score).
- si_failuremodes.png  (S5, §S6)  : brittle (bin) vs blind (tree/NN) vs clean gap (STCC), real numbers.
- si_shifttol.png      (S6, §S6a) : graded (STCC) vs brittle (bin) response to a shift.
- si_sigma.png         (S7, §S6a) : sigma is the single tolerance knob.
- si_meancentre.png    (S8, §S6b) : mean-centring turns a non-co-located peak into a negative product.
- si_zerolag.png       (S9, §S6c) : zero lag vs a global shift search (false registration).
- si_gridinvariance.png(S10, §S7) : the score is grid-invariant while >= 3 px per sigma.

Run: python3.11 make_si_figs.py
"""
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

HERE = Path(__file__).parent
VIR = "viridis"
DIV = "RdBu_r"     # red = positive, blue = negative
BLUE, ORANGE, RED, GREY = "#2c6fbb", "#e08a1e", "#b03030", "#9aa0a6"


# --------------------------------------------------------------------------- #
# shared toy-spectrum helpers (a tiny, honest version of the real pipeline)
# --------------------------------------------------------------------------- #
def render(peaks, gx, gy):
    """Histogram point peaks onto the shared grid (area-weighted = unit spikes here)."""
    img = np.zeros((gy.size, gx.size))
    for x, y, a in peaks:
        i = int(np.argmin(np.abs(gy - y)))
        j = int(np.argmin(np.abs(gx - x)))
        img[i, j] += a
    return img


def blur(img, sx, sy):
    from scipy.ndimage import gaussian_filter  # noqa
    return gaussian_filter(img, (sy, sx))


def blur_nofft(img, sx, sy):
    """Separable Gaussian without scipy (kernels in pixels)."""
    def k(s):
        r = max(1, int(round(3 * s)))
        x = np.arange(-r, r + 1)
        w = np.exp(-0.5 * (x / s) ** 2)
        return w / w.sum()
    out = np.apply_along_axis(lambda r: np.convolve(r, k(sx), "same"), 1, img)
    out = np.apply_along_axis(lambda c: np.convolve(c, k(sy), "same"), 0, out)
    return out


def centre(img):
    return img - img.mean()


# --------------------------------------------------------------------------- #
# Figure S1 : the pipeline end to end
# --------------------------------------------------------------------------- #
def fig_pipeline():
    gx = np.linspace(6, 9.5, 100)
    gy = np.linspace(108, 128, 100)
    X = [(7.0, 112, 1.0), (8.2, 120, 0.8), (7.6, 124, 0.7), (6.6, 116, 0.6)]
    # Y = same compound, small titration drift (< blur width -> high score)
    Y = [(x + 0.05, y + 0.6, a) for x, y, a in X]

    rx, ry = render(X, gx, gy), render(Y, gx, gy)
    bx, by = blur_nofft(rx, 3.5, 3.5), blur_nofft(ry, 3.5, 3.5)
    cx, cy = centre(bx), centre(by)
    prod = cx * cy
    score = float(np.sum(cx * cy) / np.sqrt(np.sum(cx**2) * np.sum(cy**2)))

    fig = plt.figure(figsize=(10.4, 5.0))
    gs = fig.add_gridspec(2, 5, width_ratios=[1, 1, 1, 0.22, 1.4], hspace=0.3, wspace=0.32)
    ext = [gx[0], gx[-1], gy[-1], gy[0]]

    def scatter(ax, peaks, col):
        for x, y, a in peaks:
            ax.scatter([x], [y], s=90 * a, color=col, edgecolor="white", linewidth=0.8, zorder=3)
        ax.set_xlim(gx[0], gx[-1]); ax.set_ylim(gy[-1], gy[0])
        ax.set_facecolor("#f3f4f6")

    def show(ax, im, cmap, sym=False):
        v = np.abs(im).max() or 1
        kw = dict(vmin=-v, vmax=v) if sym else {}
        ax.imshow(im, origin="upper", aspect="auto", cmap=cmap, extent=ext, **kw)

    axes = [[fig.add_subplot(gs[r, c]) for c in range(3)] for r in range(2)]
    scatter(axes[0][0], X, BLUE); show(axes[0][1], bx, VIR); show(axes[0][2], cx, DIV, sym=True)
    scatter(axes[1][0], Y, ORANGE); show(axes[1][1], by, VIR); show(axes[1][2], cy, DIV, sym=True)
    for r in range(2):
        for c in range(3):
            axes[r][c].set_xticks([6, 8]); axes[r][c].set_yticks([110, 120])
            axes[r][c].tick_params(labelsize=6.5)
    axes[0][0].set_title("1. Render", fontsize=9)
    axes[0][1].set_title("2. Blur", fontsize=9)
    axes[0][2].set_title("3. Mean-centre", fontsize=9)
    axes[0][0].set_ylabel("reference $X$", fontsize=8.5)
    axes[1][0].set_ylabel("same compound,\nshifted $Y$", fontsize=8.5)

    axp = fig.add_subplot(gs[:, 4])
    v = np.abs(prod).max() or 1
    axp.imshow(prod, origin="upper", aspect="auto", cmap=DIV, vmin=-v, vmax=v, extent=ext)
    axp.set_title(r"4. Score $= \widehat{G}_X\!\cdot\widehat{G}_Y$ (summed)", fontsize=9.5)
    axp.set_xticks([6, 8]); axp.set_yticks([110, 120]); axp.tick_params(labelsize=6.5)
    axp.text(0.5, -0.12, f"$S = {score:.2f}$", transform=axp.transAxes, ha="center",
             fontsize=12, color=RED, fontweight="bold")
    axp.text(0.5, -0.2, "red $+$ co-located  ·  blue $-$ mismatch", transform=axp.transAxes,
             ha="center", fontsize=7.5, color="#555")

    for x0 in (0.318, 0.516):
        fig.add_artist(FancyArrowPatch((x0, 0.70), (x0 + 0.028, 0.70), transform=fig.transFigure,
                                       mutation_scale=12, color="0.4", arrowstyle="-|>"))
        fig.add_artist(FancyArrowPatch((x0, 0.30), (x0 + 0.028, 0.30), transform=fig.transFigure,
                                       mutation_scale=12, color="0.4", arrowstyle="-|>"))
    fig.add_artist(FancyArrowPatch((0.635, 0.5), (0.70, 0.5), transform=fig.transFigure,
                                   mutation_scale=15, color=RED, arrowstyle="-|>"))
    fig.suptitle("The STCC pipeline: render → blur → mean-centre → correlate at zero lag",
                 fontsize=11.5, y=1.0)
    fig.savefig(HERE / "si_pipeline.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote si_pipeline.png  (toy score =", round(score, 3), ")")


# --------------------------------------------------------------------------- #
# Figure S2 : mean-centring = coverage penalty (1D cross-section)
# --------------------------------------------------------------------------- #
def fig_meancentre():
    x = np.linspace(0, 10, 400)
    g = lambda c, h=1.0, w=0.35: h * np.exp(-0.5 * ((x - c) / w) ** 2)
    # X has peaks at 3 and 7; Y has a peak at 3 (co-located) and at 5 (X empty there)
    X = g(3.0) + g(7.0)
    Y = g(3.05) + g(5.0)

    fig, axes = plt.subplots(3, 2, figsize=(9.4, 5.0), sharex=True)
    for col, (title, cen) in enumerate([("STCC: mean-centred", True), ("cosine: un-centred", False)]):
        Xh = X - X.mean() if cen else X
        Yh = Y - Y.mean() if cen else Y
        prod = Xh * Yh
        axes[0, col].plot(x, Xh, color=BLUE); axes[0, col].set_title(title, fontsize=10)
        axes[1, col].plot(x, Yh, color=ORANGE)
        axes[2, col].plot(x, prod, color="0.3", lw=1)
        axes[2, col].fill_between(x, prod, 0, where=prod >= 0, color=RED, alpha=0.5)
        axes[2, col].fill_between(x, prod, 0, where=prod < 0, color=BLUE, alpha=0.5)
        for r, lab in zip(range(3), [r"$\widehat{G}_X$", r"$\widehat{G}_Y$", r"product"]):
            axes[r, col].axhline(0, color="0.7", lw=0.6)
            if col == 0:
                axes[r, col].set_ylabel(lab if cen else lab, fontsize=9)
        # annotations
        if cen:
            axes[2, col].annotate("co-located\n→ $+$ (reward)", xy=(3.0, prod[np.argmin(abs(x-3))]),
                                  xytext=(1.2, 0.55), fontsize=7.5, color=RED,
                                  arrowprops=dict(arrowstyle="->", color=RED))
            axes[2, col].annotate("$X$ peak, $Y$ empty\n→ $-$ (penalty)", xy=(7.0, prod[np.argmin(abs(x-7))]),
                                  xytext=(6.2, -0.5), fontsize=7.5, color=BLUE,
                                  arrowprops=dict(arrowstyle="->", color=BLUE))
        else:
            axes[2, col].text(5, 0.18, "all products $\\geq 0$:\nmismatch not penalised",
                              fontsize=7.5, color="#555", ha="center")
    axes[2, 0].set_xlabel("chemical shift (a.u.)", fontsize=9)
    axes[2, 1].set_xlabel("chemical shift (a.u.)", fontsize=9)
    fig.suptitle("Mean-centring is the discriminating step: it turns a non-co-located peak into a negative product",
                 fontsize=10.5, y=1.0)
    fig.tight_layout()
    fig.savefig(HERE / "si_meancentre.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote si_meancentre.png")


# --------------------------------------------------------------------------- #
# Figure S3 : graded (STCC) vs brittle (bin) response to a shift
# --------------------------------------------------------------------------- #
def fig_shifttol():
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(9.4, 3.8), gridspec_kw=dict(width_ratios=[1, 1.15]))

    # (a) mechanism: a peak near a bin edge
    x = np.linspace(0, 4, 400)
    for ax_pos, c, lab in [(0.0, "#2c6fbb", "peak inside bin"), (0.0, None, None)]:
        pass
    peak = lambda c: np.exp(-0.5 * ((x - c) / 0.25) ** 2)
    edge = 2.0
    ax0.fill_between(x, peak(1.7), color=BLUE, alpha=0.35, label="peak")
    ax0.plot(x, peak(1.7), color=BLUE)
    for e in (1.0, 2.0, 3.0):
        ax0.axvline(e, color=GREY, ls="--", lw=1)
    ax0.text(1.0, 1.06, "bin edges", fontsize=7.5, color=GREY, ha="center")
    ax0.annotate("a small drift splits the\npeak's mass across the edge",
                 xy=(2.0, 0.45), xytext=(2.2, 0.82), fontsize=7.5, color="#555",
                 arrowprops=dict(arrowstyle="->", color="#555"))
    ax0.set_title("(a) Hard bin edges are brittle", fontsize=10, loc="left")
    ax0.set_ylim(0, 1.18); ax0.set_yticks([]); ax0.set_xlabel("chemical shift (a.u.)", fontsize=9)

    # (b) similarity vs displacement: smooth STCC vs discontinuous bin
    d = np.linspace(0, 1.5, 400)
    sigma = 0.4
    lcc = np.exp(-(d ** 2) / (4 * sigma ** 2))
    # bin: MONOTONE staircase -- similarity erodes slowly inside a bin, then drops discontinuously
    # each time the peak crosses an edge (mass leaves the matched bin). No recovery.
    binw = 0.5
    step = np.floor(d / binw)
    binsim = np.clip(0.95 - 0.10 * (d - step * binw) / binw - 0.33 * step, 0, 1)
    ax1.plot(d, lcc, color=BLUE, lw=2.4, label="STCC (Gaussian blur)")
    for seg in range(3):
        m = (d >= seg * binw) & (d < (seg + 1) * binw)
        ax1.plot(d[m], binsim[m], color=ORANGE, lw=2.2)
    ax1.plot([], [], color=ORANGE, lw=2.2, label="bin method (hard edges)")
    for e in (binw, 2 * binw):
        lo = np.clip(0.95 - 0.10 - 0.33 * (e / binw - 1), 0, 1)
        hi = np.clip(0.95 - 0.33 * (e / binw), 0, 1)
        ax1.annotate("", xy=(e, hi), xytext=(e, lo),
                     arrowprops=dict(arrowstyle="-", color=ORANGE, ls=":", lw=1.1))
    ax1.text(binw + 0.03, 0.72, "discontinuous\ndrop at each edge", fontsize=7.2, color=ORANGE, ha="left")
    ax1.set_title("(b) Similarity vs shift: graded vs brittle", fontsize=10, loc="left")
    ax1.set_xlabel(r"peak displacement $\Delta$ (a.u.)", fontsize=9)
    ax1.set_ylabel("similarity contribution", fontsize=9)
    ax1.set_ylim(-0.02, 1.05); ax1.legend(fontsize=8, loc="upper right")
    ax1.grid(alpha=0.25)

    fig.suptitle("Graded shift tolerance: STCC decays smoothly where the bin method jumps",
                 fontsize=10.5, y=1.0)
    fig.tight_layout()
    fig.savefig(HERE / "si_shifttol.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote si_shifttol.png")


# --------------------------------------------------------------------------- #
# Figure S1 (§S1) : baseline modes on a multiplicity-edited HSQC cross-section
# --------------------------------------------------------------------------- #
def fig_baseline():
    x = np.linspace(0, 9, 500)
    g = lambda c, h, w=0.5: h * np.exp(-0.5 * ((x - c) / w) ** 2)
    I = g(2.2, 1.0) + g(6.4, -0.8)  # CH/CH3 positive, CH2 negative

    fig, axes = plt.subplots(1, 3, figsize=(10.0, 3.2), sharey=True)
    panels = [
        ("raw edited HSQC", I, None),
        ("clip (default): $\\max(I,0)$", np.clip(I, 0, None), "CH$_2$ peak\ndeleted"),
        ("abs (edited): $|I|$", np.abs(I), "CH$_2$ peak\npreserved"),
    ]
    for ax, (title, y, note) in zip(axes, panels):
        ax.fill_between(x, y, 0, where=y >= 0, color=RED, alpha=0.45)
        ax.fill_between(x, y, 0, where=y < 0, color=BLUE, alpha=0.45)
        ax.plot(x, y, color="0.3", lw=1)
        ax.axhline(0, color="0.6", lw=0.7)
        ax.set_title(title, fontsize=9.5)
        ax.set_xlabel("chemical shift (a.u.)", fontsize=8.5)
        ax.tick_params(labelsize=7)
        ax.set_ylim(-1.05, 1.15)
    axes[0].set_ylabel("intensity", fontsize=9)
    axes[0].annotate("CH / CH$_3$ ($+$)", xy=(2.2, 1.0), xytext=(2.5, 0.55),
                     fontsize=7.5, color=RED)
    axes[0].annotate("CH$_2$ ($-$)", xy=(6.4, -0.8), xytext=(3.5, -0.9),
                     fontsize=7.5, color=BLUE, arrowprops=dict(arrowstyle="->", color=BLUE, lw=0.8))
    axes[1].annotate(panels[1][2], xy=(6.4, 0.0), xytext=(5.4, 0.55), fontsize=7.5,
                     color=RED, ha="center", arrowprops=dict(arrowstyle="->", color=RED, lw=0.8))
    axes[2].annotate(panels[2][2], xy=(6.4, 0.8), xytext=(5.2, 0.42), fontsize=7.5,
                     color="#2a7", ha="center", arrowprops=dict(arrowstyle="->", color="#2a7", lw=0.8))
    fig.suptitle("Baseline modes: clip is correct for noise, but deletes CH$_2$ peaks in edited HSQC — use abs there",
                 fontsize=10.5, y=1.02)
    fig.tight_layout()
    fig.savefig(HERE / "si_baseline.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote si_baseline.png")


# --------------------------------------------------------------------------- #
# Figure S2 (§S2) : two native resolutions rendered onto one shared grid
# --------------------------------------------------------------------------- #
def fig_sharedgrid():
    # same four peaks, sampled by two instruments at different native step sizes
    peaks = [(1.6, 6.2), (3.1, 3.4), (5.0, 8.1), (6.8, 5.0)]
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.5))

    def grid(ax, step, color, label, marker):
        ticks = np.arange(0, 8.001, step)
        for t in ticks:
            ax.axvline(t, color="0.85", lw=0.5)
            ax.axhline(t, color="0.85", lw=0.5)
        for (px, py) in peaks:
            ax.scatter([px], [py], s=70, color=color, edgecolor="white", linewidth=0.7, zorder=3)
        ax.set_xlim(0, 8); ax.set_ylim(0, 9); ax.set_aspect("auto")
        ax.set_title(label, fontsize=9.5)
        ax.tick_params(labelsize=7)

    grid(axes[0], 0.4, BLUE, "Spectrum $X$\n(native $\\Delta$ fine)", "o")
    grid(axes[1], 1.0, ORANGE, "Spectrum $Y$\n(native $\\Delta$ coarse)", "s")

    # shared grid: both rendered onto identical edges -> directly comparable
    step = 0.8
    ticks = np.arange(0, 8.001, step)
    ax = axes[2]
    for t in ticks:
        ax.axvline(t, color="0.75", lw=0.6)
        ax.axhline(t, color="0.75", lw=0.6)
    for (px, py), col, dx in [((p[0], p[1]), BLUE, -0.12) for p in peaks] + \
                             [((p[0] + 0.05, p[1] + 0.05), ORANGE, 0.12) for p in peaks]:
        # snap to shared cell centre to show both land in the SAME cell
        cx = (np.floor(px / step) + 0.5) * step
        cy = (np.floor(py / step) + 0.5) * step
        ax.scatter([cx + dx], [cy], s=70, color=col, edgecolor="white", linewidth=0.7, zorder=3)
    ax.set_xlim(0, 8); ax.set_ylim(0, 9)
    ax.set_title("Shared grid\n($X$ ●, $Y$ ● on same edges)", fontsize=9.5)
    ax.tick_params(labelsize=7)
    ax.annotate("both land in the\nsame cell → comparable", xy=(3.2, 3.6), xytext=(3.4, 7.6),
                fontsize=7.2, color="#333", ha="center",
                arrowprops=dict(arrowstyle="->", color="#333", lw=0.8))
    for a in axes:
        a.set_xlabel("$F2$ (a.u.)", fontsize=8.5)
    axes[0].set_ylabel("$F1$ (a.u.)", fontsize=8.5)
    fig.suptitle("Rendering onto shared edges makes two spectra comparable regardless of native resolution",
                 fontsize=10.5, y=1.02)
    fig.tight_layout()
    fig.savefig(HERE / "si_sharedgrid.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote si_sharedgrid.png")


# --------------------------------------------------------------------------- #
# Figure S3 (§S3) : blur width = linewidth (+) drift in quadrature
# --------------------------------------------------------------------------- #
def fig_blur():
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(9.6, 3.8), gridspec_kw=dict(width_ratios=[1, 1.2]))

    # (a) quadrature right-triangle: sigma = sqrt(l^2 + d^2)
    l, d = 3.0, 2.0
    axa.plot([0, l], [0, 0], color=BLUE, lw=2.5)
    axa.plot([l, l], [0, d], color=ORANGE, lw=2.5)
    axa.plot([0, l], [0, d], color=RED, lw=2.5)
    # right-angle mark
    axa.plot([l - 0.35, l - 0.35, l], [0, 0.35, 0.35], color="0.5", lw=0.9)
    axa.text(l / 2, -0.28, "linewidth $\\ell$", color=BLUE, fontsize=9.5, ha="center")
    axa.text(l + 0.12, d / 2, "drift $d$", color=ORANGE, fontsize=9.5, va="center")
    axa.text(l / 2 - 0.35, d / 2 + 0.3, "$\\sigma=\\sqrt{\\ell^2+d^2}$", color=RED,
             fontsize=11, ha="center", rotation=np.degrees(np.arctan2(d, l)))
    axa.set_xlim(-0.5, l + 1.2); axa.set_ylim(-0.7, d + 0.7)
    axa.set_aspect("equal"); axa.axis("off")
    axa.set_title("(a) Blur width combines two physical widths", fontsize=9.5, loc="left")

    # (b) a rendered spike blurred by two sigma -> larger sigma = wider tolerance
    x = np.linspace(-4, 4, 400)
    axb.stem([0], [1.0], linefmt="0.5", markerfmt=" ", basefmt=" ")
    axb.text(0.1, 1.02, "rendered\npeak (spike)", fontsize=7.5, color="0.4")
    for sig, col, lab in [(0.6, BLUE, "$\\sigma=\\ell$ (linewidth only)"),
                          (1.3, RED, "$\\sigma=\\sqrt{\\ell^2+d^2}$ (with drift)")]:
        y = np.exp(-0.5 * (x / sig) ** 2)
        axb.plot(x, y, color=col, lw=2.2, label=lab)
        axb.fill_between(x, y, 0, color=col, alpha=0.12)
    axb.annotate("larger $\\sigma$ →\nwider shift tolerance", xy=(1.3, 0.6), xytext=(1.7, 0.85),
                 fontsize=7.5, color=RED, arrowprops=dict(arrowstyle="->", color=RED, lw=0.8))
    axb.set_xlabel("chemical shift (a.u.)", fontsize=9)
    axb.set_ylabel("blurred intensity", fontsize=9)
    axb.set_ylim(0, 1.15); axb.legend(fontsize=8, loc="upper left")
    axb.set_title("(b) The blur sets the shift-tolerance scale", fontsize=9.5, loc="left")
    fig.suptitle("Lineshape blur $\\sigma$ = physical linewidth and expected drift added in quadrature",
                 fontsize=10.5, y=1.0)
    fig.tight_layout()
    fig.savefig(HERE / "si_blur.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote si_blur.png")


# --------------------------------------------------------------------------- #
# Figure S5 (§S6) : the two opposite failure modes (real Table S2 numbers)
# --------------------------------------------------------------------------- #
def fig_failuremodes():
    # dense 1H-15N benchmark, Table S2 (self=1 excluded from the same-range)
    same = {
        "Bin\n(Bodis)":       [0.8106, 0.7739, 0.8274, 0.8164, 0.7637, 0.7320],
        "Tree\n(Castillo)":   [0.9401, 0.8619, 0.9443, 0.9179, 0.8760, 0.8624],
        "NN\n(Pierens)":      [0.9925, 0.9981, 0.9877, 0.9848, 0.9978, 0.9851],
        "STCC\n(this work)":   [0.9662, 0.9594, 0.9485, 0.8972, 0.9558, 0.8911],
    }
    diff = {"Bin\n(Bodis)": 0.4949, "Tree\n(Castillo)": 0.8721,
            "NN\n(Pierens)": 0.9550, "STCC\n(this work)": 0.1815}
    tags = {"Bin\n(Bodis)": ("brittle", ORANGE), "Tree\n(Castillo)": ("blind", GREY),
            "NN\n(Pierens)": ("blind", GREY), "STCC\n(this work)": ("clean gap", BLUE)}

    fig, ax = plt.subplots(figsize=(8.4, 4.4))
    xs = np.arange(len(same))
    for i, m in enumerate(same):
        lo, hi = min(same[m]), max(same[m])
        col = BLUE if m.startswith("STCC") else "0.55"
        # same-compound range as a vertical band
        ax.add_patch(plt.Rectangle((i - 0.16, lo), 0.32, hi - lo, color=col, alpha=0.35, zorder=2))
        ax.plot([i - 0.16, i + 0.16], [np.mean(same[m])] * 2, color=col, lw=1.5, zorder=3)
        # different-protein marker
        ax.scatter([i], [diff[m]], marker="v", s=90, color=RED, edgecolor="white",
                   linewidth=0.8, zorder=4)
        # gap arrow + tag
        tag, tcol = tags[m]
        ax.annotate("", xy=(i + 0.28, diff[m]), xytext=(i + 0.28, lo),
                    arrowprops=dict(arrowstyle="<->", color=tcol, lw=1.1))
        ax.text(i + 0.33, (lo + diff[m]) / 2, f"{lo - diff[m]:+.2f}", color=tcol,
                fontsize=7.5, va="center")
        ax.text(i, 1.05, tag, color=tcol, fontsize=8.5, ha="center", fontweight="bold")
    ax.scatter([], [], marker="v", s=90, color=RED, edgecolor="white", label="different protein")
    ax.add_patch(plt.Rectangle((0, 0), 0, 0, color="0.55", alpha=0.35, label="same-protein range"))
    ax.set_xticks(xs); ax.set_xticklabels(list(same), fontsize=9)
    ax.set_ylim(0, 1.13); ax.set_ylabel("similarity", fontsize=9.5)
    ax.legend(fontsize=8, loc="lower left", frameon=True)
    ax.grid(axis="y", alpha=0.25)
    ax.set_title("Two opposite failures: the bin method is shift-brittle, tree/NN are shift-blind; "
                 "only STCC pushes\nthe different protein below every same-protein score",
                 fontsize=9.8)
    fig.tight_layout()
    fig.savefig(HERE / "si_failuremodes.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote si_failuremodes.png")


# --------------------------------------------------------------------------- #
# Figure S7 (§S6a) : sigma is the single tolerance knob
# --------------------------------------------------------------------------- #
def fig_sigma():
    d = np.linspace(0, 1.2, 400)
    fig, ax = plt.subplots(figsize=(7.6, 3.8))
    cols = [BLUE, "#5a8fcf", "#8bb4e0"]
    for sig, col in zip([0.15, 0.30, 0.60], cols):
        ax.plot(d, np.exp(-(d ** 2) / (4 * sig ** 2)), color=col, lw=2.4,
                label=f"$\\sigma={sig:.2f}$ ppm")
    # typical titration drift band
    ax.axvspan(0.2, 0.8, color=ORANGE, alpha=0.12)
    ax.text(0.5, 0.06, "typical $^{15}$N\ntitration drift", color="#b06a10", fontsize=7.5, ha="center")
    ax.set_xlabel(r"peak displacement $\Delta_{F1}$ (ppm)", fontsize=9.5)
    ax.set_ylabel(r"similarity contribution $\exp(-\Delta^2/4\sigma^2)$", fontsize=9.5)
    ax.set_ylim(0, 1.03); ax.legend(fontsize=8.5, title="single knob", title_fontsize=8.5)
    ax.grid(alpha=0.25)
    ax.set_title("One knob $\\sigma$ sets the shift-tolerance scale: small $\\sigma$ is strict, large $\\sigma$ is forgiving",
                 fontsize=9.8)
    fig.tight_layout()
    fig.savefig(HERE / "si_sigma.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote si_sigma.png")


# --------------------------------------------------------------------------- #
# Figure S9 (§S6c) : zero lag on purpose vs a global shift search
# --------------------------------------------------------------------------- #
def fig_zerolag():
    x = np.linspace(0, 10, 500)
    g = lambda c, h=1.0, w=0.3: h * np.exp(-0.5 * ((x - c) / w) ** 2)
    # genuinely different patterns (unequal spacing -> no rigid shift matches all), but the best
    # lag still aligns a subset and inflates similarity above the true zero-lag value
    A = g(2.0) + g(4.0) + g(6.0)
    B = g(2.6) + g(5.0, 0.9) + g(6.4, 0.85)
    Ah, Bh = A - A.mean(), B - B.mean()
    corr = np.correlate(Ah, Bh, "full") / np.sqrt((Ah @ Ah) * (Bh @ Bh))
    dx = x[1] - x[0]
    lags = (np.arange(corr.size) - (Bh.size - 1)) * dx
    zero_i = np.argmin(np.abs(lags))
    gmax_i = int(np.argmax(corr))

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(9.8, 3.8), gridspec_kw=dict(width_ratios=[1, 1.1]))
    ax0.plot(x, A, color=BLUE, lw=1.8, label="spectrum $A$")
    ax0.plot(x, B, color=ORANGE, lw=1.8, label="spectrum $B$ (different compound)")
    ax0.set_title("(a) Two different compounds", fontsize=9.5, loc="left")
    ax0.set_xlabel("chemical shift (a.u.)", fontsize=9)
    ax0.set_yticks([]); ax0.legend(fontsize=7.8, loc="upper right")

    ax1.plot(lags, corr, color="0.4", lw=1.8)
    ax1.axvline(0, color="0.7", lw=0.8, ls="--")
    ax1.scatter([0], [corr[zero_i]], color=BLUE, s=70, zorder=4, edgecolor="white")
    ax1.scatter([lags[gmax_i]], [corr[gmax_i]], color=RED, s=70, zorder=4, edgecolor="white")
    ax1.annotate(f"zero lag → {corr[zero_i]:.2f}\n(true: low)", xy=(0, corr[zero_i]),
                 xytext=(-4.6, corr[zero_i] + 0.28), fontsize=7.6, color=BLUE,
                 arrowprops=dict(arrowstyle="->", color=BLUE, lw=0.8))
    ax1.annotate(f"shift search →\n{corr[gmax_i]:.2f} (false high)", xy=(lags[gmax_i], corr[gmax_i]),
                 xytext=(lags[gmax_i] - 1.0, corr[gmax_i] + 0.12), fontsize=7.6, color=RED,
                 arrowprops=dict(arrowstyle="->", color=RED, lw=0.8))
    ax1.set_title("(b) Correlation vs lag", fontsize=9.5, loc="left")
    ax1.set_xlabel("lag $\\tau$ (a.u.)", fontsize=9)
    ax1.set_ylabel("normalized correlation", fontsize=9)
    fig.suptitle("Scoring at zero lag keeps position information; a global shift search would slide different spectra into false registration",
                 fontsize=9.8, y=1.02)
    fig.tight_layout()
    fig.savefig(HERE / "si_zerolag.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote si_zerolag.png  (zero-lag =", round(float(corr[zero_i]), 3),
          " global-max =", round(float(corr[gmax_i]), 3), ")")


# --------------------------------------------------------------------------- #
# Figure S10 (§S7) : score is grid-invariant while >= 3 px per sigma
# --------------------------------------------------------------------------- #
def fig_gridinvariance():
    # toy same-compound pair (as in the pipeline), scored across render pixel sizes
    Xp = [(7.0, 112, 1.0), (8.2, 120, 0.8), (7.6, 124, 0.7), (6.6, 116, 0.6)]
    Yp = [(x + 0.05, y + 0.6, a) for x, y, a in Xp]
    sig_x, sig_y = 0.30, 3.0  # blur width in ppm (F2, F1)

    steps = np.linspace(0.05, 1.3, 22)  # F1 pixel size in ppm; F2 kept 1/10 of it
    px_per_sigma, scores = [], []
    for st in steps:
        gx = np.arange(6, 9.5, st / 10.0)
        gy = np.arange(108, 128, st)
        rx, ry = render(Xp, gx, gy), render(Yp, gx, gy)
        sxp, syp = sig_x / (st / 10.0), sig_y / st  # blur width in pixels
        bx, by = blur_nofft(rx, sxp, syp), blur_nofft(ry, sxp, syp)
        cx, cy = centre(bx), centre(by)
        den = np.sqrt(np.sum(cx ** 2) * np.sum(cy ** 2))
        scores.append(float(np.sum(cx * cy) / den) if den else 0.0)
        px_per_sigma.append(sig_y / st)  # F1 axis is the coarse one

    px_per_sigma = np.array(px_per_sigma); scores = np.array(scores)
    order = np.argsort(px_per_sigma)
    px_per_sigma, scores = px_per_sigma[order], scores[order]

    fig, ax = plt.subplots(figsize=(7.4, 2.1))
    ax.plot(px_per_sigma, scores, "-o", color=BLUE, ms=4, lw=1.8)
    ax.axvspan(0, 3, color=RED, alpha=0.10)
    ax.axvline(3, color=RED, lw=1.2, ls="--")
    ax.text(2.9, 0.2, "under-sampled\n($<3$ px per $\\sigma$)", color=RED, fontsize=7.6, ha="right")
    ax.text(px_per_sigma.max() * 0.6, min(scores.max(), 0.99) - 0.06,
            "score flat: grid-invariant", color=BLUE, fontsize=8, ha="center")
    ax.set_xlabel("pixels per $\\sigma$ (render resolution)", fontsize=9.5)
    ax.set_ylabel("STCC score (fixed toy pair)", fontsize=9.5)
    ax.set_ylim(0, 1.02); ax.grid(alpha=0.25)
    ax.set_title("The render grid is a discretization, not a parameter: the score is flat while $\\geq 3$ px per $\\sigma$",
                 fontsize=9.5)
    fig.tight_layout()
    fig.savefig(HERE / "si_gridinvariance.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote si_gridinvariance.png")


if __name__ == "__main__":
    fig_baseline()       # S1  (§S1)
    fig_sharedgrid()     # S2  (§S2)
    fig_blur()           # S3  (§S3)
    fig_pipeline()       # S4  (§S4)
    fig_failuremodes()   # S5  (§S6)
    fig_shifttol()       # S6  (§S6a)
    fig_sigma()          # S7  (§S6a)
    fig_meancentre()     # S8  (§S6b)
    fig_zerolag()        # S9  (§S6c)
    fig_gridinvariance() # S10 (§S7)
