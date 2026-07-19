#!/usr/bin/env python3.11
"""Figure 3: sparse cluster-bootstrap CIs + the mean-centring ablation by regime.

(a) Sparse 1H-13C separation with compound-level bootstrap 95% CIs (from retrieval_13c.json),
    including the experimental Local-Contrast candidate.
(b) Mean-centring ablation: un-centred cosine vs STCC in each regime -- a large gap on the dense
    protein fingerprint (0.59->0.75), almost none on the sparse stick set (0.7445≈0.7447).

Colour by entity-class (CVD-safe, consistent with Figs 1-2): proposed = blue, its ablation =
light blue, strong baseline (bin) = orange, saturating (tree/NN) = grey. Run: python3.11 make_fig3.py
"""
import json
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).parent
BLUE, LBLUE, PURPLE, ORANGE, GREY = "#2c6fbb", "#6aa0d6", "#7b4ab0", "#e08a1e", "#9aa0a6"

# (a) sparse separation + 95% CI, from the retrieval/bootstrap harness
R = json.load(open(HERE / "retrieval_13c.json"))["rows"]
ROWS = [  # (label, key, colour) sorted high->low separation
    ("Local Contrast (experimental)", "local_contrast", PURPLE),
    ("STCC (this work)", "lcc_new", BLUE),
    ("Cosine, un-centred", "cosine_uncentred", LBLUE),
    ("Bin + 45°", "bin_rot45", ORANGE),
    ("Bin (Bodis 2009)", "bin_Bodis09", ORANGE),
    ("Quad-tree (Castillo 2013)", "tree_Castillo13", GREY),
    ("NN (Pierens 2012)", "nn_Pierens12", GREY),
]

# (b) mean-centring ablation: (regime, cosine sep, lcc sep)
ABLATION = [("dense $^1$H–$^{15}$N", 0.59, 0.7549),
            ("sparse $^1$H–$^{13}$C", 0.7444516781, 0.7446607271)]


def main() -> None:
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(9.6, 4.0), gridspec_kw=dict(width_ratios=[1.35, 1.0]))

    # ---- panel (a): point + 95% CI (horizontal) ----
    ys = np.arange(len(ROWS))[::-1]  # top row first
    for y, (label, key, col) in zip(ys, ROWS):
        sep = R[key]["sep"]
        lo, hi = R[key].get("cluster_ci95", R[key]["ci95"])  # compound-level CI (respects clustering)
        axa.plot([lo, hi], [y, y], color=col, lw=2.4, solid_capstyle="round", zorder=2)
        axa.plot([sep], [y], "o", color=col, ms=8, zorder=3,
                 markeredgecolor="white", markeredgewidth=1.2)
        axa.text(hi + 0.02, y, f"{sep:.2f}", va="center", ha="left", fontsize=8, color="#333")
    axa.set_yticks(ys)
    axa.set_yticklabels([r[0] for r in ROWS], fontsize=8.5)
    axa.set_xlim(0, 1.08)
    axa.set_xlabel("sparse $^1$H–$^{13}$C separation  (95% CI)", fontsize=9)
    axa.set_title("(a) Statistical strength", fontsize=10, loc="left")
    axa.axvline(0, color="0.7", lw=0.6)
    axa.grid(axis="x", alpha=0.25)
    # ---- panel (b): mean-centring ablation, grouped bars ----
    x = np.arange(len(ABLATION))
    w = 0.36
    cos = [a[1] for a in ABLATION]
    lcc = [a[2] for a in ABLATION]
    axb.bar(x - w / 2, cos, w, label="un-centred cosine", color=GREY)
    axb.bar(x + w / 2, lcc, w, label="STCC (mean-centred)", color=BLUE)
    for xi, (c, l) in zip(x, zip(cos, lcc)):
        axb.text(xi - w / 2, c + 0.015, f"{c:.2f}", ha="center", va="bottom", fontsize=8)
        axb.text(xi + w / 2, l + 0.015, f"{l:.2f}", ha="center", va="bottom", fontsize=8)
    # annotate the dense gap only (double-arrow between bar tops; label placed clear, above)
    gap = lcc[0] - cos[0]
    axb.annotate("", xy=(0 + w / 2, lcc[0]), xytext=(0 - w / 2, cos[0]),
                 arrowprops=dict(arrowstyle="<->", color="#b03030", lw=1.1))
    axb.annotate(f"+{gap:.2f} from\nmean-centring", xy=(0, lcc[0]), xytext=(0.02, 0.90),
                 fontsize=7.4, color="#b03030", ha="center", va="center",
                 arrowprops=dict(arrowstyle="->", color="#b03030", lw=0.8))
    axb.text(1, cos[1] + 0.05, "no gap", fontsize=7.6, color="#666", ha="center")
    axb.set_xticks(x)
    axb.set_xticklabels([a[0] for a in ABLATION], fontsize=8.5)
    axb.set_ylim(0, 1.0)
    axb.set_ylabel("separation", fontsize=9)
    axb.set_title("(b) Mean-centring ablation", fontsize=10, loc="left")
    axb.legend(fontsize=7.6, loc="upper center", frameon=False)
    axb.grid(axis="y", alpha=0.25)

    fig.tight_layout()
    fig.savefig(HERE / "fig3_stats_ablation.png", dpi=200, bbox_inches="tight")
    print("wrote fig3_stats_ablation.png")


if __name__ == "__main__":
    main()
