#!/usr/bin/env python3.11
"""Regenerate the benchmark figures from the JSON the harnesses write.

- comparison_13c.png : separation + margin per method on the sparse 1H-13C regime
                       (read live from comparison_13c.json).
- comparison_all.png : separation across BOTH regimes for every method.

The dense 1H-15N numbers come from the published protein run (method_comparison.json is
not regenerable here without the private Bruker data); the un-centred-cosine protein
separation (0.59) is the ablation value documented in the SI. The sparse 1H-13C numbers
are read live so the figures always match the current comparison_13c.json.

ponytail: no seaborn, just matplotlib; run `python3.11 make_plots.py`.
"""
import json
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).parent

# Published dense-protein 1H-15N separation/margin (from method_comparison.json; the cosine
# row is the SI ablation value, un-centred cosine on the same protein images -> sep 0.59).
PROTEIN = {
    "bin_Bodis09":      (0.2924, 0.2371),
    "bin_rot45":        (0.2497, 0.1973),
    "tree_Castillo13":  (0.0283, -0.0102),
    "nn_Pierens12":     (0.0360, 0.0298),
    "cosine_uncentred": (0.59, None),   # SI ablation: un-centred cosine, protein
    "lcc_new":          (0.7549, 0.7096),
}

LABELS = {
    "bin_Bodis09": "Bin (Bodis 2009)",
    "bin_rot45": "Bin +45°",
    "tree_Castillo13": "Tree (Castillo 2013)",
    "nn_Pierens12": "NN (Pierens 2012)",
    "cosine_uncentred": "Cosine (un-centred)",
    "lcc_new": "LCC (this work)",
}
ORDER = ["bin_Bodis09", "bin_rot45", "tree_Castillo13", "nn_Pierens12", "cosine_uncentred", "lcc_new"]


def main() -> None:
    j13 = json.load(open(HERE / "comparison_13c.json"))["rows"]

    # --- update comparison_all.json (single source of truth for the master figure) ---
    allj = {}
    for m in ORDER:
        s15, m15 = PROTEIN[m]
        r = j13.get(m, {})
        allj[m] = {"sep15": s15, "mar15": m15,
                   "sep13": r.get("separation"), "mar13": r.get("margin")}
    json.dump(allj, open(HERE / "comparison_all.json", "w"), indent=2)

    # --- Figure: sparse 1H-13C separation + margin ---
    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    x = np.arange(len(ORDER))
    seps = [j13[m]["separation"] for m in ORDER]
    mars = [j13[m]["margin"] for m in ORDER]
    ax.bar(x - 0.2, seps, 0.38, label="separation", color="#2c6fbb")
    ax.bar(x + 0.2, mars, 0.38, label="margin", color="#e08a1e")
    for xi, (s, mg) in enumerate(zip(seps, mars)):
        ax.text(xi - 0.2, s + 0.01, f"{s:.2f}", ha="center", va="bottom", fontsize=8)
        ax.text(xi + 0.2, mg + 0.01, f"{mg:.2f}", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[m] for m in ORDER], rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("score")
    ax.set_title("Sparse small-molecule $^1$H–$^{13}$C (5 same-compound pairs, 40 different)")
    ax.set_ylim(-0.05, 1.0)
    ax.axhline(0, color="0.6", lw=0.6)
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(HERE / "comparison_13c.png", dpi=200)
    plt.close(fig)

    # --- Figure: separation across both regimes (paper Figure 2) ---
    fig, ax = plt.subplots(figsize=(8.6, 4.4))
    x = np.arange(len(ORDER))
    s15 = [allj[m]["sep15"] for m in ORDER]
    s13 = [allj[m]["sep13"] for m in ORDER]
    ax.bar(x - 0.2, s15, 0.38, label="dense protein $^1$H–$^{15}$N", color="#2c6fbb")
    ax.bar(x + 0.2, s13, 0.38, label="sparse small-molecule $^1$H–$^{13}$C", color="#e08a1e")
    for xi, (a, b) in enumerate(zip(s15, s13)):
        ax.text(xi - 0.2, a + 0.01, f"{a:.2f}", ha="center", va="bottom", fontsize=8)
        ax.text(xi + 0.2, b + 0.01, f"{b:.2f}", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[m] for m in ORDER], rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("separation (mean same − mean different)")
    ax.set_title("Same/different separation across both regimes")
    ax.set_ylim(0, 0.95)
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(HERE / "comparison_all.png", dpi=200)
    plt.close(fig)
    print("wrote comparison_13c.png, comparison_all.png, comparison_all.json")


if __name__ == "__main__":
    main()
