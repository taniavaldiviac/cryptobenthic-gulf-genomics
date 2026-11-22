#!/usr/bin/env python3
"""Figure: decomposition of the eDNA deficit by in silico filter (100% stacked bars).

Input: the counts CSV written by decompose_deficit.py.
Usage: python plot_deficit_decomposition.py counts.csv --threshold 3 --out fig_prefix
"""
import argparse, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CATS = [("0_fuera_de_lista", "Outside checklist", "#c9c7c0", 1.0),
        ("1_sin_referencia", "No teleo reference", "#2f6fd6", 0.9),
        ("2_cebador", "Primer failure", "#c9c7c0", 0.6),
        ("3_no_resoluble", "Not resolvable", "#2f6fd6", 0.45),
        ("4_sin_explicacion", "Unexplained in silico", "#2aa876", 0.8)]
LABEL = {"K": "Cryptobenthic", "C": "Conspicuous"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("counts")
    ap.add_argument("--threshold", type=float, default=3)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rows = [r for r in csv.DictReader(open(a.counts)) if float(r["umbral_dmin"]) == a.threshold]
    rows.sort(key=lambda r: ["K", "C"].index(r["grupo"]))
    fig, ax = plt.subplots(figsize=(8, 2.8))
    for i, r in enumerate(rows):
        tot = int(r["n_no_detectadas"])
        left = 0
        y = len(rows) - 1 - i
        for key, name, col, alpha in CATS:
            n = int(r[key])
            w = n / tot * 100
            ax.barh(y, w, left=left, color=col, alpha=alpha, edgecolor="white", height=0.55,
                    label=name if i == 0 else None)
            if w >= 9:
                ax.text(left + w / 2, y, f"{n} ({round(w)}%)", ha="center", va="center", fontsize=9)
            elif w >= 3:
                ax.text(left + w / 2, y, f"{n}", ha="center", va="center", fontsize=9)
            left += w
        ax.text(0, y + 0.38, f"{LABEL.get(r['grupo'], r['grupo'])} ({tot} species missed)",
                fontsize=10, va="bottom")
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.5, len(rows) - 0.2)
    ax.set_yticks([])
    ax.set_xlabel("Species missed by eDNA (%)")
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.legend(ncol=5, fontsize=8, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.32))
    fig.tight_layout()
    for ext in ("png", "svg", "pdf"):
        fig.savefig(f"{a.out}.{ext}", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
