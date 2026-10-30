# Figure 2: mean chordal distance by zero proportion, loadings (top) and scores
# (bottom), one column per transformation.
#   python3 simulations/replication-figure.py   -> output/figure2.{pdf,png}

import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

infile = sys.argv[1] if len(sys.argv) > 1 else "output/replication.csv"
outfile = sys.argv[2] if len(sys.argv) > 2 else "output/figure2.pdf"

INK, MUTED, RULE = "#52514e", "#8a8984", "#e4e3df"
plt.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"], "font.size": 8,
    "axes.edgecolor": INK, "axes.labelcolor": INK, "axes.linewidth": 0.6,
    "xtick.color": INK, "ytick.color": INK, "text.color": INK,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
})

methods = [  # (column value, label, colour, line style, marker, width)
    ("Truncated copula PCA", "TLGC-PCA", "#2a78d6", "-", "o", 1.6),
    ("Standard PCA (unscaled)", "PCA (unscaled)", "#eb6834", "-", "s", 1.1),
    ("Standard PCA (scaled)", "PCA (scaled)", "#1baf7a", "--", "^", 1.1),
    ("PoissonPCA", "PoissonPCA", "#eda100", ":", "D", 1.1),
    ("ZIPFA", "ZIPFA", "#e87ba4", "-.", "v", 1.1),
]
transfs = [("scaling", "Scaling"), ("cubic", "Cubic"),
           ("exp", "Exponential"), ("mixed", "Mixed")]
rows = [("dist_V", "Loadings\nchordal distance", 1.79),
        ("dist_U", "Scores\nchordal distance", 1.98)]

d = pd.read_csv(infile)
d["zeros"] = 100 * d["trunc_prop"]
means = d.groupby(["transf", "method", "zeros"])[["dist_V", "dist_U"]].mean()

fig, axes = plt.subplots(2, 4, figsize=(6.3, 3.9), sharex=True, sharey=True)
for i, (col, ylab, rand) in enumerate(rows):
    for j, (tr, title) in enumerate(transfs):
        ax = axes[i, j]
        ax.axhline(rand, color=MUTED, lw=0.7, ls=(0, (4, 3)), zorder=1)
        for key, lab, colour, ls, mk, lw in methods:
            m = means.loc[(tr, key)][col]
            ax.plot(m.index, m.values, color=colour, ls=ls, lw=lw, marker=mk,
                    ms=3.5, mec="white", mew=0.4, label=lab,
                    zorder=3 if key.startswith("Trunc") else 2)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_ylim(0, 2.1)
        ax.set_xlim(-3, 93)
        ax.set_xticks([0, 20, 40, 60, 80])
        if i == 0:
            ax.set_title(title, fontsize=10)
        else:
            ax.set_xlabel("Zeros (%)")
        if j == 0:
            ax.set_ylabel(ylab)
axes[0, 3].text(92, 1.79 + 0.05, "random subspace", ha="right", va="bottom",
                fontsize=6.5, color=MUTED)

handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", ncol=5, frameon=False,
           fontsize=8, labelcolor="black", handlelength=2.4, columnspacing=1.6,
           bbox_to_anchor=(0.5, 1.0))
fig.tight_layout(rect=(0, 0, 1, 0.94), h_pad=0.6, w_pad=0.4)
fig.savefig(outfile)
fig.savefig(outfile.rsplit(".", 1)[0] + ".png", dpi=200)
print("Wrote", outfile)
