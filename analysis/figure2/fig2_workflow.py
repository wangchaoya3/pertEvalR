# Figure 2: workflow schematic of pertEvalR (no data needed).
# Requires Python 3 with matplotlib. Run from this folder:
#   python fig2_workflow.py
# Writes fig2_workflow.png (300 dpi) and fig2_workflow.pdf next to this script.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ink, ink2, grid, surface = "#0b0b0b", "#52514e", "#c9c7c1", "#fcfcfb"
blue, orange, green, yellow = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
plt.rcParams["font.family"] = "DejaVu Sans"

fig, ax = plt.subplots(figsize=(11, 6), dpi=300)
fig.patch.set_facecolor(surface); ax.set_facecolor(surface)
ax.set_xlim(0, 100); ax.set_ylim(0, 56); ax.axis("off")

def box(x, y, w, h, num, title, func, lines, col):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.25,rounding_size=1.2",
                                fc="white", ec=col, lw=1.6))
    ax.add_patch(FancyBboxPatch((x, y + h - 4.2), w, 4.2, boxstyle="round,pad=0.25,rounding_size=1.2",
                                fc=col, ec=col, lw=1.6))
    ax.text(x + 1.2, y + h - 2.1, f"{num}  {title}" if num else title, color="white",
            fontsize=10.5, fontweight="bold", va="center", ha="left")
    ax.text(x + 1.2, y + h - 6.3, func, color=ink, fontsize=7.8, family="DejaVu Sans Mono",
            va="center", ha="left")
    for i, t in enumerate(lines):
        ax.text(x + 1.2, y + h - 9.6 - i * 2.7, t, color=ink2, fontsize=7.3, va="center", ha="left")

def arrow(x1, y1, x2, y2, style="-|>", ls="-", col=ink2, rad=0.0):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style, mutation_scale=11,
                                 lw=1.3, color=col, linestyle=ls,
                                 connectionstyle=f"arc3,rad={rad}"))

# Row 1: data -> diagnose -> split -> predict -> score
W, H, Y1 = 18, 19, 33
xs = [0.5, 20.5, 40.5, 60.5, 80.5]
box(xs[0], Y1, W, H, "", "Data", "pe_data()", ["matrix (dense or sparse)", "SingleCellExperiment", "Seurat v5", "-> pseudobulk profiles"], ink2)
box(xs[1], Y1, W, H, "1", "Diagnose", "pe_signal()", ["systematic fraction", "effect size", "effective genes", "split-half reliability"], blue)
box(xs[2], Y1, W, H, "2", "Split", "pe_split()", ["unseen singles", "unseen combinations", "no leakage through", "combinations"], blue)
box(xs[3], Y1, W, H, "3", "Predict", "baseline_*()", ["no change, training mean,", "additive, linear;", "pe_prediction() wraps", "any external model"], orange)
box(xs[4], Y1, W, H, "4", "Score", "pe_evaluate()", ["6 metrics, each per", "held-out perturbation:", "MSE (2), PearsonDelta (3),", "discrimination (PDS)"], orange)
for a, b in zip(xs[:-1], xs[1:]):
    arrow(a + W + 0.4, Y1 + H / 2, b - 0.4, Y1 + H / 2)

# Row 2: stratify (left, data-only), compare, calibrate
Y2 = 6
box(0.5, Y2, 28, H, "7", "Stratify", "pe_subset()", ["restrict results and bounds", "to reliable perturbations;", "threshold fixed from the data", "before models are compared"], yellow)
box(35.5, Y2, 30, H, "5", "Compare", "pe_compare()", ["paired differences per perturbation", "bootstrap 95% CI, signed-rank test,", "Holm correction,", "direction-aware verdict"], green)
box(69.5, Y2, 30, H, "6", "Calibrate", "pe_bounds(), pe_calibrate()", ["0 = best uninformative baseline", "1 = split-half technical replicate", "pe_meta_metrics(): DRF and BDS", "per metric and dataset"], green)
sx, sy = xs[4] + W / 2, Y1 - 0.4
arrow(sx, sy, 35.5 + 15, Y2 + H + 0.6)
arrow(sx, sy, 69.5 + 15, Y2 + H + 0.6)
# reliability feeds stratification
arrow(xs[1] + W / 2, Y1 - 0.4, 14.5, Y2 + H + 0.6, ls=(0, (3, 2)), col=blue)
ax.text(23.2, 28.4, "reliability", color=blue, fontsize=8, style="italic", ha="left")
# stratification restricts compare and calibrate
arrow(28.9, Y2 + 8, 35.1, Y2 + 8, ls=(0, (3, 2)), col=yellow)
ax.text(32, Y2 + 9.6, "restricts", color="#a87400", fontsize=7, style="italic", ha="center")

ax.text(0.5, 2.2, "Report: dataset signal  |  calibrated next to raw scores  |  paired differences with CIs over splits  |  results on reliable perturbations",
        fontsize=8.4, color=ink, ha="left", va="center")
fig.savefig("fig2_workflow.png", dpi=300, facecolor=surface, bbox_inches="tight")
fig.savefig("fig2_workflow.pdf", facecolor=surface, bbox_inches="tight")
