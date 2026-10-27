# Workflow schematic of pertEvalR (Figure 1 of the PeerJ manuscript); no data.
# Requires Python 3 with matplotlib. Run from this folder:
#   python fig2_workflow.py
# Drawn at final print size (7.5 x 4.7 inches) so that the font sizes below
# are the sizes in the published figure; saved at 400 dpi (3000 px wide).
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ink, ink2, surface = "#0b0b0b", "#3d3c3a", "#ffffff"
blue, orange, green, yellow = "#2a78d6", "#eb6834", "#1baf7a", "#d99400"
plt.rcParams["font.family"] = "DejaVu Sans"

FIG_W, FIG_H = 7.5, 4.7
fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=400)
ax = fig.add_axes([0, 0, 1, 1])
fig.patch.set_facecolor(surface)
ax.set_xlim(0, FIG_W)
ax.set_ylim(0, FIG_H)
ax.axis("off")

TITLE_PT, FUNC_PT, BODY_PT = 9.5, 8, 8
HEAD = 0.24        # height of the coloured title band (inches)
LINE = 0.155       # line spacing of body text (inches)


def box(x, y, w, h, num, title, func, lines, col):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.06",
                                fc="white", ec=col, lw=1.3, zorder=2))
    ax.add_patch(FancyBboxPatch((x, y + h - HEAD), w, HEAD,
                                boxstyle="round,pad=0,rounding_size=0.06",
                                fc=col, ec=col, lw=1.3, zorder=3))
    label = f"{num}  {title}" if num else title
    ax.text(x + 0.07, y + h - HEAD / 2, label, color="white", fontsize=TITLE_PT,
            fontweight="bold", va="center", ha="left", zorder=4)
    ax.text(x + 0.07, y + h - HEAD - 0.15, func, color=ink, fontsize=FUNC_PT,
            family="DejaVu Sans Mono", va="center", ha="left", zorder=4)
    for i, t in enumerate(lines):
        ax.text(x + 0.07, y + h - HEAD - 0.36 - i * LINE, t, color=ink2,
                fontsize=BODY_PT, va="center", ha="left", zorder=4)


def arrow(p1, p2, col=ink2, ls="-"):
    # Dashed shafts: plain dashed line plus a short solid arrowhead, so the
    # dash pattern never cuts into the head.
    if ls != "-":
        (x1, y1), (x2, y2) = p1, p2
        d = ((x2 - x1) ** 2 + (y2 - y1) ** 2) ** 0.5
        ux, uy = (x2 - x1) / d, (y2 - y1) / d
        head = 0.09
        mx, my = x2 - ux * head, y2 - uy * head
        ax.plot([x1, mx], [y1, my], color=col, lw=1.2, linestyle=ls,
                dash_capstyle="butt", zorder=5)
        p1 = (mx - ux * 0.005, my - uy * 0.005)
    ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle="-|>", mutation_scale=9,
                                 lw=1.2, color=col, linestyle="-",
                                 joinstyle="miter", capstyle="butt",
                                 shrinkA=0, shrinkB=0, zorder=5))


def line(xs, ys, col=ink2):
    ax.plot(xs, ys, color=col, lw=1.2, solid_capstyle="butt", zorder=5)


# Row 1: Data -> 1 Diagnose -> 2 Split -> 3 Predict -> 4 Score ---------------
W1, H1, Y1 = 1.34, 1.24, 3.34
gap = (FIG_W - 0.1 - 5 * W1) / 4
xs = [0.05 + i * (W1 + gap) for i in range(5)]
box(xs[0], Y1, W1, H1, "", "Data", "pe_data()",
    ["matrix (dense or", "sparse), Single-", "CellExperiment,", "Seurat v5"], ink2)
box(xs[1], Y1, W1, H1, "1", "Diagnose", "pe_signal()",
    ["systematic fraction", "effect size", "effective genes", "split-half reliability"], blue)
box(xs[2], Y1, W1, H1, "2", "Split", "pe_split()",
    ["unseen singles", "unseen pairs", "no leakage through", "combinations"], blue)
box(xs[3], Y1, W1, H1, "3", "Predict", "baseline_*()",
    ["no change, mean,", "additive, linear;", "pe_prediction():", "external models"], orange)
box(xs[4], Y1, W1, H1, "4", "Score", "pe_evaluate()",
    ["6 metrics for each", "held-out", "perturbation", "(MSE, PearsonΔ, PDS)"], orange)
ymid = Y1 + (H1 - HEAD) / 2
for a, b in zip(xs[:-1], xs[1:]):
    arrow((a + W1 + 0.02, ymid), (b - 0.02, ymid))

# Row 2: 7 Stratify, 5 Compare, 6 Calibrate ------------------------------------
H2, Y2 = 1.30, 0.62
W7, W5, W6 = 2.2, 2.0, 2.35
x7 = 0.05
x6 = FIG_W - 0.05 - W6
x5 = x6 - 0.3 - W5
box(x7, Y2, W7, H2, "7", "Stratify", "pe_subset()",
    ["restrict results and bounds", "to reliable perturbations;", "threshold fixed from the",
     "data before comparison"], yellow)
box(x5, Y2, W5, H2, "5", "Compare", "pe_compare()",
    ["paired differences per", "perturbation, bootstrap", "95% CI, signed-rank test,",
     "Holm correction"], green)
box(x6, Y2, W6, H2, "6", "Calibrate", "pe_bounds(), pe_calibrate()",
    ["0 = best uninformative baseline", "1 = split-half replicate",
     "pe_meta_metrics(): DRF, BDS"], green)

top2 = Y2 + H2
bus = 2.80
sx = xs[4] + W1 / 2
c5, c6 = x5 + W5 / 2, x6 + W6 / 2
line([sx, sx], [Y1, bus])
line([c5, sx], [bus, bus])
arrow((c5, bus), (c5, top2 + 0.01))
arrow((c6, bus), (c6, top2 + 0.01))
ax.text(c5 + 0.06, bus + 0.05, "per-perturbation scores", color=ink2, fontsize=7.5,
        style="italic", ha="left", va="bottom")

dx = xs[1] + W1 / 2
arrow((dx, Y1), (dx, top2 + 0.01), col=blue, ls=(0, (4, 2)))
ax.text(dx + 0.06, (Y1 + top2) / 2, "reliability", color=blue, fontsize=8,
        style="italic", ha="left", va="center")

yr = Y2 + 0.45
arrow((x7 + W7 + 0.02, yr), (x5 - 0.02, yr), col=yellow, ls=(0, (4, 2)))
ax.text((x7 + W7 + x5) / 2, yr + 0.06, "restricts\n5 and 6", color="#9a6a00",
        fontsize=7, style="italic", ha="center", va="bottom")

ax.text(0.05, 0.28, "Report: dataset signal | calibrated and raw scores | "
        "paired differences with CIs | results on reliable perturbations",
        fontsize=8, color=ink, ha="left", va="center")

fig.savefig("fig2_workflow.png", dpi=400, facecolor=surface)
fig.savefig("fig2_workflow.pdf", facecolor=surface)
