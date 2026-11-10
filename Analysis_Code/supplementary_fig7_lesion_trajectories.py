"""Supplementary Figure 7: recorded lesion-positive proportions by tissue, age class and genotype.

Reads Analysis_Code/figure_inputs/Histology_descriptive.csv (five animals per group; an animal is lesion-positive
when any of its four sections shows a lesion) and redraws the 4 x 2 panel layout of the submitted figure.
usage: python3 supplementary_fig7_lesion_trajectories.py ROOT OUT_PNG
"""
import ast
import csv
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root, out = sys.argv[1:3]
rows = list(csv.DictReader(open(f"{root}/Analysis_Code/figure_inputs/Histology_descriptive.csv", encoding="utf-8-sig")))
colours = {"DD": "#0072B2", "GG": "#D55E00", "GD": "#CC79A7"}
tissues = ["Foot muscle", "Digestive gland", "Gill", "Mantle"]
ages = [("1龄", "1-year"), ("2龄", "2-year")]
plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8})
fig, axes = plt.subplots(4, 2, figsize=(7.1, 8.8))
for i, tissue in enumerate(tissues):
    for j, (code, age) in enumerate(ages):
        ax = axes[i, j]
        for g in ("DD", "GG", "GD"):
            r = next(r for r in rows if r["age_class"] == code and r["genotype"] == g and r["tissue"] == tissue)
            t = ast.literal_eval(r["time_points_h"]); v = [float(x) for x in ast.literal_eval(r["recorded_proportions"])]
            assert all(abs(x * 5 - round(x * 5)) < 1e-9 for x in v), (tissue, age, g, v)  # k of 5 animals
            ax.plot(t, v, marker="o", markersize=3.5, linewidth=1.1, color=colours[g], label=g)
        ax.set_title(f"{tissue}, {age}")
        ax.set_xticks([0, 48, 96, 168]); ax.set_ylim(-0.03, 1.05); ax.set_yticks([0, .2, .4, .6, .8, 1.0])
        ax.set_xlabel("Hours post-injection"); ax.set_ylabel("Recorded lesion-positive proportion")
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(frameon=False, ncol=3, loc="best")
fig.tight_layout()
fig.savefig(out, dpi=400, facecolor="white")
print("saved", out)
