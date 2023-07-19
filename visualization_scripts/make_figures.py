# This script makes various figures.

import os 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('agg')

working_dir = "/home/vpc/datasets/plasma_chipseq/pipeline"
PATH_ct = os.path.join(working_dir, "analysis/ct_fractions.tsv")

ct_main = pd.read_csv(PATH_ct, "\t")

ct = ct_main[["Patient", "Timepoint", "Sample", "Mutation_ctDNA_fraction", "CN_ct"]]
ct["Mutation_ctDNA_fraction"], ct["CN_ct"] = ct["Mutation_ctDNA_fraction"] * 100, ct["CN_ct"] * 100 #convert to percentage
ct.columns = ["Patient", "Timepoint", "Sample", "Mut", "CN"]
ct["Mut-CN"] = ct["Mut"] - ct["CN"] # the difference between mut and CN based ct estimates 

# =================================================================================
# Bar plot to compare mutation based and CN based ct fraction estimates
fig, ax = plt.subplots(figsize=(8, 6))
ax.bar(ct["Sample"], ct["Mut-CN"], edgecolor = "black", color="whitesmoke", zorder = 10)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_xticklabels(ct["Sample"], rotation=90)
ax.set_ylabel("Mutation based estimate - \nCN based estimate (%)")
ax.set_title('ctDNA difference')

for tick in ax.get_yticks():
    ax.axhline(y=tick, linestyle='dashed', color='black', lw = 0.5)

fig.tight_layout()
fig.savefig(os.path.join(working_dir, "figures/ct_fraction_difference.png"))
# =================================================================================

# =================================================================================
# Scatter plot showing ct % from mutation and CN estimates 
timepoint_colors = {"Baseline": "#17becf", "OT1": "lightgray", "OT2": "darkgray", "OT3": "dimgray"}
ct["Colors"] = ct["Timepoint"].map(timepoint_colors)

fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(ct["CN"], ct["Mut"], color = ct["Colors"], zorder = 10)

ax.set_xticks(range(0, 101, 10))
ax.set_yticks(range(0, 101, 10))

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_xlabel("Mutation based estimate (%)")
ax.set_ylabel("CN based estimate (%)")
ax.set_title("ctDNA estimates (%)")

# Linear regression
coefficients = np.polyfit(ct["CN"], ct["Mut"], deg=1)
poly = np.poly1d(coefficients)
x = np.linspace(0, 100, 100)
y = poly(x)

# Plot the best-fit line
ax.plot(x, y, color='red', linestyle='dashed', lw = 0.5)

# Calculate R-squared value
residuals = ct["Mut"] - poly(ct["CN"])
ss_residuals = np.sum(residuals**2)
ss_total = np.sum((ct["Mut"] - np.mean(ct["Mut"]))**2)
r_squared = 1 - (ss_residuals / ss_total)

# Add R-squared value and equation to the plot
equation = f'y = {coefficients[0]:.2f}x + {coefficients[1]:.2f}'
r_squared_text = f'R-squared = {r_squared:.2f}'
ax.text(4, 70, equation, color='red')
ax.text(4, 66, r_squared_text, color='red')
ax.text(100, 80, "Best fit line", color='red')

# legend for the points 
legend_elements = [plt.Line2D([0], [0], marker='o', color='w', label=key, markerfacecolor=value, markersize=8) for key, value in timepoint_colors.items()]
ax.legend(handles=legend_elements)

fig.tight_layout()
fig.savefig(os.path.join(working_dir, "figures/ct_fraction.png"))
# =================================================================================
