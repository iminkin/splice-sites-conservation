import os
import math
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import PathPatch

def adjust_box_widths(g, fac):
	for ax in g.axes:
		for c in ax.get_children():
			if isinstance(c, PathPatch):
				p = c.get_path()
				verts = p.vertices
				verts_sub = verts[:-1]
				xmin = np.min(verts_sub[:, 0])
				xmax = np.max(verts_sub[:, 0])
				xmid = 0.5*(xmin+xmax)
				xhalf = 0.5*(xmax - xmin)
				xmin_new = xmid-fac*xhalf
				xmax_new = xmid+fac*xhalf
				verts_sub[verts_sub[:, 0] == xmin, 0] = xmin_new
				verts_sub[verts_sub[:, 0] == xmax, 0] = xmax_new
				for l in ax.lines:
					if np.all(l.get_xdata() == [xmin, xmax]):
						l.set_xdata([xmin_new, xmax_new])

site_ends = [("d", "donor"), ("a", "acceptor")]
tr_types = [("protein_coding", "protein-coding"), ("lncRNA", "lncRNA")]
dataset = ["GENCODE", "RefSeq", "CHESS 3"]

data = pd.read_csv("../../../data/processed/model_out.csv")
data = data[(data["dataset"] != "Random") & (data["inMANE"] == 0)]

fig_idx = 0
fig_letter = "abcd"

data = data.assign(tick=0)
for i, row in data.iterrows():
	idx = dataset.index(row["dataset"])
	if row["gene_type"] == "lncRNA":
		idx += 3
	data.loc[i, "tick"] = idx

for end, end_title in site_ends:
	bx = sns.boxplot(x="tick", y="reuse", hue="conserved", data=data[data["site_type"] == end], medianprops={'color': 'red', 'ls': 'dashed', 'lw': 1}, showfliers=False, whis=(0, 95))
#	bx = sns.violinplot(x="tick", y="reuse", hue="conserved", data=data[data["site_type"] == end])
	plt.title("Sharing of less/well-supported " + end_title + " splice sites \n between different transcripts of the same gene")
	plt.ylabel("Number of transcripts sharing the site")
	plt.xlabel("Dataset")
	coding = [l + "*,\ncoding" for l in dataset]
	rna = [l + ",\nlncRNA" for l in dataset]
	plt.ylim(0, 20)
	plt.yticks(list(range(0, 20, 2)))
	print(range(0, 20, 2))
	plt.xticks(list(range(0, 6)), coding + rna)
	plt.axvline(x=2.5, color="black", linestyle="--", linewidth=.5)
	plt.subplots_adjust(bottom=0.15)

	bx.legend_.set_title("Level of support")
	bx.legend_.texts[0].set_text("Less-supported")
	bx.legend_.texts[1].set_text("Well-supported")

	fig = plt.gcf()
	adjust_box_widths(fig, 0.9)
	if fig_idx == 1:
		bx.legend().remove()

	plt.savefig("out/figure_usage_" + fig_letter[fig_idx] + ".pdf", bbox_inches="tight")
	fig_idx += 1
	plt.cla()




