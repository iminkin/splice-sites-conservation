import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

abc = "abcdef"

cons_caption = ["Less-supported (\u2265 1 less-supported site)", "Well-supported (all sites well-supported/from MANE)"]

fig_idx = 0
df = pd.read_csv("../../../data/processed/transcripts.csv")
df = df[(df["total_sites"] > 0) & (df["inMANE"] == 0)]
df["introns_count"] = df.apply(lambda x: int(min(x["total_sites"] / 2, 13)), axis=1)
df["Transcript support"] = df.apply(lambda x: cons_caption[x["conservation_status"]], axis=1)
datasets = ["GENCODE", "RefSeq", "CHESS 3"]
types = [("protein_coding", "protein-coding"), ("lncRNA", "lncRNA")]
print(df)
for type_idx, (type, type_title) in enumerate(types):
	for dataset in datasets:
		data = df[(df["gene_type"] == type) & (df["dataset"] == dataset)]
		ax = sns.countplot(data, x="introns_count", hue="Transcript support", hue_order=cons_caption)
		if fig_idx % 3 != 0:
			ax.legend_.remove()
		label = [x.get_text() for x in ax.get_xticklabels()]
		label[-1] = "> " + label[-2]
		ax.set_xticklabels(label)

		plt.xlabel("Number of introns")
		plt.ylabel("Number of transcripts")
		now_dataset = dataset + "*" if type == "protein_coding" else dataset
		plt.title("Number of less/well-supported " + type_title + "\n transcripts in " + now_dataset + " dataset")

		plt.savefig("appendix/figure_transcripts_" + abc[fig_idx] + ".pdf")
		fig_idx += 1
		plt.cla()
