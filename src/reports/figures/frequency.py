import os
import math
import string
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


minus = "*"
labels = ["Less-supported", "Well-supported"]
abc = string.ascii_lowercase

data = pd.read_csv("../../../data/processed/model_out.csv")
data = data[(data["dataset"] != "MANE") & (data["dataset"] != "Random") & (data["inMANE"] == 0)]
data["dataset_title"] = data.apply(lambda x: x["dataset"] + ("*" if x["gene_type"] == "protein_coding" else ""), axis=1)
data["conservation_legend"] = data.apply(lambda x: labels[x["conserved"]], axis=1)

site_ends = [("d", "donor"), ("a", "acceptor")]
tr_types = [("protein_coding", "protein-coding"), ("lncRNA", "lncRNA")]

fig_idx = 0
for tr_type_idx, (now_tr_type, tr_type_description) in enumerate(tr_types):
	for (now_end, end_title) in site_ends:
		now_df = data[(data["gene_type"] == now_tr_type) & (data["site_type"] == now_end)]
		now_df = now_df.melt(id_vars=["dataset_title", "conservation_legend"], value_vars=["snp_0", "snp_1"])
		now_df = now_df[now_df["value"] > 0]
		count = dict()
		for _, row in now_df.iterrows():
			now_title = row["dataset_title"]
			if not now_title in count:
				count[now_title] = [0, 0]
			count[now_title][labels.index(row["conservation_legend"])] += 1

		title_n = dict()
		for title, n in count.items():
			title_n[title] = title + "\nN=" + str(n[0]) + "/" + str(n[1])

		if now_tr_type == "protein_coding":
			plt.ylim(0, 0.21)
		else:
			plt.ylim(0, 0.14)

		print(title_n)
		print(now_df)
		now_df["dataset_title_n"] = now_df.apply(lambda x: title_n[x["dataset_title"]], axis=1)
		b = sns.boxplot(data=now_df, x="dataset_title_n", y="value", hue="conservation_legend", hue_order=labels, showfliers=False)
		if fig_idx != 0:
			b.legend().remove()
		else:
			b.legend().set_title("Splice site support level")
#			sns.move_legend(b, "upper left")
		plt.ylabel("SNP Frequency")
		plt.xlabel("Dataset")
		plt.title("Frequency of SNPs overlapping canonical dinucleotides of " + end_title + " \n splice sites from " + tr_type_description + " genes of different splice site support status")
		plt.savefig("appendix/figure_frequency_" + abc[fig_idx] + ".pdf", bbox_inches="tight")
		fig_idx += 1
		plt.cla()



