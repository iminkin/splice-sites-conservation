import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


site_ends = [("d", "donor"), ("a", "acceptor")]
tr_types = [("protein_coding", "protein-coding"), ("lncRNA", "lncRNA")]
all_datasets = ["GENCODE", "RefSeq", "CHESS 3", "MANE"]

data = result = pd.read_csv("../../../data/processed/splice_sites.csv")
data = data[data["dataset"] != "Random"]
data = data[((data["site_type"] == "d") & (data["motif"] != "GT")) | ((data["site_type"] == "a") & (data["motif"] != "AG"))]

fig_idx = 0
fig_letter = "abcd"
nn = 405

for tr_type_idx, (now_tr_type, tr_type_description) in enumerate(tr_types):
	for site_end, site_end_description in site_ends:
		legend = []
		color = []

		data_col = []
		cons_col = []

		rcnt = 0
		ry = [0 for _ in range(0, nn)]
		rx = list(range(0, nn))

		for dataset_idx, dataset in enumerate(all_datasets):
			tr_type = now_tr_type
			if dataset == "MANE" and tr_type == "lncRNA":
				continue
			if dataset == "Random" and tr_type == "lncRNA":
				tr_type = "protein_coding"

			now_set = dataset
			inMane = now_set == "MANE"

			now_df = data[(data["dataset"] == dataset) & (data["site_type"] == site_end) & (data["gene_type"] == tr_type) & (data["inMANE"] == inMane)]
			if (now_set != "MANE" and now_set != "Random") and tr_type == "protein_coding":
				now_set = now_set + "*"

			print(dataset, site_end, now_tr_type, inMane)
			print(now_df.head())

			legend.append(now_set)
			for _, row in now_df.iterrows():
				data_col.append(row["dataset"])
				cons_col.append(row["cons_GTAG"])

#		color.append(p[-1].get_color())
		plt.plot([], [])
		df = pd.DataFrame(zip(data_col, cons_col), columns=["dataset", "cons"])
#		sns.histplot(data=df, x="cons", hue="dataset", multiple="dodge", stat="probability", palette=['#ff7f0e', '#2ca02c', '#d62728', '#9467bd'])
		sns.histplot(data=df, x="cons", hue="dataset", multiple="dodge", palette=['#ff7f0e', '#2ca02c', '#d62728', '#9467bd'])

		print(color)
#		if site_end == "d":
#			plt.legend(legend, title="Dataset", loc=9)
#		print(legend)
		plt.xlabel("Number of genomes")
		plt.ylabel("Number of splice sites in a bin")
		plt.title("Conservation of " + site_end_description + " splice sites in \n" + tr_type_description + " transcripts (non-canonical dinucleotides)")
		plt.savefig("out/figure_gtag_nonc_" + fig_letter[fig_idx] + ".pdf", bbox_inches="tight")
		plt.cla()
		fig_idx += 1


