import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


site_ends = [("d", "donor"), ("a", "acceptor")]
tr_types = [("protein_coding", "protein-coding"), ("lncRNA", "lncRNA")]
all_datasets = ["Random", "GENCODE", "RefSeq", "CHESS 3", "MANE"]

data = result = pd.read_csv("../../../data/processed/splice_sites.csv")

fig_idx = 0
fig_letter = "abcd"
nn = 406

nnmax = 0

for tr_type_idx, (now_tr_type, tr_type_description) in enumerate(tr_types):
	for site_end, site_end_description in site_ends:
		legend = []
		color = dict()

		for dataset_idx, dataset in enumerate(all_datasets):
			tr_type = now_tr_type
			if dataset == "MANE" and tr_type == "lncRNA":
				continue
			if dataset == "Random" and tr_type == "lncRNA":
				tr_type = "protein_coding"

			now_set = dataset
#			inMane = now_set == "MANE"
			inMane = 1 if now_set == "MANE" else 0
			now_df = data[(data["dataset"] == dataset) & (data["site_type"] == site_end) & (data["gene_type"] == tr_type) & (data["inMANE"] == inMane)]
			if (now_set != "MANE" and now_set != "Random") and tr_type == "protein_coding":
				now_set = now_set + "*"

			print(dataset, site_end, now_tr_type, inMane)
			print(now_df.head())

			cnt = 0
			y = [0 for _ in range(0, nn)]
			for _, row in now_df.iterrows():
				nnmax = max(nnmax, row["cons_GTAG"])
				y[row["cons_GTAG"]] += 1
				cnt += 1
			print(cnt)
			total_y = sum(y)
			y = [float(yi) / total_y for yi in y]
#			y = [sum(y[0:i + 1]) for i in range(0, nn)]
			x = list(range(0, nn))
			p = plt.plot(x, y, linewidth=.5)
			color[dataset] = p[-1].get_color()
			legend.append(now_set)

		print(color)
		if site_end == "d":
			plt.legend(legend, title="Dataset", loc=9)
		print(legend)
#		if tr_type == "protein_coding":
#			plt.ylim(0, 0.04)
		plt.xlabel("Number of genomes")
		plt.ylabel("Normalized number of splice sites")
		plt.title("Conservation of " + site_end_description + " splice sites in \n" + tr_type_description + " transcripts")
		plt.savefig("out/figure_gtag_" + fig_letter[fig_idx] + ".pdf", bbox_inches="tight")
		plt.cla()
		fig_idx += 1


print(nnmax)
