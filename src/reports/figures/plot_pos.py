import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

site_ends = [("d", "donor"), ("a", "acceptor")]
tr_types = [("protein_coding", "protein-coding"), ("lncRNA", "lncRNA")]
#all_datasets = ["MANE", "CHESS 3", "GENCODE", "RefSeq"]
all_datasets = ["MANE"]

shift = {"d" : list(range(25, 37)), "a" : list(range(25, 37))}
#shift = {"d" : list(range(20, 42)), "a" : list(range(20, 42))}
data = result = pd.read_csv("../../../data/processed/splice_sites.csv")

fig_idx = 0
fig_letter = "abcd"
nn = 406

def tick(s, type):
	s -= 30
	if s == 0:
		return "0 ({})".format("G" if type == "d" else "A")
	elif s == 1:
		return "+1 ({})".format("T" if type == "d" else "G")
	return "{0:+}".format(s)

for tr_type_idx, (now_tr_type, tr_type_description) in enumerate(tr_types):
	for site_end, site_end_description in site_ends:
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
				now_set = now_set + " \\ MANE"

			print(dataset, site_end, now_tr_type, inMane)
			print(now_df.head())

			legend = []
			for now_shift in shift[site_end]:
				rel_shift = now_shift - 30
				legend.append(tick(now_shift, site_end))
				y = [0 for _ in range(0, nn)]
				x = list(range(0, nn))

				for _, row in now_df.iterrows():
					y[row["cons_" + str(rel_shift)]] += 1

				total_y = sum(y)
				y = [float(yi) / total_y for yi in y]

				if rel_shift == 0 or rel_shift == 1:
					now_style = "solid"
				elif rel_shift > 1:
					now_style = "dashed"
				else:
					now_style = "dotted"
				p = plt.plot(x, y, linestyle=now_style, linewidth=.5)
			plt.legend(legend, title="Position shift", ncol=3, loc=9)

			if fig_idx == 0:
				print(legend)
			plt.xlabel("Genomes")
			plt.ylabel("Normalized number of splice sites")
			plt.title("Conservation of certain positions around " + site_end_description + "\n splice sites in " + tr_type_description + " transcripts, " + now_set)
			plt.savefig("out/cons_pos_" + fig_letter[fig_idx] + ".pdf", bbox_inches="tight")
			plt.cla()
			fig_idx += 1


