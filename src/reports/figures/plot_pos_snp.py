import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

site_ends = [("d", "donor"), ("a", "acceptor")]
tr_types = [("protein_coding", "protein-coding"), ("lncRNA", "lncRNA")]
all_datasets = ["Random", "GENCODE", "RefSeq", "CHESS 3", "MANE"]
conservation_caption = ["less-supported", "well-supported"]
pattern = ["dashed", "solid"]

shift = list(range(24, 38))

data = result = pd.read_csv("../../../data/processed/model_out.csv")

data["phastCons"] = data.apply(lambda x: 1 if min(x["phastCons_0"], x["phastCons_1"]) > 0.5 else 0, axis=1)

fig_idx = 0
fig_letter = "abcd"

def plot(now_df, pattern, color):
	x = []
	y = []
	for now_shift in shift:
		rel_shift = now_shift - 30
		x.append(str(rel_shift))
		with_snp = now_df[now_df["snp_" + str(rel_shift)] > 0]
		fr = float(len(with_snp.index)) / len(now_df.index)
		y.append(fr)
	if color == "":
		p = plt.plot(x, y, linestyle=pattern)
	else:
		p = plt.plot(x, y, color=color, linestyle=pattern)
	plt.ylim(0.001, 0.025)
	return p[-1].get_color()

def tick(s, type):
	s -= 30
	if s == 0:
		return "0\n{}".format("G" if type == "d" else "A")
	elif s == 1:
		return "+1\n{}".format("T" if type == "d" else "G")
	return "{0:+}".format(s)

for tr_type_idx, (now_tr_type, tr_type_description) in enumerate(tr_types):
	for site_end, site_end_description in site_ends:
		legend = []
		color = dict()
		for dataset_idx, dataset in enumerate(all_datasets):
			base_pattern = "solid"
			tr_type = "protein_coding" if (now_tr_type == "lncRNA" and dataset == "Random") else now_tr_type

			if dataset == "MANE" and now_tr_type == "lncRNA":
				tr_type = "protein_coding"
				base_pattern = "dotted"

			print(dataset, now_tr_type, tr_type)
			inMane = dataset == "MANE"
			temp_data = data[(data["dataset"] == dataset) & (data["site_type"] == site_end) & (data["gene_type"] == tr_type) & (data["inMANE"] == inMane)]
			if dataset != "Random" and dataset != "MANE":
				if now_tr_type == "protein_coding":
					dataset = dataset + "*"
				for conservation in [0, 1]:
					print(temp_data)
					now_df = temp_data[temp_data["conserved"] == conservation]
#					now_df = temp_data[(temp_data["phastCons"] == conservation)]
#					now_df = now_df[now_df.apply(lambda x: x["phastCons"] != x["conserved"], axis=1)]

					now_color = color[dataset] if dataset in color else ""
					color[dataset] = plot(now_df, pattern[conservation], now_color)
					legend.append(dataset + ", "  + conservation_caption[conservation])
			else:
				plot(temp_data, base_pattern, "")
				legend.append(dataset)

		if fig_idx % 2 == 0:
			plt.legend(legend, loc='upper left', ncol=2)

		plt.xlabel("Shift")
		plt.ylabel("gnomAD SNP rate")
		plt.xticks(range(len(shift)), labels=[tick(s, site_end) for s in shift])
		plt.title("SNP rate of certain positions around conserved and non-conserved \n" + site_end_description + " splice sites of " + tr_type_description + " genes from various datasets")

		plt.savefig("out/figure_snp_" + fig_letter[fig_idx] + ".pdf", bbox_inches="tight")
		plt.cla()
		fig_idx += 1

