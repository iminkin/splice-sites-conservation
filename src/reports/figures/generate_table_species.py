import os
import sys
import math
import pandas as pd

data = pd.read_csv("../../../data/processed/splice_sites_cons.csv")
data = data[(data["cons_GTAG"] < 100) & (data["inMANE"] == 0)]
types = [("protein_coding", "Protein-Coding"), ("lncRNA", "lncRNA")]
#print(data)

species_name = dict()
for line in open("species.txt"):
	line = line.strip().split(",")
	species_name[line[0]] = line[1]

all_datasets = ["GENCODE", "RefSeq", "CHESS 3", "Random"]
for (now_type, type_title) in types:
	local_count = dict()
	global_count = dict()
	total_sites = dict()
	df = data[data["gene_type"] == now_type]
	for _, row in df.iterrows():
		dataset = row["dataset"]
		now_end = row["site_type"]
		if dataset in all_datasets:
			if not dataset in local_count:
				local_count[dataset] = dict()
				total_sites[dataset] = {"d" : 0, "a" : 0}

			total_sites[dataset][now_end] += 1
			if not isinstance(row["cons_array"], str):
				continue
			for genome in row["cons_array"].split(";"):
				if genome != "hg38":
					if not genome in local_count[dataset]:
						local_count[dataset][genome] = {"d" : 0, "a" : 0}
					if not genome in global_count:
						global_count[genome] = 0

					local_count[dataset][genome][now_end] += 1
					if genome != "Random":
						global_count[genome] += 1

	species = [(count, genome) for (genome, count) in global_count.items()]
	species.sort(reverse=True)
	species = species[:30]
	row = ["\\textit{" + type_title + "}"] + [' ' for _ in all_datasets] * 2
	print("&".join(row) + "\\\\")
	for (_, genome) in species:
		if not genome in species_name:
			print(genome, file=sys.stderr)
			continue
		row = [species_name[genome]]
		for site_end in ["d", "a"]:
			for dataset in all_datasets:
				if dataset != "Random" or now_type == "protein_coding":
					ratio = float(local_count[dataset][genome][site_end]) / total_sites[dataset][site_end]
					row.append("{:.2f}\%".format(ratio * 100))
				else:
					row.append("-")
		print("&".join(row) + "\\\\")
	if now_type == "protein_coding":
		print("\\hline")
