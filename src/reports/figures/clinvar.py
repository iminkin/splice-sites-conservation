import os
import math
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

clinvar = dict()
for line in open("/home/iminkin2/projects3/splice-sites-paper-final/data/clinvar/pathogenic.txt"):
	line = line.strip().split()
	if "pathogenic" in line[-1] or "likely_pathogenic" in line[-1]:
		if not line[0] in clinvar:
			clinvar[line[0]] = set()
		clinvar[line[0]].add(int(line[1]))

minus = "*"
inc = {"+" : +1, "-" : -1}
site_ends = [("d", "donor"), ("a", "acceptor")]
datasets = ["GENCODE", "RefSeq", "CHESS 3"]
types = [("protein_coding", "Protein-Coding"), ("lncRNA", "lncRNA")]
data = pd.read_csv("../../../data/processed/model_out.csv")

for (now_type, type_title) in types:
	print("\\hline")
	print("\\textit{" + type_title +  "} & & & & \\\\")

	for dataset in datasets:
		count = {"d" : [0, 0], "a" : [0, 0]}
		all_count = {"d" : [0, 0], "a" : [0, 0]}
		for site_end, site_end_description in site_ends:
			now_df = data[(data["dataset"] == dataset) & (data["gene_type"] == now_type) & (data["inMANE"] == 0) & (data["site_type"] == site_end)]
			for _, row in now_df.iterrows():
				cons = row["conserved"]
				chr, strand, pos = row["chr"], row["strand"], row["pos"]
				if chr in clinvar:
					if pos in clinvar[chr] or  pos + inc[strand] in clinvar[chr]:
						count[site_end][cons] += 1
				all_count[site_end][cons] += 1

		dataset_title = dataset
		if now_type == "protein_coding":
			dataset_title = dataset_title + "*"

		sites_count = []
		for cons in [1, 0]:
			for site_end, _ in site_ends:
				snp_count = count[site_end][cons]
				total_count = all_count[site_end][cons]
				ratio = float(snp_count) / total_count * 100
				ratio = "({:.2f}\%)".format(ratio)
				now_number = str(snp_count) + " " + ratio
				sites_count.append(now_number)

		print("&".join([dataset_title] + sites_count) + "\\\\")



