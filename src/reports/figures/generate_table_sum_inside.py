import os
import sys
import pandas as pd

#data = pd.read_csv("../../data/db/splice_sites.csv")
data = result = pd.read_csv("../../../data/processed/model_out.csv")
ends = ["d", "a"]
types = [("protein_coding", "Protein-Coding")]

minus = "*"
datasets = ["GENCODE", "RefSeq", "CHESS 3"]

def generate_row(df, table_row):
	n = len(df.index)
	table_row.append("{:,}".format(n))

mane_transcripts = 0
for (type, type_title) in types:
	print("\\hline")
	print("\\textit{" + type_title +  "} & & & & & & \\\\")
	for title in datasets:
		now_title = title + "*"
		row = [now_title]
		inMANE = 1 if title == "MANE" else 0
		for end in ends:
			for inside in [0, 1]:
				for cat in [0, 1]:
					frame = data[(data["dataset"] == title) & (data["gene_type"] == type) & (data["site_type"] == end) & (data["inMANE"] == inMANE) & (data["inside_mane_exon"] == inside) & (data["conserved"] == cat)]
					generate_row(frame, row)
		print("&".join(row) + "\\\\")


