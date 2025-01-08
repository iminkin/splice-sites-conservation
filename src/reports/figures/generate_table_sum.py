import os
import sys
import pandas as pd

#data = pd.read_csv("../../data/db/splice_sites.csv")
data = result = pd.read_csv("../../../data/processed/model_out.csv")
tr_data = result = pd.read_csv("../../../data/processed/transcripts.csv")
tr_data = tr_data[tr_data["total_sites"] > 0]

#test = tr_data[(tr_data["dataset"] == "CHESS 3") & (tr_data["inMANE"] == 1)]
#print(len(test.index))

ends = ["d", "a"]
types = [("protein_coding", "Protein-Coding"), ("lncRNA", "lncRNA")]

minus = "*"
datasets = ["MANE", "GENCODE", "RefSeq", "CHESS 3"]

def generate_row(df, table_row):
	n = len(df.index)
	table_row.append("{:,}".format(n))

for (type, type_title) in types:
	print("\\hline")
	print("\\textit{" + type_title +  "} & & & & & & \\\\")
	for title in datasets:
		now_title = title
		row = [now_title]
		inMANE = 1 if title == "MANE" else 0
		for cat in [0, 1]:
			for end in ends:
				if type == "protein_coding" and title != "MANE":
					row[0] = now_title + minus
				frame = data[(data["dataset"] == title) & (data["gene_type"] == type) & (data["site_type"] == end) & (data["inMANE"] == inMANE)]
				if cat == 1:
					if title != "MANE" and title != "Random":
						frame = frame[frame["conserved"] == 1]
					else:
						frame = pd.DataFrame()

				if not frame.empty:
					generate_row(frame, row)
				else:
					row.append("-")
		for cat in [0, 1]:
			frame = tr_data[(tr_data["dataset"] == title) & (tr_data["gene_type"] == type) & (tr_data["inMANE"] == inMANE)]
			if cat == 1:
				if title != "MANE" and title != "Random":
					frame = frame[frame["conservation_status"] == 1]
				else:
					frame = pd.DataFrame()

			if not frame.empty:
				generate_row(frame, row)
			else:
				row.append("-")


#		if row.count('-') < 4:
		print("&".join(row) + "\\\\")

print("\\hline")
print("\\textit{Synthetic data} & & & &\\\\")
row = ["Random"]
frame = data[(data["dataset"] == "Random") & (data["gene_type"] == "protein_coding") & (data["site_type"] == "d")]
generate_row(frame, row)
generate_row(frame, row)
print("&".join(row) + "& - & -\\\\")
print("")

