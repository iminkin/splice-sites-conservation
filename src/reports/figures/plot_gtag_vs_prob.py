import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


site_ends = [("d", "donor"), ("a", "acceptor")]
tr_types = [("protein_coding", "protein-coding"), ("lncRNA", "lncRNA")]
all_datasets = ["GENCODE", "RefSeq", "CHESS 3"]

data = result = pd.read_csv("../../../data/processed/model_out.csv")

fig_idx = 0
fig_letter = "abcd"
nn = 406
article = ["a", "an"]
for site_end, site_end_description in site_ends:
	legend = []
	color = []
	now_df = data[((data["dataset"] == "GENCODE") | (data["dataset"] == "RefSeq") | (data["dataset"] == "CHESS 3")) & (data["inMANE"] == 0) & (data["site_type"] == site_end)]
	sns.histplot(now_df, x="cons_GTAG", y="prob")
	plt.xlabel("Number of genomes")
	plt.ylabel("Probability")
	plt.title("Probability of " + article[fig_idx] + " " + site_end_description + ' splice sites being classified as \n "well-supported" and conservation of canonical dinucleotides')
	plt.savefig("out/figure_gtag_vs_prob_" + fig_letter[fig_idx] + ".pdf", bbox_inches="tight")
	plt.cla()
	fig_idx += 1
