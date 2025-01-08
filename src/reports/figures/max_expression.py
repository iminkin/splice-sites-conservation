import os
import math
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams

fig_letter = "abcd"
#cons = {"00" : "Neither site conserved", "10" : "Only donor conserved", "01" : "Only acceptor conserved", "11" : "Both conserved"}
cons = {"00" : "Neither site is well-suported", "10" : "Only donor is well-supported", "01" : "Only acceptor is well-supported", "11" : "Both are well-suported"}

all_tissues = []
coverage_path = '/ccb/salz8-1/avaraby/chess3_rerun_31102021/step1/'
for tissue in os.listdir(coverage_path):
	tissue_path = os.path.join(coverage_path, tissue)
	if os.path.isdir(tissue_path):
		all_tissues.append(tissue)

tr_types = [("protein_coding", "protein-coding"), ("lncRNA", "lncRNA")]
df = pd.read_csv("../../../data/processed/introns.csv")
df = df[df["inMANE"] == 0]
df["dataset_title"] = df.apply(lambda x: x["dataset"] + ("*" if x["gene_type"] == "protein_coding" else ""), axis=1)
df["donor_cons_res"] = df.apply(lambda x: x["donor_mane"] or x["donor_cons"], axis=1)
df["acceptor_cons_res"] = df.apply(lambda x: x["acceptor_mane"] or x["acceptor_cons"], axis=1)
df["conservation"] = df.apply(lambda x: cons[str(x["donor_cons_res"]) + str(x["acceptor_cons_res"])], axis=1)
df["max_coverage"] = df.apply(lambda x: np.nanmax(np.array(x[all_tissues], dtype="float32")), axis=1)
#df["max_coverage"] = df.apply(lambda x: max(x[all_tissues]), axis=1)

for tr_type_idx, (now_tr_type, tr_type_description) in enumerate(tr_types):
	now_df = df[(df["gene_type"] == now_tr_type)]
	now_df = now_df.melt(id_vars=["conservation", "dataset_title"], value_vars=["max_coverage"]).dropna()
	print(now_df)
	rcParams['figure.figsize'] = 10,8

	b = sns.boxplot(data=now_df, x="dataset_title", y="value", hue="conservation", hue_order=cons.values(), showfliers=False, whis=1.5)
	if tr_type_idx == 1:
		b.legend().remove()
	else:
		b.legend().set_title("Splice site support")
		sns.move_legend(b, "upper center")

	plt.ylabel("Coverage")
	plt.xlabel("Dataset")
	plt.title("Maximum coverage of introns from " + tr_type_description + " genes across \n different splice site support categories")
	plt.savefig("out/figure_max_coverage_" + fig_letter[tr_type_idx] + ".pdf", bbox_inches="tight")
	plt.cla()


