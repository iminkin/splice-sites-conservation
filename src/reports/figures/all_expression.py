import os
import math
import string
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams

fig_letter = "abcd"
abc = string.ascii_lowercase
#cons = {"00" : "Neither site conserved", "10" : "Only donor conserved", "01" : "Only acceptor conserved", "11" : "Both conserved"}
#cons = {"00" : "Neither site conserved", "10" : "Donor only", "01" : "Acceptor only", "11" : "Both conserved"}
cons = {"00" : "Neither site is well-suported", "10" : "Only donor is well-supported", "01" : "Only acceptor is well-supported", "11" : "Both are well-suported"}

all_tissues = []
coverage_path = '/ccb/salz8-1/avaraby/chess3_rerun_31102021/step1/'
for tissue in os.listdir(coverage_path):
	tissue_path = os.path.join(coverage_path, tissue)
	if os.path.isdir(tissue_path):
		all_tissues.append(tissue)

tr_types = [("protein_coding", "protein-coding"), ("lncRNA", "lncRNA")]
df = pd.read_csv("/home/iminkin2/projects3/splice-sites-paper-final/data/model/introns.csv")
df = df[df["inMANE"] == 0]
df["dataset_title"] = df.apply(lambda x: x["dataset"] + ("*" if x["gene_type"] == "protein_coding" else ""), axis=1)
df["conservation"] = df.apply(lambda x: cons[str(x["donor_cons"]) + str(x["acceptor_cons"])], axis=1)

for tr_type_idx, (now_tr_type, tr_type_description) in enumerate(tr_types):
	for tissue_idx, tissue in enumerate(all_tissues):
		now_df = df[(df["gene_type"] == now_tr_type)]
		now_df = now_df.melt(id_vars=["conservation", "dataset_title"], value_vars=[tissue]).dropna()
#		sns.set(rc={'figure.figsize':(10,8)})
		rcParams['figure.figsize'] = 10,8
		b = sns.boxplot(data=now_df, x="dataset_title", y="value", hue="conservation", hue_order=cons.values(), showfliers=False, whis=1.5)
		if tissue_idx != 0:
			b.legend().remove()
		else:
			b.legend().set_title("Splice site support")
#			sns.move_legend(b, "upper center")
		plt.ylabel("Coverage")
		plt.xlabel("Dataset")
		tissue_str = " ".join(tissue.split("_")).lower()
		plt.title("Coverage of introns from " + tr_type_description + " involving \n splice site of different support categories in " + tissue_str + " tissue")
		plt.savefig("appendix/coverage_" + now_tr_type + "_" + str(tissue_idx) + ".pdf", bbox_inches="tight")
		plt.cla()


