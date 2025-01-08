import os
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

def make_roc_plots(data):
	fig, ax = plt.subplots()
	fsize = 18
	plt.tick_params(labelsize=fsize)
	matplotlib.rcParams.update({'font.size': fsize})
	ax.set_xlabel("False positive rate", fontsize=fsize)
	ax.set_ylabel("True positive rate", fontsize=fsize)

	curve = []
	legend = []
	for label, fpr, tpr, thr, auc, pattern, color in data:
		curve.append(ax.plot(fpr, tpr, pattern, label=label, color=color))
		legend.append(label)
	ax.legend(loc=4)
	plt.gca().set_ylim([0.8, 1.0])
	fig.set_size_inches(9.5, 9.5)
	plt.title("Receiver operating characteristics")
	plt.savefig("out/roc.pdf")

data = []
base = "../../../data/processed/roc"
color = ["orange", "green"]
prefix = ["roc_0_", "roc_"]
pr_title = [["GT only", "AG only "], ["full model", "full model "]]
pattern = ["--", "-"]

for pidx, p in enumerate(prefix):
	for eidx, (end_label, end) in enumerate([("Donor", "d"), ("Acceptor", "a")]):
		file = os.path.join(base, p + end + ".txt")
		handle = open(file)
		auc = float(handle.readline())
		fpr, tpr, thr = [], [], []

		for line in handle:
			line = line.split()
			fpr.append(float(line[0]))
			tpr.append(float(line[1]))
			thr.append(float(line[2]))
		data.append((end_label + " sites, " + pr_title[pidx][eidx] + ", AUC=%.3f" % auc, fpr, tpr, thr, auc, pattern[pidx], color[eidx]))

make_roc_plots(data)
