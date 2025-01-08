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
	for label, fpr, tpr, thr, auc in data:
		curve.append(ax.plot(fpr, tpr, label=label))
		legend.append(label)
	ax.legend(loc=4)
	plt.gca().set_ylim([0.8, 1.0])
	fig.set_size_inches(8.5, 8.5)
	plt.title("Receiver operating characteristics")
	plt.savefig("out/roc.pdf")

data = []
base = "../../../data/processed/roc"

for end_label, end in [("Donor", "d"), ("Acceptor", "a")]:
	file = os.path.join(base, "roc_" + end + ".txt")
	handle = open(file)
	auc = float(handle.readline())
	fpr, tpr, thr = [], [], []

	for line in handle:
		line = line.split()
		fpr.append(float(line[0]))
		tpr.append(float(line[1]))
		thr.append(float(line[2]))
	data.append((end_label + " splice sites, AUC=%.2f" % auc, fpr, tpr, thr, auc))

make_roc_plots(data)
