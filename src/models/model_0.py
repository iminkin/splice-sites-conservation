import os
import numpy as np
import pandas as pd

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score, roc_curve, recall_score, precision_score, f1_score
from sklearn.model_selection import train_test_split
motif = 30

def clipped_range(a, b):
	return ["cons_" + str(i - motif) for i in range(a, b) if i != 30 and i != 31] + ["cons_GTAG"]

feature = {"d": ["cons_GTAG"], "a" : ["cons_GTAG"]}
print(feature)

def train_model(data, end):
	data = data[data["site_type"] == end]
	data = data[((data["dataset"] == "MANE"))  | (data["dataset"] == "Random")]
	print(data.head())

	y = data["inMANE"].values
	X = data[feature[end]].values

	X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
	regressor = LogisticRegression(max_iter=50000)
	regressor.fit(X_train, y_train)

	y_predict_proba = regressor.predict_proba(X_test)[::,1]

	decision = 0.5
	print("Decision=" + str(decision))
	auc = roc_auc_score(y_test, y_predict_proba)
	fpr, tpr, thr = roc_curve(y_test,  y_predict_proba)

	roc_handle = open("../data/model/roc_0_" + end + ".txt", "w")
	print(auc, file=roc_handle)
	for fp, tp, th in zip(fpr, tpr, thr):
		print(fp, tp, th, file=roc_handle)

	y_test_pred = [(1 if p >= decision else 0) for p in y_predict_proba]
	accuracy = [regressor.score(X_train, y_train), regressor.score(X_test, y_test)]
	recall = recall_score(y_test, y_test_pred)
	precision = precision_score(y_test, y_test_pred)
	f1 = f1_score(y_test, y_test_pred)

	print("Intercept=" + str(regressor.intercept_))
	print("AUC=" + str(auc), regressor.coef_[0])
	print("Recall=" + str(recall), "Precision=" + str(precision), "F1=" + str(f1))
	print(f'Confusion Matrix: \n{confusion_matrix(y_test, y_test_pred)}')
	return regressor, decision

def run_model(data, end):
	new_data = data[data["site_type"] == end]
	all_X = new_data[feature[end]].values
	predict_all = regressor.predict_proba(all_X)[::, 1]
	new_data = new_data.assign(prob=predict_all)
	conserved = [(1 if p >= decision else 0) for p in predict_all]
	new_data = new_data.assign(conserved=conserved)

	return new_data


end_labels = {"a": "Acceptor", "d": "Donor"}
type_labels = {"protein_coding": "Coding", "lncRNA": "lncRNA"}

out = []
data = result = pd.read_csv("../data/db/splice_sites.csv")

for end, end_label in end_labels.items():
	regressor, decision = train_model(data, end)
	out_data = run_model(data, end)
	out.append(out_data)

#pd.concat(out).to_csv("../data/model/model_out_0.csv", index=False)
pd.concat(out).to_csv("../data/model/model_out.csv", index=False)
