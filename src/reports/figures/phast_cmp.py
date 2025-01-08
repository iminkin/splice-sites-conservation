import os
import math
import scipy.stats
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch

data = pd.read_csv("../../../data/processed/model_out.csv")
data = data[(data["inMANE"] == 0) & (data["dataset"] != "Random") & (data["phastCons_0"].notnull())]
print(data)
data["phastCons"] = data.apply(lambda x: min(x["phastCons_0"], x["phastCons_1"]), axis=1)
print(data["phastCons"].corr(data["prob"]))




