import os
import sys
import pandas as pd

csv = []
for db in sys.argv[1:-1]:
	table = pd.read_csv(db)
	csv.append(table)

result = pd.concat(csv)
result.to_csv(sys.argv[-1], index=False)
