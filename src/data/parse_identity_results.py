import os
import sys
import shutil
import statistics
sys.path.append("src/lib")
from gtf_parse import getline

result = []
exons_dir = sys.argv[1]
out_file = sys.argv[2]
genome = out_file.split("/")[-1]
for chr in os.listdir(exons_dir):
	for line in getline(open(os.path.join(exons_dir, chr))):
		match = line.attr["id"].split(",")
		length = float(line.end - line.start + 1)
		if line.src != "RND":
			for p in match:
				p = p.split(":")
				if len(p) == 1 or p[0] != genome:
					continue
				now_id = int(p[1]) / length
				if now_id > 0:
					result.append(now_id)

handle = open(out_file, "w")
avg = sum(result) / len(result)
stdev = statistics.stdev(result)
print(avg, stdev, file=handle)
