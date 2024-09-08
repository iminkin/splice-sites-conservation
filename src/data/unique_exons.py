import os
import sys
sys.path.append("src/lib")
from gtf_parse import getline


def intron(strand, start, end):
	return " ".join((strand, str(start), str(end)))

seen = set()
for dir in sys.argv[1:]:
	for file in os.listdir(dir):
		for rec in getline(open(os.path.join(dir, file))):
			signature = intron(rec.strand, rec.start, rec.end)
			if not signature in seen:
				seen.add(signature)
				print(rec.line)

