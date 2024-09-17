import os
import sys
sys.path.append("src/lib")
from gtf_parse import getline


def intron(strand, start, end):
	return " ".join((strand, str(start), str(end)))

for db in sys.argv[1:]:
	seen = set()
	if os.path.isfile(db):
		for rec in getline(open(db)):
			signature = intron(rec.strand, rec.start, rec.end)
			if not signature in seen:
				seen.add(signature)
				print(rec.line)

