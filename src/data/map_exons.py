import os
import sys
import math
import bisect
from Bio import AlignIO
from multiprocessing import Pool
sys.path.append("src/lib")
from gtf_parse import getline

maf_path = sys.argv[1]
exons_path = sys.argv[2]
out_path = sys.argv[3]

def fm(n):
	return "{:,}".format(n)

all_elements = []
for line in getline(open(exons_path)):
	all_elements.append(line)


elements_pointer = 0
all_elements.sort(key=lambda x: x.start)
all_elements_map_left = [dict() for _ in all_elements]
all_elements_map_right = [dict() for _ in all_elements]
all_elements_start = [l.start for l in all_elements]
all_elements_count = [dict() for l in all_elements]

sweep = []

def species(sid):
	return sid.split(".")[0]

match_cnt = 0
for bidx, block in enumerate(AlignIO.parse(maf_path, "maf")):
	genome = [row.id.split(".")[:2] for row in block]
	pos = [row.annotations["start"] for row in block]
	inc = [row.annotations["strand"] for row in block]
	assert(genome[0][0] == "hg38" and block[0].annotations["strand"] == 1)
#	print(bidx, match_cnt, file=sys.stderr)
	seen = False
	for j in range(0, len(block[0].seq)):
		if block[0].seq[j] != '-':
			pos0 = pos[0] + 1
			new_sweep = [e for e in sweep if pos0 <= all_elements[e].end]
			while elements_pointer < len(all_elements) and pos0 == all_elements[elements_pointer].start:
				new_sweep.append(elements_pointer)
				elements_pointer += 1
			sweep = new_sweep
			for e in sweep:
				assert(pos0 >= all_elements[e].start and pos0 <= all_elements[e].end)
				match_cnt += 1
				for i in range(1, len(block)):
					target_genome, target_chr = genome[i]
					if block[i].seq[j] != '-':
						now_coord = (pos0, target_chr,  pos[i], inc[i])
						if not target_genome in all_elements_map_left[e]:
							all_elements_map_left[e][target_genome] = all_elements_map_right[e][target_genome] = now_coord
						else:
							coord_left = all_elements_map_left[e][target_genome]
							coord_right = all_elements_map_right[e][target_genome]
							if pos0 < coord_left[0]:
								all_elements_map_left[e][target_genome] = now_coord
							if pos0 > coord_right[0]:
								all_elements_map_right[e][target_genome] = now_coord

					if block[0].seq[j] == block[i].seq[j]:
						now_coord = (pos0, target_chr,  pos[i], inc[i])
						if not target_genome in all_elements_count[e]:
							all_elements_count[e][target_genome] = 0
						all_elements_count[e][target_genome] += 1

		for i in range(0, len(block)):
			if block[i].seq[j] != '-':
				pos[i] += inc[i]
#	if match_cnt > 1000:
#		break

out_file = open(out_path, "w")
for idx, line in enumerate(all_elements):
	line = line.line
	if line[-1] != ";":
		line = line + ";"

	line = line + ' map "'
	record = []
	for genome, coord_left in all_elements_map_left[idx].items():
		coord_right = all_elements_map_right[idx][genome]
		now_record = [genome]
		for coord in (coord_left, coord_right):
			now_record.append("&".join((str(x) for x in coord)))
		record.append(":".join(now_record))
	line = line + ",".join(record) + '"; id "'
	record = []
	for genome, cnt in all_elements_count[idx].items():
		record.append(genome + ":" + str(cnt))
	print(line + ",".join(record) + '";', file=out_file)

