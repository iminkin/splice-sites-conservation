import os
import sys
import math
import gzip
import shutil
import bisect
from Bio import AlignIO
from Bio import SeqIO
from collections import namedtuple
sys.path.append("src/lib")
from gtf_parse import getline

shift = 30

hg38_path = sys.argv[1]
hg38_stats = sys.argv[2]
mane_path = sys.argv[3]
query_path = sys.argv[4]
jobs_file = sys.argv[5]
genome = jobs_file.split("/")[-1]
Coord = namedtuple('Coord', ['human_pos', 'chr', 'pos', 'strand'])

def get_maps(line):
	ret = dict()
	now_map = line.attr["map"].split(",")
	if len(now_map) > 1:
		for now_gen_map in now_map:
			point = []
			genome,left,right = now_gen_map.split(":")
			for point_coord in (left, right):
				coord = point_coord.split("&")
				coord_tuple = Coord(int(coord[0]),coord[1],int(coord[2]),int(coord[3]))
				point.append(coord_tuple)
			human_start, human_end = line.start, line.end
			target_start, target_end = point[0].pos, point[1].pos
			ret[genome] = point
	return ret

def translate(index_ann, query_pos, idx):
	pos = index_ann[idx]
	idx = bisect.bisect_left(query_pos, pos)
	assert(query_pos[idx] == pos)
	return idx

def cutline(line):
	line = line.split(";")[:-2]
	return ";".join(line)

def range(line):
        return " ".join((line.strand, str(line.start), str(line.end)))


def generate_jobs(mane_path, query_path, seq, prefix, genome, jobs_file):
	strand = {"+" : +1, "-" : -1}
	mane_range = set()
	for trid in os.listdir(mane_path):
		for line in getline(open(os.path.join(mane_path, trid))):
			mane_range.add(range(line))
	mane = []
	query = []
	all_genome = { genome : [] }
	for line in getline(open(query_path)):
		gene_map = get_maps(line)
		if range(line) in mane_range:
			mane.append((line, gene_map))
			all_genome[genome].append(line.start)
		else:
			query.append((line, gene_map))

	mane.sort(key=lambda x: x[0].start)
	mane_pos = [x[0].start for x in mane]

	cnt = 0
	handle = open(jobs_file, "w")
	for idx, (line, gene_map) in enumerate(query):
		missing = 0
		collinear = 0
		for genome, ann_index in all_genome.items():
			if not genome in gene_map.keys():
				missing += 1
				new_idx = bisect.bisect_left(ann_index, line.start)
				if new_idx > 0 and new_idx < len(ann_index):
					left_ortholog = translate(ann_index, mane_pos, new_idx - 1)
					right_ortholog = translate(ann_index, mane_pos, new_idx)
					left_line = mane[left_ortholog][0]
					right_line = mane[right_ortholog][0]
					left_gap = line.start - left_line.start
					right_gap = right_line.start - line.start
					if left_gap < threshold and right_gap < threshold and genome in mane[left_ortholog][1] and genome in mane[right_ortholog][1]:
						left_ortholog_coord = mane[left_ortholog][1][genome][1]
						right_ortholog_coord = mane[right_ortholog][1][genome][0]
						target_strand = (left_ortholog_coord.strand, right_ortholog_coord.strand)
						target_strand_minus = (-target_strand[0], -target_strand[1])
						cnt += 1
						human_strand = (strand[left_line.strand], strand[right_line.strand])
						if left_ortholog_coord.chr == right_ortholog_coord.chr:
							collinear += 1
							chr = genome
							print(cutline(line.line), file=handle)
							pos = (line.start - prefix, line.end + prefix)
							now_seq = seq[pos[0] - 1:pos[1]]
							print(now_seq, file=handle)
							print(cutline(left_line.line), file=handle)
							print(cutline(right_line.line), file=handle)
							print(left_ortholog_coord.chr, left_ortholog_coord.pos, left_ortholog_coord.strand, right_ortholog_coord.pos, right_ortholog_coord.strand, file=handle)
							print("", file=handle)
							handle.close()


mol = dict()
for line in open(sys.argv[2]):
	line = line.strip()
	if line[0] != '#':
		line = line.split('\t')
		if line[-1] != "na":
			mol[line[6]] = line[-1]

threshold = 100000
hg38_handle = gzip.open(hg38_path, mode='rt')
for record in SeqIO.parse(hg38_handle, "fasta"):
	record_id = mol[record.id]
	mane_path_now = os.path.join(mane_path, record_id)
	query_path_now = os.path.join(query_path, record_id)
	if os.path.isdir(mane_path_now) and os.path.isfile(query_path_now):
		generate_jobs(mane_path_now, query_path_now, record.seq, shift, genome, jobs_file)

