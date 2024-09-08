import os
import sys
import math
import shutil
import bisect
from Bio import AlignIO
from Bio import SeqIO
from collections import namedtuple
sys.path.append("src/lib")
from gtf_parse import getline

exons_dir = sys.argv[1]
mane_exons_path = sys.argv[2]
genome = sys.argv[3]
jobs_file = sys.argv[4]
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

def generate_jobs(mane_path, query_path, gene_type, seq, prefix, jobs_dir):
	strand = {"+" : +1, "-" : -1}
	mane = []
	all_genome = dict()
	for line in getline(open(mane_path)):
		gene_map = get_maps(line)
		for genome in gene_map.keys():
			if not genome in all_genome:
				all_genome[genome] = []
		mane.append((line, gene_map))

	mane.sort(key=lambda x: x[0].start)
	mane_pos = [x[0].start for x in mane]

	for idx, (line, gene_map) in enumerate(mane):
		for genome in gene_map.keys():
			all_genome[genome].append(line.start)

	query = []
	for line in getline(open(query_path)):
		gene_map = get_maps(line)
		query.append((line, gene_map))

	cnt = 0
	for idx, (line, gene_map) in enumerate(query):
		gene_id = line.attr["gene_id"]
		now_gene_type = gene_type[gene_id]
		if now_gene_type in gene_types:
			missing = 0
			collinear = 0
			for genome, ann_index in all_genome.items():
				if not genome in gene_map.keys():
#					print(line.line)
#					print(genome, gene_map)
#					continue

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
								handle = open(os.path.join(jobs_dir, chr), "a+")
								print(cutline(line.line), file=handle)
								pos = (line.start - prefix, line.end + prefix)
								now_seq = seq[line.chr][pos[0] - 1:pos[1]]
								print(now_seq, file=handle)
								print(cutline(left_line.line), file=handle)
								print(cutline(right_line.line), file=handle)
								print(left_ortholog_coord.chr, left_ortholog_coord.pos, left_ortholog_coord.strand, right_ortholog_coord.pos, right_ortholog_coord.strand, file=handle)
								print("", file=handle)
								handle.close()


threshold = 100000

job_handle = dict()
job_handle[genome] = open(jobs_file)

mane_path = "exons_out/mane/"
query_path = os.path.join("exons_out", db)
query_gene_path = os.path.join("/home/iminkin2/projects3/splice-sites-paper-final/data/db/", db, "genes", "genes.gtf")

gene_types = ["protein_coding", "lncRNA"]
gene_type = dict()
for line in getline(open(query_gene_path)):
	attr = line.attr
	if db == "refseq":
		gene_type[attr["gene_id"]] = attr["gene_biotype"]
	else:
		gene_type[attr["gene_id"]] = attr["gene_type"]

seq = dict()
for record in SeqIO.parse("/home/iminkin2/projects3/splice-sites-paper-final/data/hg38/hg38.fa", "fasta"):
	seq[record.id] = record.seq
	mane_path_now = os.path.join(mane_path, record.id)
	query_path_now = os.path.join(query_path, record.id)
	if os.path.isfile(mane_path_now) and os.path.isfile(query_path_now):
		generate_jobs(mane_path_now, query_path_now, gene_type, seq, 30, jobs_dir)

