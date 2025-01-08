import os
import re
import sys
import gzip
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from functools import partial
sys.path.append("../../lib")
from gtf_parse import getline
from gtf_parse import getline_str
from gene_type import parse_type
from gene_type import get_type

#skip = "HLsarHar2 HLpseCor1 HLpseOcc1 HLpseCup1 HLgraAgi1 HLmunMun1 HLantFla1 HLtriVul1 monDom5".split()
skip = []

annotations_dir = "../../../data/raw/annotation"
mapped_exons_dir = "../../../data/interim/realignment/mapped_exons"
threshold_results_dir = "../../../data/interim/realignment/identity"
alignment_results_dir = "../../../data/interim/realignment/results"

exons_missed = dict()
exons_covered = dict()
exons_recovered = dict()
allowed_types = ["protein_coding", "lncRNA"]

def parse_report(report_path, threshold, recovered, gene_type):
	buffer = []
	for line in open(report_path):
		line = line.strip()
		if line == "":
			record = getline_str(buffer[1])
			human_seq = buffer[2]
			match_cnt = buffer[3].count("|")
			human_cnt = len(human_seq) - human_seq.count("-")
			ratio = float(match_cnt) / human_cnt
			if record.src != "RND":
				human_chr = record.chr
				if not human_chr in recovered:
					recovered[human_chr] = {"protein_coding" : 0, "lncRNA" : 0}
				if ratio >= threshold:
					now_type = get_type(record, gene_type)
					if now_type in recovered[human_chr]:
						recovered[human_chr][now_type] += 1
			buffer = []
		else:
			buffer.append(line)

def is_complete_chr(chr):
	return not "_" in chr


gene_type = dict()
for ann in os.listdir(annotations_dir):
	if ".gtf.gz" in ann:
		gtf = os.path.join(annotations_dir, ann)
		now_gene_type = parse_type(gzip.open(gtf, 'rt'))
		for gene, type in now_gene_type.items():
			gene_type[gene] = type

all_genomes = set()
for genome in os.listdir(threshold_results_dir):
	id_file = open(os.path.join(threshold_results_dir, genome))
	identity = id_file.readline().strip().split()
	if identity == []:
		continue
#	print(genome, identity)
	avg, std = float(identity[0]), float(identity[1])
	threshold = avg - std
	report = os.path.join(alignment_results_dir, genome)
	if os.path.isfile(report) and not genome in skip:
		parse_report(report, threshold, exons_recovered, gene_type)
		all_genomes.add(genome)


for chr in os.listdir(mapped_exons_dir):
	for record in getline(open(os.path.join(mapped_exons_dir, chr))):
		if record.src == "RND":
			continue
		now_type = get_type(record, gene_type)
		if not now_type in allowed_types:
			continue
		if not chr in exons_covered:
			exons_missed[chr] =  {"protein_coding" : 0, "lncRNA" : 0}
			exons_covered[chr] = {"protein_coding" : 0, "lncRNA" : 0}
		genome_map_list = record.attr["map"].split(",")
		genome_set = set()
		for now_map in genome_map_list:
			genome = now_map.split(":")[0]
			genome_set.add(genome)
		for genome in all_genomes:
			if genome in genome_set:
				exons_covered[chr][now_type] += 1
			else:
				exons_missed[chr][now_type] += 1

def fn(n):
	return "{:,}".format(n)


print(exons_covered)
print(exons_missed)
print(exons_recovered)

all_chr = list(exons_covered.keys())
all_chr.sort()
for chr in all_chr:
	if not "_" in chr:
		if not chr in exons_recovered:
			exons_recovered[chr] = {"protein_coding" : 0, "lncRNA" : 0}
		val = [chr]
		for d in [exons_covered, exons_missed, exons_recovered]:
			for now_type in allowed_types:
				val.append(fn(d[chr][now_type]))

		print("&".join(val))
