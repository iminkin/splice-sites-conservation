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
			if record.src == "RND":
				human_chr = record.chr
				if not human_chr in recovered:
					recovered[human_chr] = 0
				if ratio >= threshold:
					recovered[human_chr] += 1
			buffer = []
		else:
			buffer.append(line)

def is_complete_chr(chr):
	return not "_" in chr


gene_type = dict()
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
		if record.src != "RND":
			continue
		if not chr in exons_covered:
			exons_missed[chr] = 0
			exons_covered[chr] = 0
		genome_map_list = record.attr["map"].split(",")
		genome_set = set()
		for now_map in genome_map_list:
			genome = now_map.split(":")[0]
			genome_set.add(genome)
		for genome in all_genomes:
			if genome in genome_set:
				exons_covered[chr] += 1
			else:
				exons_missed[chr] += 1

def fn(n):
	return "{:,}".format(n)

all_chr = list(exons_covered.keys())
all_chr.sort()
for chr in all_chr:
	if not "_" in chr:
		if not chr in exons_recovered:
			exons_recovered[chr] = 0
		val = [chr]
		for d in [exons_covered, exons_missed, exons_recovered]:
			val.append(fn(d[chr]))

		print("&".join(val))
