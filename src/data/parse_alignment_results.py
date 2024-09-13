import os
import re
import sys
import gzip
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from functools import partial

def parse_report(hg38, report_path, threshold, out_file):
	buffer = []
	cons = dict()
	report_path_parts = report_path.split("/")
	for line in open(report_path):
		line = line.strip()
		if line == "":
			line = buffer[1].split()
			human_chr, human_start, human_end = line[0], int(line[3]), int(line[4])
			human_chr_seq = hg38[human_chr]

			human_seq = buffer[2]
			match_cnt = buffer[3].count("|")
			human_cnt = len(human_seq) - human_seq.count("-")
			ratio = float(match_cnt) / human_cnt
			if ratio < threshold:
				buffer = []
				continue

			target_seq = buffer[4]
			result = buffer[5].split()
			query_strand = result[2]

			inc_pos = +1
			human_pos = human_start - prefix
			if query_strand == "-":
				inc_pos = -1
				human_pos = human_end + prefix
			for i in range(len(human_seq)):
				if human_seq[i] != "-":
					check_ch = human_chr_seq[human_pos - 1].upper()
					if query_strand == "-":
						check_ch = Seq(check_ch).complement()
					if human_seq[i] != check_ch:
						print(human_seq)
						print(target_seq)
						print(line)
						print(i, human_seq[i], check_ch, query_strand, human_pos)
					assert(human_seq[i] == check_ch)
					if human_seq[i] == target_seq[i]:
						if not human_chr in cons:
							cons[human_chr] = []
						cons[human_chr].append(human_pos)
					human_pos += inc_pos
			buffer = []
		else:
			buffer.append(line)

	for chr, val in cons.items():
		print(" ".join([chr] + [str(x) for x in val]), file=out_file)

mol = dict()
for line in open(sys.argv[2]):
	line = line.strip()
	if line[0] != '#':
		line = line.split('\t')
		if line[-1] != "na":
			mol[line[6]] = line[-1]

prefix = 30
hg38 = dict()
for record in SeqIO.parse(gzip.open(sys.argv[1], "rt"), "fasta"):
	hg38[mol[record.id]] = Seq(record.seq.upper())

id_file = open(sys.argv[4])
identity = id_file.readline().strip().split()
avg, std = float(identity[0]), float(identity[1])
out_file = open(sys.argv[5], "w")
parse_report(hg38, sys.argv[3], avg - std, out_file)
