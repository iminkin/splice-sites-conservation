import os
import sys
import json
import gzip
import edlib
import shutil
from Bio import SeqIO
from Bio.Seq import Seq


def realign(job_path, genome_fagz, report_path):
	job_path_part = job_path.split("/")
	if not os.path.isfile(genome_fagz):
		return

	seq = dict()
	for record in SeqIO.parse(gzip.open(genome_fagz, "rt"), "fasta"):
		rid = record.id.split(".")[0]
		seq[rid] = record.seq

	alias = dict()
	refseq = dict()
#	if os.path.isfile(genome_js):
#		for line in open(genome_js):
#			mp = json.loads(line)
#			gid = mp["genbankAccession"].split(".")[0]
#			alias[gid] = mp["sequenceName"]
#			if "refseqAccession" in mp:
#				refseq[gid] = mp["refseqAccession"].split(".")[0]

	all = 0
	match = 0
	buffer = []
	out_file = open(report_path, "w")
	for line in open(job_path):
		line = line.strip()
		if line == "":
			coord_query = buffer[-1].split()
			chr, start, end = coord_query[0], int(coord_query[1]), int(coord_query[3])
			if chr in seq:
				now_seq = seq[chr]
			elif chr in alias and alias[chr] in seq:
				now_seq = seq[alias[chr]]
			elif chr in refseq and refseq[chr] in seq:
				now_seq = seq[refseq[chr]]
			else:
				buffer = []
				continue

			target = now_seq[start:end].upper()
			query = buffer[1].upper()

			lq = float(len(query))
			result = edlib.align(query=query, target=target, mode="HW", task="path")
			dist = result["editDistance"]
			strand = "+"

			rev_query = str(Seq(query).reverse_complement())
			rev_result = edlib.align(query=rev_query, target=target, mode="HW", task="path")
			if rev_result["editDistance"] < dist:
				result = rev_result
				query = rev_query
				strand = "-"

			all += 1
			ratio = float(result["editDistance"]) / len(query)
			if ratio < .5:
				loc = result["locations"][0]
				match += 1
				print(match, all, file=out_file)
				print(buffer[0], file=out_file)
				al = edlib.getNiceAlignment(result, query, target)
				al = (str(l) for l in al.values())
				print("\n".join(al), file=out_file)
				print(result["editDistance"], chr, strand, loc[0], loc[1], file=out_file)
				print("", file=out_file)

			buffer = []
		else:
			buffer.append(line)


realign(sys.argv[1], sys.argv[2], sys.argv[3])

