import os
import sys
import math
import random
import gzip
from Bio import SeqIO
sys.path.append("src/lib")
from gtf_parse import getline

suffix = ["_d"]
alpha = "ACGT"

mane_path = sys.argv[1]
all_introns = dict()
for root, dirs, files in os.walk(mane_path):
	for file in files:
		intron_idx = 0
		h = open(os.path.join(root, file))
		for rec in getline(h):
			if rec.strand == '+':
				if not rec.chr in all_introns:
					all_introns[rec.chr] = set()
				all_introns[rec.chr].add((rec.start - 1, rec.end))

top_str = "N" * 20

cnt = dict()
mol = dict()

genome_stats_path = sys.argv[2]
for line in open(genome_stats_path):
	line = line.strip()
	if line[0] != '#':
		line = line.split('\t')
		mol[line[6]] = line[-1]

limit = int(sys.argv[4])
motif_idx = 0
human_path = sys.argv[3]
for seq_record in SeqIO.parse(gzip.open(human_path, "rt"), "fasta"):
	if not seq_record.id in mol:
		continue

	chr = mol[seq_record.id]
	if not chr in all_introns:
		continue

	site = set()
	introns = all_introns[chr]
	for (a, d) in introns:
		a += 100
		d -= 100
		if a > d:
			continue
		for _ in range(10):
			if motif_idx > limit:
				break

			br = random.randint(a, d)
			cr = random.randint(a, d)
			b = min(br, cr)
			c = max(br, cr)
			if c - b < 10:
				continue

			if (b in site) or (c in site):
				continue

			site.add(b)
			site.add(c)

			attr = 'transcript_id "TR.' + str(motif_idx) + '"; gene_type "protein_coding"'

			trans_record =  [chr, "RND", "transcript", str(a + 1), str(d + 1), ".", "+", ".", attr]
			exon_record_1 = [chr, "RND", "exon",       str(a + 1), str(b + 1), ".", "+", ".",  attr]
			exon_record_2 = [chr, "RND", "exon",       str(c + 1), str(d + 1), ".", "+", ".",  attr]
			for rec in [trans_record, exon_record_1, exon_record_2]:
				print("\t".join(rec))
			motif_idx += 1

