import os
import sys
import gzip
import shutil
from Bio import SeqIO
from Bio.Seq import Seq

sys.path.append("src/lib")
from gtf_parse import getline
from gene_type import parse_type
from gene_type import get_type

def gen_introns(transcript, exons, out):
	tr_start = transcript.start
	size = transcript.end - tr_start + 1
	val = [0 for _ in range(size)]
	for e in exons:
		for i in range(e.start, e.end + 1):
			val[i - tr_start] = 1
	start = 0
	while start < size:
		end = start
		while end < size and val[end] == 0:
			end += 1
		if end > start:
			intron = [transcript.chr, transcript.src, "intron", str(tr_start + start), str(tr_start + end - 1), ".", transcript.strand, ".", transcript.attr_orig]
			print("\t".join(intron), file=out)
			start = end
		else:
			start += 1

gtf = sys.argv[1]
chr = sys.argv[2]
out_file = sys.argv[3]

gene_type = parse_type(gzip.open(gtf, 'rt'))
exon = dict()
transcript = dict()

pass_gene_types = ["protein_coding", "lncRNA"]
with gzip.open(gtf, 'rt') as handle:
	for record in getline(handle):
		if record.chr != chr:
			continue
		if record.type == "transcript":
			trid = record.attr["transcript_id"]
			if get_type(record, gene_type) in pass_gene_types:
				transcript[trid] = record
		elif record.type == "exon":
			trid = record.attr["transcript_id"]
			if get_type(record, gene_type) in pass_gene_types:
				if not trid in exon:
					exon[trid] = []
				exon[trid].append(record)

out_handle = open(out_file, "w")
for trid, record in transcript.items():
	gen_introns(record, exon[trid], out_handle)
out_handle.close()
