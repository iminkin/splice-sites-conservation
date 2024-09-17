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

gtf = sys.argv[1]
chr = sys.argv[2]
out_file = sys.argv[3]

gene_type = parse_type(gzip.open(gtf, 'rt'))
exon = dict()
transcript = dict()

out_handle = open(out_file, "w")
pass_gene_types = ["protein_coding", "lncRNA"]

if not os.path.isfile(gtf):
	exit()

with gzip.open(gtf, 'rt') as handle:
	for record in getline(handle):
		if record.chr != chr:
			continue
		elif record.type == "exon":
			trid = record.attr["transcript_id"]
			if get_type(record, gene_type) in pass_gene_types:
				print(record.line, file=out_handle)
