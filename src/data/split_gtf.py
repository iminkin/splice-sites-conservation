import os
import sys
import gzip
import shutil
from Bio import SeqIO
from Bio.Seq import Seq

sys.path.append("src/lib")
from gtf_parse import getline

gtf = sys.argv[1]
main_dir = sys.argv[2]
chromosomes = set(sys.argv[3].split())

ex_dir = os.path.join(main_dir, "exons/")
tr_dir = os.path.join(main_dir, "transcripts/")
introns_dir = os.path.join(main_dir, "introns")
introns_out_dir = os.path.join(main_dir, "introns_all")

try:
	shutil.rmtree(main_dir)
except:
	pass

os.mkdir(main_dir)
os.mkdir(ex_dir)
os.mkdir(tr_dir)
os.mkdir(introns_dir)
os.mkdir(introns_out_dir)
os.mkdir(os.path.join(main_dir, "query"))
os.mkdir(os.path.join(main_dir, "query_all"))
os.mkdir(os.path.join(main_dir, "index"))

for chr in chromosomes:
	os.mkdir(os.path.join(ex_dir, chr))
	os.mkdir(os.path.join(tr_dir, chr))

transcript = dict()

gene_type = dict()
with gzip.open(gtf, 'rt') as handle:
	for line in getline(handle):
		if line.type == "gene" or (line.type == "transcript" and "gene_id" in line.attr):
			attr = line.attr
			if "gene_biotype" in attr:
				gene_type[attr["gene_id"]] = attr["gene_biotype"]
			elif "gene_type" in attr:
				gene_type[attr["gene_id"]] = attr["gene_type"]

def get_type(line, gene_type):
	if "gene_type" in line.attr:
		return line.attr["gene_type"]
	return gene_type[line.attr["gene_id"]]

pass_gene_types = ["protein_coding", "lncRNA"]
with gzip.open(gtf, 'rt') as handle:
	for record in getline(handle):
		chr, type, strand = record.chr, record.type, record.strand
		if not chr in chromosomes:
			continue
		ex_path = ex_dir + chr
		tr_path = tr_dir + chr
		attr = record.attr
		if type == "transcript":
			if get_type(record, gene_type) in pass_gene_types:
				trid = record.attr["transcript_id"]
				h = open(tr_path + "/" + trid, "w")
				print(record.line, file=h)
				h.close()
		elif type == "exon":
			if get_type(record, gene_type) in pass_gene_types:
				trid = record.attr["transcript_id"]
				h = open(ex_path + "/" + trid, "a")
				print(record.line, file=h)
				h.close()


