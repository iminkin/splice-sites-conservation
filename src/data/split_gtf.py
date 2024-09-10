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

transcript = dict()


gene_type = dict()
with gzip.open(gtf, 'rt') as handle:
	for line in getline(handle):
		if line.type == "gene":
			attr = line.attr
			if "gene_biotype" in attr:
				gene_type[attr["gene_id"]] = attr["gene_biotype"]
			elif "gene_type" in attr:
				gene_type[attr["gene_id"]] = attr["gene_type"]


pass_gene_types = ["protein_coding", "lncRNA"]
with gzip.open(gtf, 'rt') as handle:
	for record in getline(handle):
		chr, type, strand = record.chr, record.type, record.strand
		if not chr in chromosomes:
			continue
		ex_path = ex_dir + chr
		tr_path = tr_dir + chr
		if not os.path.isdir(ex_path):
			os.mkdir(ex_path)
			os.mkdir(tr_path)

		attr = record.attr
		if type == "transcript":
			trid = record.attr["transcript_id"]
			if ("gene_id" in attr and gene_type[attr["gene_id"]] in pass_gene_types) or ("gene_type" in attr and attr["gene_type"] in pass_gene_types):
				h = open(tr_path + "/" + trid, "w")
				print(record.line, file=h)
				h.close()
		elif type == "exon":
			trid = record.attr["transcript_id"]
			if ("gene_id" in attr and gene_type[attr["gene_id"]] in pass_gene_types) or ("gene_type" in attr and attr["gene_type"] in pass_gene_types):
				h = open(ex_path + "/" + trid, "a")
				print(record.line, file=h)
				h.close()


