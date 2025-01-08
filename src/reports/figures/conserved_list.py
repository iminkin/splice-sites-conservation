import os
import sys
import math
import pandas as pd
from Bio import SeqIO
sys.path.append("/home/iminkin2/projects3/splice-sites-paper/gtf_parse_1")
from gtf_parse import getline

def get_type(dataset, gtf):
	gene_type = dict()
	trid_to_geneid = dict()
	transcript_type = dict()
	if dataset == "RefSeq":
		type_attr = "gene_biotype"
	else:
		type_attr = "gene_type"

	for rec in getline(open(gtf)):
		if rec.type == "gene":
			gene_id = rec.attr["gene_id"]
			gene_type[gene_id] = rec.attr[type_attr]

		if rec.type == "transcript":
			transcript_id = rec.attr["transcript_id"]
			if "gene_id" in rec.attr:
				gene_id = rec.attr["gene_id"]
				if gene_id in gene_type:
					transcript_type[transcript_id] = gene_type[gene_id]
				else:
					transcript_type[transcript_id] = rec.attr[type_attr]
			elif type_attr in rec.attr:
				transcript_type[transcript_id] = rec.attr[type_attr]

        return transcript_type


base = "/home/iminkin2/projects3/splice-sites-paper-final/data/db"
for entry in ["refseq", "gencode", "chess"]:
	db_path = os.path.join(base, entry)
	if os.path.isdir(db_path):
		for file in os.listdir(db_path):
			if ".gtf" in file:
				gtf = file
		title = open(os.path.join(db_path, "title.txt").readline().strip()
		type = get_type(title, os.path.join(db_path, gtf)

		transcripts_path = os.path.join(db_path, "transcripts")
		for chr in os.listdir(transcripts_path):
			chr_path = os.path.join(transcripts_path, chr)
			for trid in os.listdir(chr_path):
				trid_path = os.path.join(chr_path, trid)
				for line in getline(open(trid_path)):
