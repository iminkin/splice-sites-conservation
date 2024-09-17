import os
import sys
import gzip
import shutil
from Bio import SeqIO
from Bio.Seq import Seq

sys.path.append("src/lib")
from gtf_parse import getline

def parse_type(handle):
	gene_type = dict()
	for line in getline(handle):
		if line.type == "gene" or (line.type == "transcript" and "gene_id" in line.attr):
			attr = line.attr
			if "gene_biotype" in attr:
				gene_type[attr["gene_id"]] = attr["gene_biotype"]
			elif "gene_type" in attr:
				gene_type[attr["gene_id"]] = attr["gene_type"]
	return gene_type

def get_type(line, gene_type):
	if "gene_type" in line.attr:
		return line.attr["gene_type"]
	return gene_type[line.attr["gene_id"]]

