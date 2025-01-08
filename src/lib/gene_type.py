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

def parse_transcript_coords(handle):
	coords = dict()
	for line in getline(handle):
		if line.type == "transcript":
			attr = line.attr
			coords[attr["transcript_id"]] = (line.chr, line.strand, line.start, line.end)
	return coords

def get_type(line, gene_type):
	if "transcript_type" in line.attr:
		return line.attr["transcript_type"]
	if "gene_type" in line.attr:
		return line.attr["gene_type"]

	gene_type = gene_type[line.attr["gene_id"]]
	if gene_type == "protein_coding":
		if "transcript_biotype" in line.attr:
			if line.attr["transcript_biotype"] == "mRNA":
				return gene_type
		else:
			return gene_type

	if gene_type == "lncRNA":
		if "transcript_biotype" in line.attr:
			if line.attr["transcript_biotype"] == "lnc_RNA":
				return gene_type
		else:
			return gene_type
	return "NA"

