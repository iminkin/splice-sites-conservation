import os
import sys
from Bio import SeqIO

main_dir = sys.argv[1]
chr = sys.argv[2]

query_dir = os.path.join(main_dir, "query", chr)
query_all = os.path.join(main_dir, "query_all", chr)
out_handle = open(query_all, "w")

for site in os.listdir(query_dir):
	in_handle = open(os.path.join(query_dir, site))
	for record in SeqIO.parse(in_handle, "fasta"):
		record.id = site + " " + record.id
		SeqIO.write(record, out_handle, "fasta")
