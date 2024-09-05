import os
import sys
from Bio import AlignIO
from Bio.AlignIO import MafIO

dir = sys.argv[1]
maf = sys.argv[2]
chr = sys.argv[3]

index_path = os.path.join(dir, "index", chr)
idx = MafIO.MafIndex(index_path, maf, "hg38." + chr)
shift = 30

sys.path.append("src/lib")
from gtf_parse import getline

query_dir = os.path.join(dir, "query", chr)
try:
    	shutil.rmtree(query_dir)
except:
       	pass

os.mkdir(query_dir)

prev_tr = ""
for rec in getline(open(os.path.join(dir, "introns_all", chr))):
	tr = rec.attr["transcript_id"]
	if tr != prev_tr:
		intron_idx = 0

	start = rec.start - 1
	end = rec.end

	if rec.strand == '+':
		start_center = start
		end_center = end - 2

		multiple_alignment = idx.get_spliced([start_center - shift], [start_center + shift + 2], +1)
		AlignIO.write(multiple_alignment, "$".join((query_dir + "/" + tr, str(intron_idx), "d")), "fasta")
		multiple_alignment = idx.get_spliced([end_center - shift], [end_center + shift + 2], +1)
		AlignIO.write(multiple_alignment, "$".join((query_dir + "/" + tr, str(intron_idx), "a")), "fasta")
	else:
		start_center = end - 2
		end_center = start

		multiple_alignment = idx.get_spliced([start_center - shift], [start_center + shift + 2], -1)
		AlignIO.write(multiple_alignment, "$".join((query_dir + "/" + tr, str(intron_idx), "d")), "fasta")
		multiple_alignment = idx.get_spliced([end_center - shift], [end_center + shift + 2], -1)
		AlignIO.write(multiple_alignment, "$".join((query_dir + "/" + tr, str(intron_idx), "a")), "fasta")
	intron_idx += 1

os.remove(index_path)

