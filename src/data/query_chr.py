import os
import sys
from Bio import AlignIO
from Bio.AlignIO import MafIO
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

introns = sys.argv[1]
maf = sys.argv[2]
index_path = sys.argv[3]
chr = sys.argv[4]
out = sys.argv[5]

if os.path.isfile(index_path):
	os.remove(index_path)

idx = MafIO.MafIndex(index_path, maf, "hg38." + chr)
shift = 30

def to_seq(a, trid, idx, type):
	return [SeqRecord(Seq(r.seq), id="{}${}${}#{}".format(trid, idx, type, r.id), description="") for r in a]

sys.path.append("src/lib")
from gtf_parse import getline

out_handle = open(out, "w")
prev_tr = ""
for rec in getline(open(sys.argv[1])):
	tr = rec.attr["transcript_id"]
	if tr != prev_tr:
		intron_idx = 0
		prev_tr = tr

	start = rec.start - 1
	end = rec.end

	if rec.strand == '+':
		start_center = start
		end_center = end - 2

		multiple_alignment = idx.get_spliced([start_center - shift], [start_center + shift + 2], +1)
		multiple_alignment = to_seq(multiple_alignment, tr, intron_idx, "d")
		SeqIO.write(multiple_alignment, out_handle, "fasta")

		multiple_alignment = idx.get_spliced([end_center - shift], [end_center + shift + 2], +1)
		multiple_alignment = to_seq(multiple_alignment, tr, intron_idx, "a")
		SeqIO.write(multiple_alignment, out_handle, "fasta")
	else:
		start_center = end - 2
		end_center = start

		multiple_alignment = idx.get_spliced([start_center - shift], [start_center + shift + 2], -1)
		multiple_alignment = to_seq(multiple_alignment, tr, intron_idx, "d")
		SeqIO.write(multiple_alignment, out_handle, "fasta")

		multiple_alignment = idx.get_spliced([end_center - shift], [end_center + shift + 2], -1)
		multiple_alignment = to_seq(multiple_alignment, tr, intron_idx, "d")
		SeqIO.write(multiple_alignment, out_handle, "fasta")

	intron_idx += 1

os.remove(index_path)

