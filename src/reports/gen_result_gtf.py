import os
import sys
import gzip
import shutil
import pandas as pd
sys.path.append("src/lib")
from gtf_parse import getline
from gene_type import parse_type
from gene_type import get_type
from gene_type import parse_transcript_coords

gtf = sys.argv[1]
transcripts = sys.argv[2]
status = sys.argv[3]
output = sys.argv[4]

out_tr = set()
df = pd.read_csv(transcripts)
for _, row in df.iterrows():
	if row["well_supported"] == int(status):
		out_tr.add(row["transcript_id"])
output_handle = open(output, "w")
for record in getline(gzip.open(gtf, "rt")):
	if record.type == "gene" or ("transcript_id" in record.attr and record.attr["transcript_id"] in out_tr):
		print(record.line, file=output_handle)

output_handle.close()
