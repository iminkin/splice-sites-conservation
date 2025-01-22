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

def site(chr, strand, start):
	return "&".join((chr, strand, str(start)))

def donor_acceptor(chr, strand, start, end):
	if strand == '+':
		donor_coords = "&".join((chr, strand, str(start)))
		acceptor_coords = "&".join((chr, strand, str(end - 1)))
	else:
		donor_coords = "&".join((chr, strand, str(end)))
		acceptor_coords = "&".join((chr, strand, str(start + 1)))
	return (donor_coords, acceptor_coords)


all_types = ["protein_coding", "lncRNA"]
model_out = sys.argv[1]
out_file = sys.argv[2]
tr_out_file = sys.argv[3]
coverage_path = sys.argv[4]
db_list = sys.argv[5:]

coverage = dict()
for bed in os.listdir(coverage_path):
	bed_path = os.path.join(coverage_path, bed)
	tissue = bed.split(".")[0]
	coverage[tissue] = dict()
	handle = gzip.open(bed, "rt")
	handle.readline()
	now_coverage = coverage[tissue]
	for line in handle:
		line = line.strip().split("\t")
		chr, strand, start, end, cov = line[0], line[-1], int(line[1]), int(line[2]), int(line[-2])
		start += 1
		donor, acceptor = donor_acceptor(chr, strand, start, end)
		if not donor in now_coverage:
			now_coverage[donor] = dict()
		now_coverage[donor][acceptor] = cov

mane_sites = set()
logit_conserved = set()
phast_conserved = set()
df = pd.read_csv(model_out)
df_mane = df[df["inMANE"] == 1]

for _, row in df_mane.iterrows():
	mane_sites.add(site(row["chr"], row["strand"], row["pos"]))

df_non_mane = df[df["inMANE"] == 0]
df_logit = df[df["well_supported"] == 1]
df_phast = df_non_mane[df_non_mane.apply(lambda x: min(x["phastCons_0"], x["phastCons_1"]) > 0.5, axis=1)]

for _, row in df_logit.iterrows():
	logit_conserved.add(site(row["chr"], row["strand"], row["pos"]))

for _, row in df_phast.iterrows():
	phast_conserved.add(site(row["chr"], row["strand"], row["pos"]))

out = open(out_file, "w")
tr_out = open(tr_out_file, "w")

def read_introns(handle):
	trid = ""
	introns = []
	for record in getline(handle):
		now_trid = record.attr["transcript_id"]
		if now_trid != trid:
			if introns != []:
				yield (trid, introns)
				introns = []
			trid = now_trid
		introns.append(record)
	if introns != []:
		yield (trid, introns)

chess_transcript = set()

all_tissues = list(coverage.keys())
all_tissues.sort()
header = ["dataset", "gene_type", "inMANE", "chr", "strand", "start", "end", "donor_mane", "acceptor_mane", "donor_well_supported", "acceptor_well_supported"] + all_tissues
print(",".join(header), file=out)
tr_header = ["dataset", "transcript_id", "gene_type", "inMANE", "chr", "strand", "start", "end", "total_sites", "mane_sites", "well_supported_non_mane_sites", "well_supported"]
print(",".join(tr_header), file=tr_out)
mane_intron = set()
mane_transcript = set()
for db in db_list:
	title, introns_dir, gtf = db.split("!")
	if title[0] == '"':
		title = title[1:-1]
	handle = gzip.open(gtf, "rt")
	transcript_type = parse_type(handle)
	handle.close()
	transcript_coords = parse_transcript_coords(gzip.open(gtf, "rt"))
	for chr in os.listdir(introns_dir):
		if "chrY" in chr:
			continue
		seen = set()
		introns_handle = open(os.path.join(introns_dir, chr))
		for (trid, introns_list) in read_introns(introns_handle):
			mane_sites_count = 0
			total_sites_count = 0
			conserved_non_mane_sites_count = 0
			trid_type = get_type(introns_list[0], transcript_type)
			if not trid_type in all_types:
				continue
			tr_id = introns_list[0].attr["transcript_id"]
			now_coords = transcript_coords[trid]
			tr_signature = [str(x) for x in now_coords]
			for rec in introns_list:
				all_coverage = []
				donor_coords, acceptor_coords = donor_acceptor(rec.chr, rec.strand, rec.start, rec.end)
				intron_coords = "&".join((donor_coords, acceptor_coords))
				tr_signature.append(intron_coords)
				if title == "MANE":
					mane_intron.add(intron_coords)

				all_coverage = []

				acceptor_mane = acceptor_coords in mane_sites
				acceptor_logit_conserved = acceptor_coords in logit_conserved
				if acceptor_mane:
					mane_sites_count += 1
				elif acceptor_logit_conserved:
					conserved_non_mane_sites_count += 1

				donor_mane = donor_coords in mane_sites
				donor_logit_conserved = donor_coords in logit_conserved
				if donor_mane:
					mane_sites_count += 1
				elif donor_logit_conserved:
					conserved_non_mane_sites_count += 1

				total_sites_count += 2

				cons = [donor_mane, acceptor_mane, donor_logit_conserved, acceptor_logit_conserved]
				cons = [str(int(x)) for x in cons]
				for tissue, now_coverage in coverage.items():
					if donor_coords in now_coverage and acceptor_coords in now_coverage[donor_coords]:
						cov = now_coverage[donor_coords][acceptor_coords]
						all_coverage.append(str(cov))
					else:
						all_coverage.append("")
				in_mane = "1" if intron_coords in mane_intron else "0"
				data = [title, trid_type, in_mane, rec.chr, rec.strand, str(rec.start), str(rec.end)] + cons + all_coverage
				if not intron_coords in seen:
					print(",".join(data), file=out)
					seen.add(intron_coords)

			tr_signature = "!".join(tr_signature)
			if title == "MANE":
				mane_transcript.add(tr_signature)
				print(tr_signature)
			cons_status = "1" if total_sites_count == mane_sites_count + conserved_non_mane_sites_count else "0"
			now_mane_transcript = "1" if tr_signature in mane_transcript else "0"
			print(",".join((title, trid, trid_type, now_mane_transcript, now_coords[0], now_coords[1], str(now_coords[2]), str(now_coords[3]),  str(total_sites_count), str(mane_sites_count), str(conserved_non_mane_sites_count), cons_status)), file=tr_out)


