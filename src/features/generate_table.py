import os
import sys
import math
import gzip

from Bio import SeqIO
from Bio.Seq import Seq
sys.path.append("src/lib")
from gtf_parse import getline

nn = 62
main_pos = nn / 2 - 1
half_motif = int((nn - 2) / 2)

gen_header_idx = ["cons_" + str(idx - half_motif) for idx in range(nn)] + ["snp_" + str(idx - half_motif) for idx in range(nn)]
gen_header = ",".join(gen_header_idx)

suffix = ["a", "d"]
alpha = "ACGT"

def donor_acceptor(chr, strand, start, end):
	if strand == '+':
		donor_coords = "&".join((chr, strand, str(start)))
		acceptor_coords = "&".join((chr, strand, str(end - 1)))
	else:
		donor_coords = "&".join((chr, strand, str(end)))
		acceptor_coords = "&".join((chr, strand, str(start + 1)))
	return (donor_coords, acceptor_coords)

def parse_variants(path, chr):
	variants = dict()
	file = os.path.join(path, chr)
	if os.path.isfile(file):
		line = open(file).readline()
		line = line.strip().split(",")
		variants[chr] = dict()
		for idx, af in enumerate(line):
			af = float(af)
			if af > 0:
				variants[chr][idx] = af
	return variants

def get_type(dataset, gtf, chr):
	gene_type = dict()
	trid_to_geneid = dict()
	transcript_type = dict()
	transcripts_per_gene = dict()
	if dataset == "RefSeq":
		type_attr = "gene_biotype"
	else:
		type_attr = "gene_type"

	for rec in getline(gzip.open(sys.argv[1], "rt")):
		if rec.chr != chr:
			continue

		if rec.type == "gene":
			gene_id = rec.attr["gene_id"]
			gene_type[gene_id] = rec.attr[type_attr]
#			print(rec.attr, file=sys.stderr)

		if rec.type == "transcript":
			transcript_id = rec.attr["transcript_id"]
			if "gene_id" in rec.attr:
				gene_id = rec.attr["gene_id"]
				if gene_id in gene_type:
					transcript_type[transcript_id] = gene_type[gene_id]
				else:
					transcript_type[transcript_id] = rec.attr[type_attr]

				if not gene_id in transcripts_per_gene:
					transcripts_per_gene[gene_id] = 0
				transcripts_per_gene[gene_id] += 1
				trid_to_geneid[transcript_id] = gene_id
			elif type_attr in rec.attr:
				transcript_type[transcript_id] = rec.attr[type_attr]

	return (gene_type, trid_to_geneid, transcript_type, transcripts_per_gene)

def pack_header(trid, intron_idx, suffix):
	return "$".join((trid, str(intron_idx), suffix))

def unpack_header(header):
	header = header.split("$")
	return (header[0], int(header[1]), header[2])

def get_coords_set(path):
	coords_set = {"a" : set(), "d" : set()}
	for root, dirs, files in os.walk(path):
		for trid in files:
			for rec in getline(open(os.path.join(root, trid))):
				donor_coords, acceptor_coords = donor_acceptor(rec.chr, rec.strand, rec.start, rec.end)
				coords_set["d"].add(donor_coords)
				coords_set["a"].add(acceptor_coords)
	return coords_set


def get_coords(introns_dir, chr):
	site_coords = dict()
	coords_use_rate = {"a" : dict(), "d" : dict()}
	root = os.path.join(introns_dir, chr)
	for trid in os.listdir(root):
		intron_idx = 0
		h = open(os.path.join(root, trid))
		for rec in getline(h):
			donor_coords, acceptor_coords = donor_acceptor(rec.chr, rec.strand, rec.start, rec.end)
			if not donor_coords in coords_use_rate["d"]:
				coords_use_rate["d"][donor_coords] = 0
			coords_use_rate["d"][donor_coords] += 1

			if not acceptor_coords in coords_use_rate["a"]:
				coords_use_rate["a"][acceptor_coords] = 0
			coords_use_rate["a"][acceptor_coords] += 1

			site_coords[pack_header(trid, str(intron_idx), "d")] = donor_coords
			site_coords[pack_header(trid, str(intron_idx), "a")] = acceptor_coords
			intron_idx += 1
	return (site_coords, coords_use_rate)

def get_extra_cons(extra_dir, chr_dir):
	cons = dict()
	for genome in os.listdir(extra_dir):
		cons[genome] = dict()
		for line in open(os.path.join(extra_dir, genome)):
			line = line.strip().split()
			if line[0] == chr_dir:
				cons[genome][line[0]] = set([int(p) for p in line[1:]])
	return cons

def parse_phast(phast_base, chr):
	phast = dict()
	phast[chr] = dict()
	file_path = os.path.join(phast_base, chr + ".wigFix")
	if os.path.isfile(file_path):
		for line in gzip.open(file_path, "rt"):
			line = line.strip().split()
			if len(line) > 1:
				start = int(line[2].split('=')[1])
			else:
				phast[chr][start] = float(line[0])
				start += 1
	return phast

def parse_clinvar(clinvar):
	ret = dict()
	handle = open(clinvar)
	handle.readline()
	for line in handle:
		line = line.strip().split(",")
		chr, pos = "chr" + line[0], int(line[1])
		if not chr in ret:
			ret[chr] = set()
		ret[chr].add(pos)
	return ret


#hg38 = dict()
#for record in SeqIO.parse("/home/iminkin2/projects3/splice-sites-paper-final/data/hg38/hg38.fa", "fasta"):
#	hg38[record.id] = Seq(record.seq.upper())

gtf = sys.argv[1]
dataset = sys.argv[2]
introns_dir = sys.argv[3]
query_dir = sys.argv[4]
mane_introns = sys.argv[5]
snp_base = sys.argv[6]
clinvar_file = sys.argv[7]
extra_cons_dir = sys.argv[8]
phast_dir = sys.argv[9]
all_genomes = set(sys.argv[10].split())
clinvar = parse_clinvar(clinvar_file)

limit = 180000 if dataset == "Random" else sys.maxsize

count = {"a" : 0, "d" : 0}
all_header = "dataset,transcript_id,intron_index,site_type,gene_type,inMANE,chr,strand,pos,cons_GTAG,motif," + gen_header + ",reuse,phastCons_0,phastCons_1,clinvar_0,clinvar_1"
print(all_header)
printed = set()

all_transcripts = dict()
for chr_dir in os.listdir(query_dir):
	extra_cons = get_extra_cons(extra_cons_dir, chr_dir)
	gene_type, trid_to_geneid, transcript_type, transcripts_per_gene = get_type(dataset, gtf, chr_dir)
	site_coords, coords_use_rate = get_coords(introns_dir, chr_dir)

	site_seen = {"a" : dict(), "d" : dict()}
	chr_path = os.path.join(query_dir, chr_dir)
	phast = parse_phast(phast_dir, chr_dir)
	variants = parse_variants(snp_base, chr_dir)
	mane_coords = get_coords_set(os.path.join(mane_introns, chr_dir))
	for site_id in os.listdir(chr_path):
		rec = dict()
		score = dict()
		trid, site_idx, suffix = unpack_header(site_id)
		type = transcript_type[trid]
		if type != "lncRNA" and type != "protein_coding":
			continue
#		print(trid, file=sys.stderr)
		coords = site_coords[site_id]
		if coords in site_seen[suffix]:
			site_signature, site_usage = site_seen[suffix][coords]
			site_seen[suffix][coords] = (site_signature, site_usage + 1)
			continue
		else:
			site_seen[suffix][coords] = (site_id, 1)

		gtat_conserved = dict()
		for genome in all_genomes:
			gtat_conserved[genome] = 1

		chr, strand, now_pos = coords.split("&")
		now_pos = int(now_pos)

		for record in SeqIO.parse(os.path.join(chr_path, site_id), "fasta"):
			genome = record.id.split(".")[0]
			seq = record.seq.upper()
			rec[genome] = seq

		idx = 0
		var_freq = [0] * nn
		cons_count = [0] * nn
		if strand == '+':
			genome_pos = now_pos - half_motif
			check_start = genome_pos - 1
#			check = hg38[chr][check_start:check_start + nn]
			inc = +1
		else:
			genome_pos = now_pos + half_motif
			check_start = genome_pos - 1 - nn + 1
#			check = hg38[chr][check_start:check_start + nn].reverse_complement()
			inc = -1

		site = ''.join((c for c in rec["hg38"] if c != "-"))
		now_motif = ''
		for pos, c in enumerate(rec["hg38"]):
			if c != '-':
#				original_c = hg38[chr][genome_pos - 1]
#				if strand == "-":
#					original_c = Seq(original_c).complement()

#				if c != original_c:
#					print(site, now_pos, genome_pos, strand, site_id, half_motif, file=sys.stderr)
#					print(check, file=sys.stderr)

#				assert (c == original_c) or c == "N" or original_c == "N"
				if idx == half_motif or idx == half_motif + 1:
					now_motif = now_motif + c

				for genome in all_genomes:
					if genome in rec:
						seq = rec[genome]
						if seq[pos] == c:
							cons_count[idx] += 1
						elif idx == half_motif or idx == half_motif + 1:
							gtat_conserved[genome] = 0
					else:
						if genome in extra_cons and chr in extra_cons[genome] and pos in extra_cons[genome][chr]:
							cons_count[idx] += 1
						elif idx == half_motif or idx == half_motif + 1:
							gtat_conserved[genome]

				if chr in variants and (genome_pos - 1) in variants[chr]:
					var_freq[idx] = variants[chr][genome_pos - 1]

				idx += 1
				genome_pos += inc

		if count[suffix] < limit:
			cons_array = [str(c) for c in cons_count]
			freq_array = [str(f) for f in var_freq]
			gtat_cons_count = str(list(gtat_conserved.values()).count(1))
			if trid in trid_to_geneid and trid_to_geneid[trid] in transcripts_per_gene:
				gene_id = trid_to_geneid[trid]
				reuse = str(coords_use_rate[suffix][coords])
			else:
#				print(trid, trid in trid_to_geneid, trid in transcripts_per_gene, file=sys.stderr)
				reuse = ""

			mane = "1" if coords in mane_coords[suffix] else "0"
			if mane == "1" and dataset == "Random":
				continue

			phast_pos_0 = now_pos
			if strand == "+":
				phast_pos_1 = phast_pos_0 + 1
			else:
				phast_pos_1 = phast_pos_0 - 1

			phast_chr = phast[chr] if chr in phast else dict()
			phast_0 = phast_chr[phast_pos_0] if phast_pos_0 in phast_chr else math.nan
			phast_1 = phast_chr[phast_pos_1] if phast_pos_1 in phast_chr else math.nan

			assert dataset != "MANE" or (dataset == "MANE" and mane == "1")
			clinvar_chr = clinvar[chr] if chr in clinvar else dict()
			clinvar_0 = "1" if phast_pos_0 in clinvar_chr else "0"
			clinvar_1 = "1" if phast_pos_1 in clinvar_chr else "0"
			val = [dataset, trid, str(site_idx), suffix, type, mane, chr, strand, str(now_pos), gtat_cons_count, now_motif] + cons_array + freq_array + [reuse, str(phast_0), str(phast_1), clinvar_0, clinvar_1]
			row = ",".join(val)
			print(row)

			if row.count(",") != all_header.count(","):
				print(all_header.split(","))
				print(row.split(","))
				print(row.count(","), all_header.count(","), len(cons_array), len(freq_array))

			assert row.count(",") == all_header.count(",")
			count[suffix] += 1

