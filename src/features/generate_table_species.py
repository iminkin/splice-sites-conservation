import os
import sys
import math
import gzip
import bisect

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
sys.path.append("src/lib")
from gtf_parse import getline

nn = 62
main_pos = nn / 2 - 1
half_motif = int((nn - 2) / 2)

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

def get_type(dataset, gtf):
	gene_type = dict()
	trid_to_geneid = dict()
	transcript_type = dict()
	transcripts_per_gene = dict()
	if dataset == "RefSeq":
		type_attr = "gene_biotype"
	else:
		type_attr = "gene_type"

	for rec in getline(gzip.open(sys.argv[1], "rt")):
#		if rec.chr != chr:
#			continue

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
	for rec in getline(open(path)):
		donor_coords, acceptor_coords = donor_acceptor(rec.chr, rec.strand, rec.start, rec.end)
		coords_set["d"].add(donor_coords)
		coords_set["a"].add(acceptor_coords)
	return coords_set

def get_pos_array(path):
	ret = []
	for rec in getline(open(path)):
		for i in range(rec.start, rec.end + 1):
			ret.append(i)
	ret.sort()
	return ret

def get_coords(introns_file):
	site_coords = dict()
	coords_use_rate = {"a" : dict(), "d" : dict()}
	prev_trid = ""
	intron_idx = 0
	for rec in getline(open(introns_file)):
		trid = rec.attr["transcript_id"]
		if trid != prev_trid:
			intron_idx = 0
			prev_trid = trid
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

def get_extra_cons(extra_dir, all_genomes, chr):
	cons = dict()
	for genome in all_genomes:
		cons[genome] = dict()
		for line in open(os.path.join(extra_dir, genome + "!" + chr)):
			line = line.strip().split()
			if line[0] == chr:
				cons[genome][line[0]] = set([int(p) for p in line[1:]])
	return cons

def bin_search(a, x):
	i = bisect.bisect_left(a, x)
	if i != len(a) and a[i] == x:
        	return True
	return False

def parse_minor_introns(csv):
	ret = dict()
	handle = open(csv)
	handle.readline()
	for line in handle:
		line = line.strip().split(",")
		coord = line[3].split("_")[2:]
		chr = "chr" + coord[0]
		start, end = int(coord[1]), int(coord[2])
		strand = "+" if coord[3] == "F" else "-"
		donor, acceptor = donor_acceptor(chr, strand, start, end)
		if not chr in ret:
			ret[chr] = {"d": set(), "a" : set()}
		if line[13] == "minor":
			ret[chr]["d"].add(donor)
			ret[chr]["a"].add(acceptor)
	return ret

def site_in_exon(mane_exon_coords, pos, strand, site_type):
	inc = +1 if strand == "+" else -1
	return bin_search(mane_exon_coords, pos) and bin_search(mane_exon_coords, pos + inc)

mol = dict()
for line in open(sys.argv[15]):
	line = line.strip()
	if line[0] != '#':
		line = line.split('\t')
		if line[-1] != "na":
			mol[line[6]] = line[-1]

hg38 = dict()
#hg38_handle = gzip.open(sys.argv[14], mode='rt')
#for record in SeqIO.parse(hg38_handle, "fasta"):
#	record_id = mol[record.id]
#	hg38[record_id] = record.seq.upper()

gtf = sys.argv[1]
dataset = sys.argv[2]
introns_dir = sys.argv[3]
query_dir = sys.argv[4]
minor_introns = sys.argv[5]
mane_exons = sys.argv[6]
mane_introns = sys.argv[7]
snp_base = sys.argv[8]
clinvar_file = sys.argv[9]
extra_cons_dir = sys.argv[10]
phast_dir = sys.argv[11]
all_genomes = set(sys.argv[12].split())
chromosomes = sys.argv[13].split()


limit = 180000 if dataset == "Random" else sys.maxsize

count = {"a" : 0, "d" : 0}
all_header = "dataset,transcript_id,intron_index,site_type,gene_type,inMANE,chr,strand,pos,cons_GTAG,motif,inside_mane_exon,minor_intron,cons_array"
print(all_header)
printed = set()

def parse_query_batch(handle):
	prev = ""
	batch = []
	for record in SeqIO.parse(handle, "fasta"):
		head = record.id.split("#")
		if head[0] != prev and batch != []:
			ret = prev.split("$")
			yield(ret[0], ret[1], ret[2], batch)
			batch = []
		prev = head[0]
		genome = head[1].split(".")[0]
		batch.append((genome, record.seq.upper()))
	ret = prev.split("$")
	if len(ret) == 3:
		yield(ret[0], ret[1], ret[2], batch)

minor_introns_set = parse_minor_introns(minor_introns)
gene_type, trid_to_geneid, transcript_type, transcripts_per_gene = get_type(dataset, gtf)
all_transcripts = dict()
for chr in chromosomes:
	query_path = os.path.join(query_dir, chr)
	if not os.path.isfile(query_path):
		continue
	init = False
	for (trid, site_idx, suffix, seq_batch) in parse_query_batch(open(query_path)):
		if not init:
			extra_cons = get_extra_cons(extra_cons_dir, all_genomes, chr)
			mane_coords = get_coords_set(os.path.join(mane_introns, chr))
			mane_exons_pos = get_pos_array(os.path.join(mane_exons, chr))
			site_coords, coords_use_rate = get_coords(os.path.join(introns_dir, chr))
			site_seen = {"a" : dict(), "d" : dict()}
			init = True

		rec = dict()
		score = dict()
		type = transcript_type[trid]
		if type != "lncRNA" and type != "protein_coding":
			continue

		site_id = pack_header(trid, site_idx, suffix)
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

		for (genome, seq) in seq_batch:
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
#
#				if c != original_c:
#					print(site, now_pos, genome_pos, strand, site_id, half_motif, file=sys.stderr)
#					print(check, file=sys.stderr)

#				assert (c == original_c) or c == "N" or original_c == "N"
				if idx == half_motif or idx == half_motif + 1:
					now_motif = now_motif + c

				for genome in all_genomes:
					if genome in rec:
						seq = rec[genome]
						if pos >= len(seq):
							print(genome)
						if seq[pos] == c:
							cons_count[idx] += 1
						elif idx == half_motif or idx == half_motif + 1:
							gtat_conserved[genome] = 0
					else:
						if genome in extra_cons and chr in extra_cons[genome] and pos in extra_cons[genome][chr]:
							cons_count[idx] += 1
						elif idx == half_motif or idx == half_motif + 1:
							gtat_conserved[genome] = 0

				idx += 1
				genome_pos += inc

		cons_array_species = ";".join([genome for genome in all_genomes if gtat_conserved[genome] == 1])
		if count[suffix] < limit:
			gtat_cons_count = str(list(gtat_conserved.values()).count(1))
			mane = "1" if coords in mane_coords[suffix] else "0"
			if mane == "1" and dataset == "Random":
				continue

			assert dataset != "MANE" or (dataset == "MANE" and mane == "1")
			inside_mane_exon = "1" if site_in_exon(mane_exons_pos, now_pos, strand, suffix) else "0"
			is_minor_intron = "1" if chr in minor_introns_set and coords in minor_introns_set[chr][suffix] else "0"
			val = [dataset, trid, str(site_idx), suffix, type, mane, chr, strand, str(now_pos), gtat_cons_count, now_motif, str(inside_mane_exon), is_minor_intron, cons_array_species]
			row = ",".join(val)
			print(row)

			if row.count(",") != all_header.count(","):
				print(all_header.split(","))
				print(row.split(","))
				print(row.count(","), all_header.count(","), len(cons_array), len(freq_array))

			assert row.count(",") == all_header.count(",")
			count[suffix] += 1

