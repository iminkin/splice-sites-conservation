import os
import sys
import gzip

def parse_vcf_line(line):
	ret = dict()
	if line[0] != "#":
		line = line.split()
		ret["chr"] = line[0]
		ret["pos"] = int(line[1])
		ret["attr"] = dict()
		for item in line[-1].split(";"):
			item = item.split("=")
			if len(item) == 2:
				ret["attr"][item[0]] = item[1]
			else:
				ret["attr"][item[0]] = None
	return ret

pathogenic = set()
for sig in open(sys.argv[1]):
	pathogenic.add(sig.strip())

print(pathogenic)

chrom = dict()
for line in open(sys.argv[2]):
        line = line.strip().split()
        chrom[line[0]] = line[-1]


src_vcf = sys.argv[3]
out = sys.argv[4]

out_handle = open(out, "w")
significance_field = "CLNSIG"
print("chr,pos,type", file=out_handle)

with gzip.open(src_vcf, "rt") as handle:
	for line in handle:
		line = parse_vcf_line(line.strip())
		if "attr" in line and significance_field in line["attr"]:
			now_significance = line["attr"][significance_field].lower()
			if now_significance in pathogenic:
				print("{},{},{}".format(line["chr"], str(line["pos"]), now_significance), file=out_handle)

