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

source_dir = sys.argv[1]
target = sys.argv[2]

chr = target.split("/")[-1].split(".")[0]
src_vcf = "gnomad.genomes.v4.0.sites.{}.vcf.bgz".format(chr)
src_vcf = os.path.join(source_dir, src_vcf)

out_handle = open(target, "w")
print("pos,AF", file=out_handle)
with gzip.open(src_vcf, "rt") as handle:
	for line in handle:
		line = parse_vcf_line(line.strip())
		if line != {} and int(line["attr"]["nhomalt"]) > 0:
			print("{},{}".format(str(line["pos"]), line["attr"]["AF"]), file=out_handle)

