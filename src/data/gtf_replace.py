import gzip
import sys

cnt = dict()
mol = dict()

for line in open(sys.argv[1]):
	line = line.strip()
	if line[0] != '#':
		line = line.split('\t')
		if line[-1] != "na":
			mol[line[6]] = line[-1]

out_handle = open(sys.argv[3], "w")
for line in gzip.open(sys.argv[2], 'rt'):
	line = line.strip()
	orig = line
	if line[0] != '#':
		line = line.split('\t')
		if len(line) > 2 and (line[2] == "mRNA" or line[2] == "lnc_RNA"):
			line[2] = "transcript"
		if line[0] in mol:
			line[0] = mol[line[0]]
			line = '\t'.join(line)
			print(line, file=out_handle)
	else:
		print(line, file=out_handle)

