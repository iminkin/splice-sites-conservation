import sys
from collections import namedtuple

gtf_record = namedtuple("gtf_record", ["type", "src", "chr", "strand", "start", "end", "attr", "attr_orig", "line"])

def getline(handle):
	for line in handle:
		orig = line.strip()
		if line[0] == "#":
			continue
		attr = dict()
		line = line.strip().split("\t")
		tag_array = line[-1].split(";")
		for tag in tag_array:
			tag = tag.strip().split(" ")
			if len(tag) == 2:
				attr[tag[0]] = tag[1][1:-1]
		yield gtf_record(line[2], line[1], line[0], line[6], int(line[3]), int(line[4]), attr, line[-1], orig)

def print_line(line, handle):
	print("\t".join((line.type, line.src, line.chr, line.strand, line.start, line.end, line.attr_orig)), file=handle)

def getline_filter(handle, filter):
	for line in getline(handle):
		if line.type == filter:
			yield line
