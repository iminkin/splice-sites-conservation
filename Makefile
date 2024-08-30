#################################################################################
# GLOBALS                                                                       #
#################################################################################


PROJECT_NAME = splice-sites-conservation
PYTHON_INTERPRETER = python
CONFIGURATION = debug

ifeq ($(CONFIGURATION),debug)
chromosomes := chr21 chr22
ucsc_genomes := panTro6 panPan3
zoo_genomes := HLmacFus1 HLallNig1
ncbi_genomes := HLhylMol2
#HLsemEnt1
else
chromosomes =
endif

maf_files := $(foreach chr,$(chromosomes),data/raw/maf/$(chr).maf)
gnomad_files := $(foreach chr,$(chromosomes),data/raw/gnomad/gnomad.genomes.v4.0.sites.$(chr).vcf.bgz)
ucsc_genomes_files := $(foreach genome,$(ucsc_genomes),data/raw/genomes/$(genome).fa.gz)
zoo_genomes_files = $(foreach genome,$(zoo_genomes),data/raw/genomes/$(genome).fa.gz)
ncbi_genomes_files = $(foreach genome,$(ncbi_genomes),data/raw/genomes/$(genome).fa.gz)

#################################################################################
# COMMANDS                                                                      #
#################################################################################

all: requirements $(maf_files) $(gnomad_files) $(ucsc_genomes_files) $(zoo_genomes_files) $(ncbi_genomes_files) data/raw/clinvar/clinvar.vcf.gz

## Install requirements

requirements:

## Make Dataset

data/raw/clinvar/clinvar.vcf.gz:
	wget --directory-prefix=data/raw/clinvar https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz

$(maf_files):
	wget --directory-prefix=data/raw/maf https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz470way/maf/$(@F)

$(gnomad_files):
	wget --directory-prefix=data/raw/gnomad https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.0/vcf/genomes/$(@F)

$(ucsc_genomes_files):
	wget --directory-prefix=data/raw/genomes https://hgdownload.soe.ucsc.edu/goldenPath/$(basename $(basename $(@F)))/bigZips/$(@F)

$(zoo_genomes_files):
	src/data/download_zoo.sh src/data/dna_zoo.txt $(@F) data/raw/genomes

$(ncbi_genomes_files):
	src/data/download_ncbi.sh src/data/ncbi.txt $(@F) data/raw/genomes data/tmp
