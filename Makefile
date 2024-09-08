#################################################################################
# GLOBALS                                                                       #
#################################################################################

PROJECT_NAME = splice-sites-conservation
PYTHON_INTERPRETER = python
CONFIGURATION = debug


ifeq ($(CONFIGURATION),debug)
chromosomes := chr21 chr22
rand_limit := 20000
ucsc_genomes := panTro6 panPan3
zoo_genomes := HLmacFus1 HLallNig1
ncbi_genomes := HLhylMol2 HLsemEnt1
else
chromosomes =
endif

gencode_version := 45
refseq_version := 110
mane_version := 1.3
chess_version := 3.1.0

#################################################################################
# RAW FILES                                                                     #
#################################################################################
human_version := GCF_000001405.40_GRCh38.p14
human_genome := data/raw/human/$(human_version)_genomic.fna.gz
human_stats := data/raw/human/$(human_version)_assembly_report.txt
phast_wig := $(foreach chr,$(chromosomes),data/raw/phast/$(chr).phastCons470way.wigFix.gz)
maf_files := $(foreach chr,$(chromosomes),data/raw/maf/$(chr).maf)
gnomad_files := $(foreach chr,$(chromosomes),data/raw/gnomad/gnomad.genomes.v4.0.sites.$(chr).vcf.bgz)
ucsc_genomes_files := $(foreach genome,$(ucsc_genomes),data/raw/genomes/$(genome).fa.gz)
zoo_genomes_files := $(foreach genome,$(zoo_genomes),data/raw/genomes/$(genome).fa.gz)
ncbi_genomes_files := $(foreach genome,$(ncbi_genomes),data/raw/genomes/$(genome).fa.gz)
all_genomes := $(ucsc_genomes) $(zoo_genomes) $(ncbi_genomes)

gencode_gtf := data/raw/annotation/gencode.v$(gencode_version).annotation.gtf.gz
refseq_filename := $(human_version)_genomic.gtf.gz
refseq_gtf := data/raw/annotation/$(refseq_filename)
mane_gtf := data/raw/annotation/MANE.GRCh38.v$(mane_version).ensembl_genomic.gtf.gz
chess_gtf := data/raw/annotation/chess$(chess_version).GRCh38.gtf.gz
random_gtf := data/raw/annotation/random.gtf.gz

#################################################################################
# INTERIM FILES                                                                 #
#################################################################################

gnomad_tracks := $(foreach chr,$(chromosomes),data/interim/gnomad/$(chr).csv)

gencode_dir := data/interim/gencode/$(gencode_version)
gencode_introns := $(foreach chr,$(chromosomes),$(gencode_dir)/introns_all/$(chr))
gencode_query := $(foreach chr,$(chromosomes),$(gencode_dir)/query_all/$(chr))

refseq_dir := data/interim/refseq/$(refseq_version)
refseq_introns := $(foreach chr,$(chromosomes),$(refseq_dir)/introns_all/$(chr))
refseq_query := $(foreach chr,$(chromosomes),$(refseq_dir)/query_all/$(chr))

mane_dir := data/interim/mane/$(mane_version)
mane_introns := $(foreach chr,$(chromosomes),$(mane_dir)/introns_all/$(chr))
mane_query := $(foreach chr,$(chromosomes),$(mane_dir)/query_all/$(chr))

chess_dir := data/interim/chess/$(chess_version)
chess_introns := $(foreach chr,$(chromosomes),$(chess_dir)/introns_all/$(chr))
chess_query := $(foreach chr,$(chromosomes),$(mane_dir)/query_all/$(chr))

random_dir := data/interim/random/1.0
random_introns := $(foreach chr,$(chromosomes),$(random_dir)/introns_all/$(chr))
random_query := $(foreach chr,$(chromosomes),$(random_dir)/query_all/$(chr))

realignment_dir := data/interim/realignment
unique_exons := $(foreach chr,$(chromosomes),$(realignment_dir)/unique_exons/$(chr))
mapped_exons := $(foreach chr,$(chromosomes),$(realignment_dir)/mapped_exons/$(chr))
alignment_jobs := $(foreach chr,$(all_genomes),$(realignment_dir)/jobs/$(genome))
alignment_result := $(foreach chr,$(all_genomes),$(realignment_dir)/alignment/$(genome))

#################################################################################
# COMMANDS                                                                      #
#################################################################################

#all: requirements $(maf_files) $(gnomad_tracks) $(ucsc_genomes_files) $(zoo_genomes_files) $(ncbi_genomes_files) $(gencode_query) $(refseq_query) $(mane_query) $(chess_query) $(random_query) data/interim/clinvar/clinvar.csv

all: $(alignemnt_jobs)

clean:
	rm -rf data/raw/annotation/*
	rm -rf $(gencode_dir)
	rm -rf $(chess_dir)
	rm -rf $(mane_dir)
	rm -rf $(refseq_dir)
	rm -rf $(random_dir)

## Install requirements

requirements:

## Realign exons to other genomes

$(alignemnt_jobs): $(mapped_exons)
	$(PYTHON_INTERPRETER) src/data/map_exons.py

## Map exons to other genomes

$(mapped_exons): $(unique_exons) $(maf_files)
	$(PYTHON_INTERPRETER) src/data/map_exons.py data/raw/maf/$(@F).maf $(realignment_dir)/unique_exons/$(@F) $(realignment_dir)/mapped_exons/$(@F)

## Generate unique exons from all datasets:
$(unique_exons): $(random_gtf) $(gencode_gtf) $(refseq_gtf) $(chess_gtf)
	$(PYTHON_INTERPRETER) src/data/unique_exons.py $(gencode_dir)/exons/$(@F) $(refseq_dir)/exons/$(@F) $(chess_dir)/exons/$(@F) $(mane_dir)/exons/$(@F) $(random_dir)/exons/$(@F) > $(realignment_dir)/unique_exons/$(@F)

## Get GRCh38:

$(human_stats):
	wget --directory-prefix=data/raw/human https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/$(human_version)/$(human_version)_assembly_report.txt

$(human_genome):
	wget --directory-prefix=data/raw/human https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/$(human_version)/$(human_version)_genomic.fna.gz

## Query random dataset

$(random_query): $(random_introns) $(maf_files)
	$(PYTHON_INTERPRETER) src/data/query_chr.py data/interim/random/1.0 data/raw/maf/$(@F).maf $(@F)
	$(PYTHON_INTERPRETER) src/data/collect_queries.py data/interim/random/1.0/ $(@F)

$(random_introns): $(random_gtf)
	src/data/make_introns.sh data/interim/random/1.0 $(@F)

## Generate random dataset

$(random_gtf): $(human_genome) $(human_stats) $(mane_introns)
	$(PYTHON_INTERPRETER) src/data/generate_random.py data/interim/mane/$(mane_version)/introns $(human_stats) $(human_genome) $(rand_limit) > $(basename $(random_gtf))
	gzip $(basename $(random_gtf))
	$(PYTHON_INTERPRETER) src/data/split_gtf.py $(random_gtf) data/interim/random/1.0/ "$(chromosomes)"

## Query GENCODE

$(gencode_query): $(gencode_introns) $(maf_files)
	$(PYTHON_INTERPRETER) src/data/query_chr.py data/interim/gencode/$(gencode_version) data/raw/maf/$(@F).maf $(@F)
	$(PYTHON_INTERPRETER) src/data/collect_queries.py data/interim/gencode/$(gencode_version)/ $(@F)

$(gencode_introns): $(gencode_gtf)
	src/data/make_introns.sh data/interim/gencode/$(gencode_version) $(@F)

## Get GENCODE

$(gencode_gtf):
	wget --directory-prefix=data/raw/annotation https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$(gencode_version)/gencode.v$(gencode_version).annotation.gtf.gz
	$(PYTHON_INTERPRETER) src/data/split_gtf.py $(gencode_gtf) data/interim/gencode/$(gencode_version)/ "$(chromosomes)"

## Query RefSeq
$(refseq_query): $(refseq_introns) $(maf_files)
	$(PYTHON_INTERPRETER) src/data/query_chr.py data/interim/refseq/$(refseq_version) data/raw/maf/$(@F).maf $(@F)
	$(PYTHON_INTERPRETER) src/data/collect_queries.py data/interim/refseq/$(refseq_version)/ $(@F)

$(refseq_introns): $(refseq_gtf)
	src/data/make_introns.sh data/interim/refseq/$(refseq_version) $(@F)

## Get RefSeq:

$(refseq_gtf): $(human_stats)
	wget --directory-prefix=data/raw/annotation https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/110/$(human_version)/$(refseq_filename)
	$(PYTHON_INTERPRETER) src/data/gtf_replace.py $(human_stats) $(refseq_gtf) data/raw/annotation/temp.gtf
	rm $(refseq_gtf)
	gzip data/raw/annotation/temp.gtf
	mv data/raw/annotation/temp.gtf.gz data/raw/annotation/$(refseq_filename)
	$(PYTHON_INTERPRETER) src/data/split_gtf.py $(refseq_gtf) data/interim/refseq/$(refseq_version)/ "$(chromosomes)"


## Query MANE
$(mane_query): $(mane_introns) $(maf_files)
	$(PYTHON_INTERPRETER) src/data/query_chr.py data/interim/mane/$(mane_version) data/raw/maf/$(@F).maf $(@F)
	$(PYTHON_INTERPRETER) src/data/collect_queries.py data/interim/mane/$(mane_version)/ $(@F)

$(mane_introns): $(mane_gtf)
	src/data/make_introns.sh data/interim/mane/$(mane_version) $(@F)

## Get MANE:

$(mane_gtf):
	wget --directory-prefix=data/raw/annotation https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_$(mane_version)/MANE.GRCh38.v$(mane_version).ensembl_genomic.gtf.gz
	$(PYTHON_INTERPRETER) src/data/split_gtf.py $(mane_gtf) data/interim/mane/$(mane_version)/ "$(chromosomes)"


## Query CHESS
$(chess_query): $(chess_introns) $(maf_files)
	$(PYTHON_INTERPRETER) src/data/query_chr.py data/interim/chess/$(chess_version) data/raw/maf/$(@F).maf $(@F)
	$(PYTHON_INTERPRETER) src/data/collect_queries.py data/interim/chess/$(chess_version)/ $(@F)

$(chess_introns): $(chess_gtf)
	src/data/make_introns.sh data/interim/chess/$(chess_version) $(@F)

## Get CHESS:

$(chess_gtf):
	wget --directory-prefix=data/raw/annotation https://github.com/chess-genome/chess/releases/download/v.$(chess_version)/chess$(chess_version).GRCh38.gtf.gz
	$(PYTHON_INTERPRETER) src/data/split_gtf.py $(chess_gtf) data/interim/chess/$(chess_version)/ "$(chromosomes)"


## Get and parse the clinvar data

data/interim/clinvar/clinvar.csv: data/raw/clinvar/clinvar.vcf.gz $(human_stats)
	$(PYTHON_INTERPRETER) src/data/parse_clinvar.py src/data/clinvar_pathogenic.txt $(human_stats) data/raw/clinvar/clinvar.vcf.gz data/interim/clinvar/clinvar.csv

data/raw/clinvar/clinvar.vcf.gz:
	wget --directory-prefix=data/raw/clinvar https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz

## Get the alignmnet

$(maf_files):
	wget --directory-prefix=data/raw/maf https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz470way/maf/$(@F)

## Get and parse the gnomAD data

$(gnomad_tracks): $(gnomad_files)
	$(PYTHON_INTERPRETER) src/data/parse_gnomad.py data/raw/gnomad $@

$(gnomad_files):
	wget --directory-prefix=data/raw/gnomad https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.0/vcf/genomes/$(@F)

## Get the phastCons data

$(phast_wig):
	wget --directory-prefix=data/raw/phast https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons470way/hg38.470way.phastCons/$(@F)

## Get the genomes from the alignment hosted by UCSC

$(ucsc_genomes_files):
	wget --directory-prefix=data/raw/genomes https://hgdownload.soe.ucsc.edu/goldenPath/$(basename $(basename $(@F)))/bigZips/$(@F)

## Get the genomes from the alignment hosted by DNA Zoo

$(zoo_genomes_files):
	src/data/download_zoo.sh src/data/dna_zoo.txt $(@F) data/raw/genomes

## Get the genomes from the alignment hosted by NCBI

$(ncbi_genomes_files):
	src/data/download_ncbi.sh src/data/ncbi.txt $(@F) data/raw/genomes data/tmp
