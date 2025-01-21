Splice site conservation analysis pipeline
==========================================

Code necessary to reproduce the analyses from the paper 
"Conservation assessment of human splice site annotation based on a 470-genome alignment"
by Ilia Minkin and Steven L. Salzberg ([biorxiv link](https://www.biorxiv.org/content/10.1101/2023.12.01.569581v2)).

The documentation is still in progress, here is a brief description of the resulting data files:

* data/processed/splice_sites.csv.gz: the input file for the regression models; summarizes all
 unique splice from all the datasets.

* data/processed/model_out.csv.gz: the CSV table containing all unique splice sites from
 each annotation we analyzed, along with their properties and support status. The description
 is in a subsection below.

* data/processed/model_out_0.csv.gz: the same as above, except the results from the model that
 uses conservation of the GT/AG nucleotides only.

* data/processed/introns.csv.gz: the CSV table containing information about the unique introns 
 from each annotation analyzed. The table is described in the section below.

* data/processed/transcripts.csv.gz: the CSV table containing one row per transcript in each
 annotation with their support status (well-supported or less-supported).

* data/processed/{gencode,refseq,chess}/{gencode_version,refseq_version,chess_version}_pos.gtf:
 a GTF file containing a subset of the corresponding annotation in which all transcripts are 
 well-supported.

* data/processed/{gencode,refseq,chess}/{gencode_version,refseq_version,chess_version}_neg.gtf:
 a GTF file containing a subset of the corresponding annotation in which all transcripts are
 less-supported.

Intermediate files and directories, absent in the repository and filled by the pipeline: 

* data/interim/{gencode,refseq,chess,mane,random}: directories containing split GTF annotations
  into introns and exons, plus result of the querying MAF alignment.

* data/interim/clinvar: directory containing processed ClinVar data.

* data/interim/gnomad: directory containing processed GnomAD data.

* data/interim/unique_exons: directory containing all the unique exons from the input databases.

* data/interim/mapped_exons: all the unique exons mapped to the other genomes from the alignment.

* data/interim/jobs: directory containing alignment jobs, one job per genome per its chromosome.

* data/interim/results: directory containing the results of the alignment jobs. Each file contains
  human exons that were realigned to this genome. The file contains a sequence of records; each record
  contains 6 lines, and records are separated by a blank line. Forma of the record, line-by-line:

   1. A pair of integers</li>
   2. The GTF line corresponding to the exon in the original annotation</li>
   3. The sequence of the human exon</li>
   4. The alignment string indicating matching basepairs</li>
   5. The sequence of the target genome where the exon mapped</li>
   6. Alignment score, chromosome, strand, start and end coordinates in the target genome, separated by whitespaces</li>
  
  Note that the alignment contains extra 30 basepairs both upstream and downstream of the human exon as
  these alignments are required by the model. For convenience and peer review, the alignments are archived at
  the FTP: ftp://ftp.ccb.jhu.edu/pub/iminkin2/splice-sites-pub/alignment_results/.

* data/interim/identity: the average edit distance and its standard error for all the exons mapped to
  a particular genome. 

Raw data files:

* genomes.txt: a list of 405 genomes we used for the analysis; the list is a subset of all 470 genomes
  used in the original UCSC alignment. This file is only included for the reference; this list is also
  included in the Makefile (variable "all_genomes")

* data/raw/annotation: the directory containing the annotations GTF files.

* data/raw/clinvar: original ClinVar VCF file.

* data/raw/gnomad: original GnomAD VCF files.

* data/raw/human: GRCh38 genome sequence and assembly report.

* data/raw/maf: the UCSC whole genome alignment.

* data/raw/phast: phastCons track downloaded from the UCSC Genome Browser website.

* data/raw/genomes: sequence of 405 genomes we used to realign the exons to. Originally, these files are
  populated and formatted by the Makefile from three sources: UCSCS Genome Browser FTP, NCBI, and DnaZoo.
  Unfortunately, the availability of these genomes change over time, and for the sake of reproducibility,
  we archived all of them at FTP: ftp://ftp.ccb.jhu.edu/pub/iminkin2/splice-sites-pub/genomes/

Source files:

* Makefile: the make pipeline to generate the analysis of the splice sites
 and introns/transcripts. Downloads the human genome, gene annotations, the
 alignment, the sequences involved in it, gnomAD data, ClinVar data, GTEx 
 data and other necessary prerequisites. Then it patches the alignment using
 our original method described in the paper, generates a CSV table containing
 the information about all the splice sites, trains, and runs the classification
 model, and produces an annotated table of splice sites with their support status,
 a CSV table with information about introns (conservation of splice site support, expression 
 support), and a CSV table containing the support status for the transcripts.
 In addition, it produces two files for each annotation: a "pos" file that contains
 all transcripts that are well-supported, and a "neg" file containing "less-supported"
 ones. These files are described below. Only protein-coding and lncRNA genes were included
in the analysis.

* data/src/models/model.py: the python script generating the model_out.csv and containing the
 splice site classification model.

* data/src/model/model_0.py: the python script generating the the model_out_0.csv and containing
 the splice site classification model that uses conservation of the GT/AG nucleotides only.

* src/features/annotate_introns.py: the python script used to generate files introns.csv and
 transcripts.csv

* src/features/generate_table.py: the python script used to generate the input table
 splice_sites.csv for the file implementing the models.


* data/src/reports/*: a collection of scripts generating figures for the manuscript. They
 are not called by the Makefile (yet), and each script is supposed to be run from inside
 its directory.


Description of the fields of the table data/processed/model_out(0).csv
======================================================================

* dataset: the name of the dataset from where the splice site comes
* transcript_id: the id of the transcript that that the splice site belongs to;
 if a site shared by  multiple transcripts, it will appear in the table only
 once
* index: the number of the intron that the site belongs to, introns ordered
 by their starting coordinate on the + strand
* site_type: "a" for acceptor, "d" for donor sites
* gene_type: either "protein_coding" or "lncRNA" depending on the source gene
* inMANE: 1 if the site appears in a MANE transcript, 0 otherwise
* chr: chromosome of the site location
* strand: strand of the site location
* pos: position of the site, 1-based
* cons_GTAG: the number of species in which the canonical dinucleotides GT/AG
 are conserved in 470-species whole-genome alignment
* cons_X: the number of species in which the position with the shift X is
 conserved in 470-species whole-genome alignment. The shift is defined as 
 follows: shifts +0 and +1 correspond to the canonical dinucleotides GT/AT, e.g.
 for donor sites +0 is G and +1 is T. Positive shifts correpond to positions
 downstream of the splice sites, and negative shifts to downsteram positions.
* snp_X: the number of homozygous samples that an SNP from gnomAD v4.0.0
 database contains located at position with the shift X relative to the site.
 Shifts are defined analogously to the previous category. This value is 0 if
 there are no SNPs at this position or it has 0 homozygous samples
* reuse: the number of isoforms of a gene that share this particular site
* well_supported: 0 or 1 depending on whether the site is deemed well_supported by
 the model (1 is conserved)
* prob: probabibility of the site being well_supported, as calculated by the
 model

Description of the fields of the table data/processed/introns.csv
======================================================================
* dataset: the name of the dataset from where the splice site comes
* transcript_id: the id of the transcript
* gene_type: either "protein_coding" or "lncRNA" depending on the source gene
* inMANE: 1 if the whole intron appears in a MANE transcript, 0 otherwise
* chr: chromosome of the intron location
* strand: strand of the intron location
* start, end: starting and ending position of the intron (1-based)
* donor(acceptor)_mane: 1 if the donor (acceptor) site of the intron is in MANE,
 0 otherwise
* donor(acceptor)_well_supported: 1 if the donor (acceptor) site of the intron
 is well_supported, 0 otherwise
* Adipose_Tissue, Adrenal_Gland, Bladder,Blood, Blood_Vessel, Bone_Marrow,Brain,
 Breast, Cervix_Uteri, Colon, Esophagus, Fallopian_Tube, Heart, Kidney, Liver,Lung, 
 Muscle, Nerve, Ovary, Pancreas, Pituitary, Prostate, Salivary_Gland, Skin,
 Small_Intestine, Spleen, Stomach, Testis,Thyroid, Uterus, Vagina: the number of
 reads spanning this particular intron junction in GTEx data, aggregated
 across all samples of a certain tissue

Description of the fields of the table data/processed/transcripts.csv
=====================================================================
* dataset: the name of the dataset from where the splice site comes
* transcript_id: the id of the transcript
* gene_type: either "protein_coding" or "lncRNA" depending on the source gene
* inMANE: 1 if the whole transcript appears in a MANE transcript, 0 otherwise
* chr: chromosome of the transcript location
* strand: strand of the transcript location
* start, end: starting and ending position of the transcript (1-based)
* total_sites: number of splice sites in the transcript
* mane_sites: number of splice sites that also appear in MANE
* well_supported_non_mane_sites: numbe of splice site not appearing in MANE, but
 that are well-supported
* well_supported: 1 if each splice site of the transcript is either well-supported
 or appears in MANE, 0 if at least one site is less-supported (and does not appear
 in MANE)

