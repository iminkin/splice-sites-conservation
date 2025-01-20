Splice site conservation analysis pipeline
==========================================

Code necessary reproduce the analyses from the paper 
"Conservation assessment of human splice site annotation based on a 470-genome alignment"
by Ilia Minkin and Steven L. Salzberg ([biorxiv link](https://www.biorxiv.org/content/10.1101/2023.12.01.569581v2))

The documentation is still in progress, here is a brief description of the key files:

* Makefile: the make pipeline to generate the analysis of the splice sites
 and introns/transcripts. Downloads the human genome, gene annotations, the
 alignment, the sequences involved in it, gnomAD data, ClinVar data, GTEx 
 data and other necessary prerequisites. Then it patches the alignment using
 our original method described in the paper, generates a CSV table containing
 the information about all the splice sites, trains and runs the classification
 model, and produces an annotated table of splice sites with their support status,
 a CSV table with information about introns (splice site support, expression 
 support), and a CSV table containing the support status for the transcripts.
 In addition, it produces two files for each annotation: a "pos" file that contains
 all transcripts that are well supported, and "neg" file containing "less-supported"
 ones. These files are described below. Only protein-coding and lncRNA genes were included
  in the analysis.

* data/processed/splice_sites.csv.gz: the input file for the regression models; summarizes all
 unique splice from all the dataset.

* data/processed/model_out.csv.gz: the CSV table containing all unique splice sites from
 each annotation we analyzed, along with their properties and support status. The description
 is in a subsection below.

* data/src/models/model.py: the python script generating the model_out.csv and containing the
 splice site classification model.

* data/processed/model_out_0.csv.gz: the same as above, except the results from the model that
 uses conservation of the GT/AG nucleotides only.

* data/src/model/model_0.py: the python script generating the the model_out_0.csv and containing
 the splice site classification model that uses conservation of the GT/AG nucleotides only.

* data/processed/introns.csv.gz: the CSV table containing information about the unique introns 
 form each annotation analyzed. The table is described in a section below.

* data/processed/transcripts.csv.gz: the CSV table containing one row per transcript in each
 annotation with their support status (well-supported or less-supported).

* src/features/annotate_introns.py: the python script used to generate files introns.csv and
 transcripts.csv

* src/features/generate_table.py: the python script used to generate the input table
 splice_sites.csv for the file implementing the models.

* data/processed/{gencode,refseq,chess}/{gencode_version,refseq_version,chess_version}_pos.gtf:
 a GTF file containing a subset of corresponding annotation in which all transcripts are 
 well-supported.

* data/processed/{gencode,refseq,chess}/{gencode_version,refseq_version,chess_version}_neg.gtf:
 a GTF file containing a subset of corresponding annotation in which all transcripts are
 less-supported.

* data/src/reports/*: a collection of scripts generating figures for the manuscript. They
 are not called by the Makefile (yet), and each script is supposed to be ran from inside
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
* Adipose_Tissue,Adrenal_Gland,Bladder,Blood,Blood_Vessel,Bone_Marrow,Brain,
 Breast,Cervix_Uteri,Colon,Esophagus,Fallopian_Tube,Heart, Kidney,Liver,Lung,
 Muscle,Nerve,Ovary,Pancreas,Pituitary,Prostate,Salivary_Gland,Skin,
 Small_Intestine,Spleen,Stomach,Testis,Thyroid,Uterus,Vagina: the number of
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

