Splice site conservation analysis pipeline
==========================================

Code for the paper.

TBD.

Description of the fields of the table
======================================

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
* conserved: 0 or 1 depending on whether the site is deemed conserved by
 the model (1 is conserved)
* prob: probabibility of the site being conserved, as calculated by the
 model
