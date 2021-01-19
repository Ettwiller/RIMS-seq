## About RIMS-seq
Here is a set of programs to analyze data derived from RIMS-seq (Rapid Identification of Methylase Specificity). 
RIMS-seq is a new method to simultaneously sequence bacterial genomes and determine m5C methylase specificities. It uses a simple experimental protocol that closely resembles the DNA-seq protocol for Illumina. For more details on RIMS-seq, please read our preprint (link)

## RIMS-seq pipeline 
#1. Data pre-processing
######Adapters need to be removed and the reads mappd to the reference genome. Mapping should be done using BWA-mem

#2. split_mapped_reads.pl

#3. get_motif_all.pl
