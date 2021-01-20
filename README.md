# About RIMS-seq
Here is a set of programs to analyze data derived from RIMS-seq (Rapid Identification of Methylase Specificity). 
RIMS-seq is a new method to simultaneously sequence bacterial genomes and determine m5C methylase specificities. It uses a simple experimental protocol that closely resembles the DNA-seq protocol for Illumina. For more details on RIMS-seq, please read our preprint (link).

# Notes before starting
The pipeline requires the follwing programs to be installed on your system:
- BWA
- Samtools (version xx)
- Bedtools (version xx)
- SPAdes
- MoSDi (https://bitbucket.org/tobiasmarschall/mosdi/src/master/)

# RIMS-seq pipeline options
The RIMS-seq pipeline is composed of 3 mains scripts:
- mapping.pl
- split_mappes_reads.pl
- get_motif_all.pl
  According to your needs, 2 options can be used to run the analysis.  The option 1 performs the genome assembly while the option 2 allows you to use your own reference genome for the analysis (see below for more details).

## Option 1: automatic analysis using RIMS-seq.pl (easy)
The RIMS-seq.pl script automatically runs the 3 previous scripts and outputs a list of mehylase specifities. It includes a genome assembly step using SPAdes and uses this assembled genome for the downstream analysis.

## Option 2: use your own reference genome (custom)
If you do not want to perform a genome assembly and want to use your own reference genome, you will need to run the 3 main scripts detailed above seprately.
First, adapters need to be trimmed and the reads mapped to the reference genome using BWA-mem. Please see the next section 'scripts description' for more details.

# Scripts description
## 1. mapping.pl
Adapters need to be removed and the reads mapped to the reference genome. If not using this script, the mapping should be done using BWA-mem.

## 2. split_mapped_reads.pl
This script splits the reads from the bam file generated by BWA. Reads are separated between Read 1 and Read 2 and two files are created accordingly.

## 3. get_motif_all.pl
This script takes the two mpileup files previsouly generated and identifies the methylase specificity using a program called MosDi. 

## 4. RIMS-seq.pl
This script performs a genome assembly and uses the assembled genome to perform the analysis. It automatically runs the 3 previous scripts (see option 1).
