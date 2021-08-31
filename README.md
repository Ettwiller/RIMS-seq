# About RIMS-seq
RIMS-seq (*Rapid Identification of Methylase Specificity*) simultaneously sequences bacterial genomes and determines m5C methylase specificities using a simple experimental protocol that closely resembles the DNA-seq protocol for Illumina. For more details on RIMS-seq, please read our preprint https://www.biorxiv.org/content/10.1101/2021.03.08.434449v1. \
Note that the sequencing library should be done using a protocol based on adapter ligation.\
(:exclamation:RIMS-seq does not work if the library preparation protocol is based on tagmentation or other technologies:exclamation:)

# Notes before starting
The pipeline requires the following programs to be installed on your system:\
Note: if you have a reference genome and perform you own mapping, SPAdes and BWA are not necessary (see option 2).

- SPAdes (https://github.com/ablab/spades) (Optional)
- BWA (https://github.com/lh3/bwa) (Optional)
- Samtools (http://www.htslib.org/download/)
- Bedtools (https://bedtools.readthedocs.io/en/latest/content/installation.html)
- MoSDi (https://bitbucket.org/tobiasmarschall/mosdi/src/master/)

# RIMS-seq pipeline options

According to your needs, 2 options can be used to run the analysis.\
The option1 performs the genome assembly while the option2 allows you to use your own reference genome for the analysis (see below for more details).\
:exclamation:In case of option2, your reference genome needs to be **exactly** the same as the bacteria sequenced using RIMS-seq. If the exact reference genome is not available, we recommend to use option1, where the genome is derived from the assembly done using RIMS-seq reads.

## Option1: automatic analysis using RIMS-seq.pl (easy)
The `RIMS-seq.pl` script takes paired end reads and outputs a list of mehylase specifities.\
It includes a genome assembly step using SPAdes and uses this assembled genome for the downstream analysis.

## Option2: use your own reference genome (custom)
The RIMS-seq custom pipeline is composed of 3 mains scripts:
- `mapping.pl` 
- `split_mapped_reads.pl`
- `get_motif_all.pl`

Please see the next section 'scripts description and command-line usage' for more details.

# Scripts description and command-line usage
## mapping.pl
This script uses BWA-mem to map the reads to a reference genome. Adapters need to be removed prior using this script.\
Alternatively, you can use your own script for mapping. We have optimized the next analytical steps based on the BWA-mem mapping alorithm - we therefore recommend to use BWA-mem. Reads needs to be mapped in paired-end mode. RIMS-seq requires that the sequencing is done in paired-end mode.  

Usage:
```
perl mapping.pl --fq1 fastq1 --fq2 fastq2 --genome fasta_file --out output_bam_name
```
Required arguments:

`--fq1` : fastq containing Read1

`--fq2` : fastq containing Read2

`--genome` : genome of reference (in fasta format). :exclamation: The reference genome needs to be exactly the same as used in RIMS-seq. We recommend that the genome is derived from the assembly done using RIMS-seq reads to avoid possible issue.

`--out` : name for the output (.bam format)

## split_mapped_reads.pl
This script splits the reads from the bam file generated by BWA. Reads are separated between Read 1 and Read 2 and two files are created accordingly. This is done using samtools mpileup.

Usage:
```
perl split_mapped_reads.pl -bam bam_file -genome genome_file.fasta -mpileup1 output_file1 -mpileup2 output_file2 -Q 10 (DEFAULT 0) -q 20 (DEFAULT 10)
```
Required arguments:

`-bam` : bam file from the mapping step

`-genome` : genome of reference (in fasta format) used for mapping

`-mpileup1` : name for the mpileup1 file generated from Read1

`-mpileup2` : name for the mpileup2 file generated from Read2

Optional arguments : 

`-Q` : minimum base quality for a base to be considered (default=10)

`-q` : minimum mapping quality for an alignment to be used (default=0)
 

## get_motif_all.pl
This script takes the two mpileup files previsouly generated and identifies the methylase specificity using a program called MosDi. 

Usage example:
```
perl get_motif_all.pl --mpileup1 mpileup1 --mpileup2 mpileup2 --qualityscore 35 (DEFAULT 30) --outdir directory --genome genome_file.fasta --flank 8 (DEFAULT 7) --significance 1e-50
```
Required arguments:

`--mpileup1`: mpileup generated using the `split_mapped_reads.pl` program using R1 reads

`--mpileup2`: mpileup generated using the `split_mapped_reads.pl` program using R2 reads

`--genome` : genome of reference (in fasta format) used for mapping 

`--outdir` : output directory

Optional arguments : 

`--qualityscore` : minimum base calling quality score (default=30). It is recommendended that you use a base quality score between 30-35 for Illumina sequencing. Including lower quality base quality scores will increase the background noise :thinking:

`--flank` : the regions around the C to T read variants (default 7bp) <---7bp---C---7bp---> 

`--significance` : motif significance (default = 1e-50) 


## RIMS-seq.pl
This script performs a genome assembly and uses the assembled genome to perform the analysis. It automatically runs the 3 previous scripts (see option 1).

Usage:
```
perl RIMS-seq.pl --fq1 fastq1 --fq2 fastq2
```

Required arguments:

`--fq1`: fastq containing Read1

`--fq2`: fastq containing Read2



## Test the RIMS-seq pipeline
Raw fastq.gz data (5 million reads) have been added to the ```data``` folder. These data correspond to a RIMS-seq experiment on E. coli K12 MG1655 using a 3h incubation. You can use these paired-end files to test the analysis pipeline. The methylated motif found by the pipeline should be C**C**WGG.
```
ecoli_rims_3h_5M.1.fastq.gz
ecoli_rims_3h_5M.2.fastq.gz
```

## How to cite RIMS-seq
RIMS-seq is now published in NAR :
https://doi.org/10.1093/nar/gkab705



