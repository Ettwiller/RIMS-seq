# About RIMS-seq
RIMS-seq (*Rapid Identification of Methylase Specificity*) simultaneously sequence bacterial genomes and determine m5C methylase specificities using a simple experimental protocol that closely resembles the DNA-seq protocol for Illumina. 

Here is a set of programs to analyze data derived from RIMS-seq (Rapid Identification of Methylase Specificity).\
RIMS-seq is a new method to simultaneously sequence bacterial genomes and determine m5C methylase specificities. It uses a simple experimental protocol that closely resembles the DNA-seq protocol for Illumina. For more details on RIMS-seq, please read our preprint (link).

# Notes before starting
The pipeline requires the following programs to be installed on your system: If you have a reference genome and perform you own mapping SPAdes and BWA are not necessary (see option 2).

- SPAdes (https://github.com/ablab/spades) (Optional)
- BWA (https://github.com/lh3/bwa) (Optional)
- Samtools (http://www.htslib.org/download/)
- Bedtools (https://bedtools.readthedocs.io/en/latest/content/installation.html)
- MoSDi (https://bitbucket.org/tobiasmarschall/mosdi/src/master/)

# RIMS-seq pipeline options

According to your needs, 2 options can be used to run the analysis.\
The option 1 performs the genome assembly while the option 2 allows you to use your own reference genome for the analysis (see below for more details).

## Option 1: automatic analysis using RIMS-seq.pl (easy)
The `RIMS-seq.pl` script takes paired end reads and outputs a list of mehylase specifities.\
It includes a genome assembly step using SPAdes and uses this assembled genome for the downstream analysis.

## Option 2: use your own reference genome (custom)
The RIMS-seq custom pipeline is composed of 3 mains scripts:
- `mapping.pl` (optional)
- `split_mappes_reads.pl`
- `get_motif_all.pl`

Please see the next section 'scripts description and command-line usage' for more details.

# Scripts description and command-line usage
## mapping.pl
This script uses BWA-mem to map the reads to a reference genome. Adapters need to be removed prior using this script.\
You can use your own script for mapping. We have optimized the next analytical steps based on the BWA-mem mapping alorithm - we therefore recommend to use BWA-mem.

Usage:
```
perl mapping.pl --fq1 fastq1 --fq2 fastq2 --genome fasta_file --out output_bam_name
```

## split_mapped_reads.pl
This script splits the reads from the bam file generated by BWA. Reads are separated between Read 1 and Read 2 and two files are created accordingly.

Usage:
```
perl split_mapped_reads.pl --bam bam_file -genome genome_file.fasta -mpileup1 output_file1 -mpileup2 output_file2 -Q 10 (DEFAULT 0) -q 20 (DEFAULT 10)
```

## get_motif_all.pl
This script takes the two mpileup files previsouly generated and identifies the methylase specificity using a program called MosDi. 

Usage example:
```
perl get_motif_all.pl --mpileup1 mileup1 --mpileup2 mileup2 --qualityscore 35 (DEFAULT 30) --outdir directory --genome genome_file.fasta --flank 8 (DEFAULT 7) --significance 1e-50
```
Required arguments:

`--mpileup1`: mpileup generated using the `split_mapped_reads.pl` program using R1 reads.

`--mpileup2`: mpileup generated using the `split_mapped_reads.pl` program using R2 reads. 

`--genome` : Genome (in fasta format) used for the mapping. 

`--outdir` : output directory. 

Optional arguments : 

`--qualityscore` : Minimum base calling quality score (default =30). It is recommendended that you use a base quality score between 30-35 for Illumina sequencing. Including lower quality base quality scores will increases the background noise.

`--flank` : The regions around the C to T variants (default 7bp) <---7bp---C---7bp--->. 

`--significance` : Motif significance (default = 1e-50). 


## RIMS-seq.pl
This script performs a genome assembly and uses the assembled genome to perform the analysis. It automatically runs the 3 previous scripts (see option 1).

Usage:
```
perl RIMS-seq.pl --fq1 fastq1 --fq2 fastq2
```

## How to cite RIMS-seq



