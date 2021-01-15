use strict;
use Cwd;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);

my $error_sentence = "USAGE : perl $0 --fq1 fastq1 --fq2 fastq2 --genome fasta_file --out bam\n";

# declare the options upfront :
my $fq1;
my $fq2;
my $genome;
my $bam;

#get all options :
GetOptions (    "fq1=s" => \$fq1,    # the Read1 file 
		"fq2=s" => \$fq2,# the read 2 file
		"genome=s" => \$genome, #the genome_file
		"out=s" => \$bam #output file name.
    ) or die $error_sentence;

#=================================
#if something went wrong in getting the option, notify :
if (!$fq1 ||!$fq2|| !$bam || !$genome) {die $error_sentence}
#=================================


my $DIR = getcwd;
my $generic = $bam;
$generic =~ s/\.bam//;
my $sam = $generic.".sam";
my $bam = $generic.".bam";

my $index = $genome.".bwt";
if (!$index){my $command0 = "bwa index $genome"; system($command0);}

my $command1 = "bwa mem $genome $fq1 $fq2  > $sam";
my $command2 = "samtools view -bS $sam | samtools sort -o $bam";
my $command3 = "samtools index $bam";

system($command1);
system($command2);
system($command3);





