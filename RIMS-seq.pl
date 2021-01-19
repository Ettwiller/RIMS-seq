use strict;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);
use Cwd;
use File::Spec;


my($vol,$dir,$file) = File::Spec->splitpath($0);

my $script_dir = $dir;
print STDERR "$script_dir\n";

my $error_sentence = "\n\nUSAGE : perl $0 --fq1 fastq1.fq --fq2 fastq2.fq\n with fastq1.fq and fastq2.fq being the first and second paired-end read files respectively\n\n";

# declare the options upfront :
my $fq1;
my $fq2;

#get all options :
GetOptions (    "fq1=s" => \$fq1,    # the Read1 file 
		"fq2=s" => \$fq2 # the read 2 file
    ) or die $error_sentence;

#=================================
#if something went wrong in getting the option, notify :
if (!$fq1 ||!$fq2) {die $error_sentence}
#=================================

my $phred = 35;





#programs :

my $mapping_program = $script_dir."/mapping.pl";
my $mpileup_program = $script_dir."/split_mapped_reads.pl";
my $motif_program = $script_dir."/get_motif_all.pl";


#create the necessary directory
my $mapping_dir = "mapping"; my $mpileup_dir = "mpileup"; my $motif_dir = "motif";; my $assembly_dir = "assembly";
if ( ! -d $mapping_dir){system("mkdir $mapping_dir");}
if ( ! -d $mpileup_dir){system("mkdir $mpileup_dir");}
if ( ! -d $motif_dir){system("mkdir $motif_dir");}
if ( ! -d $assembly_dir){system("mkdir $assembly_dir");}

my $generic = $fq1;
$generic =~ s/\.fastq//;
$generic =~ s/.*\///g;
my $assembly_file= $assembly_dir."/".$generic."_assembly";
my $assembly_tmp= $assembly_file."/scaffolds.fasta";
my $genome= $assembly_file."/".$generic."_assembly.fasta";
#my $c0 = "/home/ettwiller/exe/SPAdes-3.13.0-Linux/bin/spades.py -1 $fq1 -2 $fq2 -o  $assembly_file";

my $c0 = "spades.py -1 $fq1 -2 $fq2 -o  $assembly_file"; 
my $c01 = "mv $assembly_tmp $genome";
my $c02 = "bwa index $genome";


my $generic_genome = $genome;
$generic_genome=~ s/.*\///g; $generic_genome=~ s/\.fasta//;
my $bam = $mapping_dir."/".$generic."_".$generic_genome.".bam";
my $mpileup1 = $mpileup_dir."/".$generic."_".$generic_genome."_R1.mpileup";
my $mpileup2 = $mpileup_dir."/".$generic."_".$generic_genome."_R2.mpileup";

my $c1 = "perl $mapping_program --fq1 $fq1 --fq2 $fq2 --genome $genome --out $bam";
my $c2 = "perl $mpileup_program -bam $bam -genome $genome -mpileup1 $mpileup1 -mpileup2 $mpileup2";
my $c3 = "perl $motif_program --mpileup1 $mpileup1 --mpileup2 $mpileup2 --qualityscore $phred --outdir $motif_dir --genome $genome --flank 7\n";

#system($c0);
#system($c01);
#system($c02);
system($c1);
system($c2);
system($c3);
