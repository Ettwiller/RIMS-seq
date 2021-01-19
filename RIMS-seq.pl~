use strict;
use Cwd;


my $DIR = getcwd;
#my $path_to_program ="/mnt/home/ettwiller/laurence/projects/5mC_detection_by_damage/programs/";

my $fq1 = $ARGV[0];
my $fq2 = $ARGV[1];
my $path_to_program = $ARGV[2];
my $phred = 35;

#my @fastq1 = <reads/*val_1.fq>;
#my @fastq2 = <reads/*val_2.fq>;
#programs :

my $mapping_program = $path_to_program."mapping.pl";
my $mpileup_program = $path_to_program."split_mapped_reads.pl";
my $motif_program = $path_to_program."get_motif_all.pl";

#my $assembly_program  = "/home/ettwiller/exe/SPAdes-3.13.0-Linux/bin/spades.py";

#create the necessary directory
my $mapping_dir = "mapping"; my $mpileup_dir = "mpileup"; my $motif_dir = "motif_test4";my $info_dir = "info"; my $assembly_dir = "assembly";  my $context_dir = "context";
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
my $c01 = "mv $assembly_tmp $assembly_genome";
my $c02 = "bwa index $assembly_genome";


my $generic_genome = $genome;
$generic_genome=~ s/.*\///g; $generic_genome=~ s/\.fasta//;
my $bam = $mapping_dir."/".$generic."_".$generic_genome.".bam";
my $mpileup1 = $mpileup_dir."/".$generic."_".$generic_genome."_R1.mpileup";
my $mpileup2 = $mpileup_dir."/".$generic."_".$generic_genome."_R2.mpileup";

my $c1 = "perl $mapping_program --fq1 $fq1 --fq2 $fq2 --genome $genome --out $bam";
my $c2 = "perl $mpileup_program -bam $bam -genome $genome -mpileup1 $mpileup1 -mpileup2 $mpileup2";
my $c3 = "perl $motif_program --mpileup1 $mpileup1 --mpileup2 $mpileup2 --qualityscore $phred --outdir $motif_dir --genome $genome --flank 7\n";

system($c0);
system($c01);
system($c02);
system($c1);
system($c2);
system($c3);
