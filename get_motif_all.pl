#!/usr/bin/perl
use strict;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);

#this program takes two mpileup files (from RIMS-seq) and identifiea the methylase specificity. 
#this program requires bedtools (version xx) and mosdi (https://bitbucket.org/tobiasmarschall/mosdi/src/master/) to be installed on your system.

my $error_sentence = "USAGE : perl $0 --mpileup1 mileup1 --mpileup2 mileup2 --qualityscore 35 (DEFAULT 30) --outdir directory --genome genome_file.fasta\n";

# declare the options :

my $mpileup1; #mpileup file from R1 (output of the split_mapped_reads.pl program)
my $mpileup2; #mpileup file from R2 (output of the split_mapped_reads.pl program)
my $QUALITY_CUTOFF = 30; # base quality (either direct from Illumina or recalibrated). Please keep the quality score high otherwise the noise is higher than the signal. between 30 and 35 is fine. 
my $outdir;  #output dir name. 
my $genome; #genome in fasta format. 
my $flank = 7;
my $significance = "1e-50"; # motif significance (good range between 1e-50 to 1e-100)

#get options :
GetOptions (    "mpileup1=s" => \$mpileup1,
		"mpileup2=s" => \$mpileup2,
	        "qualityscore=s" => \$QUALITY_CUTOFF,
	        "outdir=s" => \$outdir,
		"genome=s" => \$genome,
		"flank=s" => \$flank,
		"significance=s" => \$significance
    ) or die $error_sentence;

#=================================
#if something went wrong in getting the option, notify :
if (!$mpileup1 || !$mpileup2 || !$genome) {die $error_sentence}
#=================================

my $bed1 = $mpileup1; $bed1 =~ s/\.mpileup//; $bed1 =~ s/.*\///g; $bed1 = $outdir."/".$bed1.".bed";
my $bed1_rev = $mpileup1; $bed1_rev =~ s/\.mpileup//; $bed1_rev =~ s/.*\///g; $bed1_rev = $outdir."/".$bed1_rev."rev.bed"; 
my $bed_forground = $bed1."forground.bed";

my $bed2 = $mpileup2; $bed2 =~ s/\.mpileup//; $bed2 =~ s/.*\///g; $bed2 = $outdir."/".$bed2.".bed";
my $bed2_rev = $mpileup2; $bed2_rev =~ s/\.mpileup//; $bed2_rev =~ s/.*\///g; $bed2_rev = $outdir."/".$bed2_rev."rev.bed";
my $bed_background = $bed2."background.bed";

my $out1 = $mpileup1; $out1 =~ s/\.mpileup//; $out1 =~ s/.*\///g; $out1 = $outdir."/".$out1.".fasta";
my $out2 = $mpileup2; $out2 =~ s/\.mpileup//; $out2 =~ s/.*\///g; $out2 = $outdir."/".$out2.".fasta";



#forground
get_bed($mpileup1, $bed1);
get_bed_reverse($mpileup2, $bed1_rev);
my $concatenation1 = "cat $bed1 $bed1_rev > $bed_forground";
system($concatenation1);
get_seq($bed_forground, $out1, $genome, $flank);



#background
get_bed($mpileup2, $bed2);
get_bed_reverse($mpileup1, $bed2_rev);
my $concatenation2 = "cat $bed2 $bed2_rev > $bed_background";
system($concatenation2);
get_seq($bed_background, $out2, $genome, $flank);


#run the motif discovery pipeline
run_motif($out1, $out2,$significance);

sub get_seq {
    my ($file, $out, $genome)=@_;
    my $fai = $genome.".fai";
    my $command = "bedtools slop -i $file -g $fai -b $flank |  bedtools getfasta -s -fi $genome -bed - -fo $out";
    system($command);
}


sub run_motif {
    my ($fasta1, $fasta2,$significance)= @_;
    my $generic = $fasta1;
    $generic =~ s/\.fasta//;

    my @result;

    my $bc = $generic."_background";
    my $final_output = $generic."_signal";
    my $motif_output = $generic."_result";
    #mosdi-utils count-qgrams -p 0.01 -A
    my $c1 = "mosdi-utils count-qgrams -A \"dna\" $fasta2 4 > $bc";
    my $c2 = "mosdi-discovery  -v discovery -q $bc -i -T $significance -M 8,1,0,4 8 occ-count $fasta1 > $final_output";
  
    system($c1);
    clean_error($bc);
    system($c2);
    
    my $motif = qx(grep '^>> ' $final_output | sort -g -k 5 | tail -n 1 | cut -d ' ' -f 2);
    my $motif_ext = qx(grep '^>> ' $final_output | sort -g -k 5 | tail -n 1);
    $motif =~ s/\n//;

    push @result, $motif_ext;
    my $i=1;
    my $tmp1 = $generic.".tmp"; my $tmp1_bc = $generic.".bctmp";
    my $c3 = "cp $fasta1 $tmp1"; my $c4 = "cp $fasta2 $tmp1_bc";
    system($c3); system($c4);
    # if the program keep finding motif, enter this loop :
    while($motif =~/\S+/)
    {
	$i++;
	my $tmp2 = $generic.".tmp2";
	my $final_output_secondaries = $generic."_final_motif".$i;
	my $c1 = "mosdi-utils cut-out-motif -M X $tmp1 $motif > $tmp2";
	my $c2 = "mv $tmp2 $tmp1";
	my $tmp2_bc = $generic.".bctmp2";
	my $c3 = "mosdi-utils cut-out-motif -M X $tmp1_bc $motif > $tmp2_bc";
	my $c4 = "mosdi-utils count-qgrams -A \"dna\" $tmp2_bc 4 > $bc";
	
	my $c5 = "mv $tmp2_bc $tmp1_bc";

	my $c6 = "mosdi-discovery  -v discovery -q $bc -i -T $significance -M 8,2,0,4 8 occ-count $tmp1 > $final_output_secondaries";

	system($c1);
	system($c2);
	system($c3);
	system($c4);
	clean_error($bc);
	system($c5);
	system($c6);
	
	$motif = qx(grep '^>> ' $final_output_secondaries | sort -g -k 5 | tail -n 1 | cut -d ' ' -f 2);
	$motif =~ s/\n//;
	$motif_ext = qx(grep '^>> ' $final_output_secondaries | sort -g -k 5 | tail -n 1);
	
	push @result, $motif_ext;
    }
    my $l = @result;
    open(OUT, ">$motif_output") or die;
    for (my $i=0; $i<$l; $i++)
    {
	my $motif_result = $result[$i];
	print OUT "motif $i = $motif_result\n";
    }
    close OUT;
  
}

sub clean_error {
    my ($bg) = @_;
    open (BG, $bg) or die "can't open $bg\n";
    my $tmp_file = $bg.".tmp";
    open (OUT, ">$tmp_file") or die "can't save into $tmp_file\n";
    foreach my $line (<BG>)
    {
	chomp $line;
	my ($m, $c) = split /\t/, $line;
	if ($c == 0) {print OUT "$m\t1\n";}
	else {print OUT "$line\n";}
	
    }
    close BG;
    close OUT;
    rename $tmp_file, $bg;
}

sub get_nt{
    my ($line) =@_;
    my ($chr,$loc,$ref, $number, $bases,$q1, $q2, $pos) = split /\t/, $line;
    return ($ref);
}

sub get_bed {
    my ($file, $out)=@_; 
    open (FILE, $file) or die;
    open (OUT, ">$out") or die;
    while(my $line = <FILE>)
    {
	$line =~ s/\r|\n//g;
	

	my ($chr,$loc,$ref, $number, $dbases,$q1, $q2, $pos) = split /\t/, $line;
	my $bases = clean_bases($dbases, $number); 
	my $length_bases = length($bases);	    
	my @qualities =split //, $q1;	    
	my @poss =split /,/, $pos;
	my $length_quality = length($q1);
	
	$ref = uc($ref);
	if ($number == $length_bases && $ref =~/[ACTG]/)
	{
	    
	    my @base=split (//,$bases);
	    my $snp = check_snp($bases, $q1, $QUALITY_CUTOFF);
	    if ($snp ==0)
	    {
		for(my $i=0;$i<@base;$i++){
		    my $quality_score = ord($qualities[$i]) - 33;
		    my $pos_read = $poss[$i];
		    if ($quality_score > $QUALITY_CUTOFF && $pos_read != 1 )
		    {
			my $ch=$base[$i];
			my $start_bed = $loc-1;
			my $end_bed = $loc;
			my $type = $ref."_".$ch;
#			if ($ch =~/[atcg]/){print STDERR "$type\n";}
			if ($type eq "C_T"){print OUT "$chr\t$start_bed\t$end_bed\tC_T\t$quality_score\t+\n";}
			if ($type eq "G_a"){
			    print OUT "$chr\t$start_bed\t$end_bed\tG_a\t$quality_score\t-\n";
			}
		    }
		    
		}#end the loop
	    }
	    
	}
    }
    close OUT;
    close FILE;
}


sub get_bed_reverse {
    my ($file, $out)=@_;
    open (FILE, $file) or die;
    open (OUT, ">$out") or die;
    while(my $line = <FILE>)
    {
        $line =~ s/\r|\n//g;
	
	
        my ($chr,$loc,$ref, $number, $dbases,$q1, $q2, $pos) = split /\t/, $line;
        my $bases = clean_bases($dbases, $number);
	my $length_bases = length($bases);
	my @qualities =split //, $q1;
	my @poss =split /,/, $pos;
	my $length_quality = length($q1);
	
	$ref = uc($ref);
        if ($number == $length_bases && $ref =~/[ACTG]/)
        {
	    
            my @base=split (//,$bases);
            my $snp = check_snp($bases, $q1, $QUALITY_CUTOFF);
            if ($snp ==0)
            {
                for(my $i=0;$i<@base;$i++){
                    my $quality_score = ord($qualities[$i]) - 33;
                    my $pos_read = $poss[$i];
                    if ($quality_score > $QUALITY_CUTOFF && $pos_read != 1 )
                    {
			my $ch=$base[$i];
			my $start_bed = $loc-1;
			my $end_bed = $loc;
			my $type = $ref."_".$ch;
#                       if ($ch =~/[atcg]/){print STDERR "$type\n";}                                                              
			if ($type eq "C_t"){print OUT "$chr\t$start_bed\t$end_bed\tC_T\t$quality_score\t+\n";}
			if ($type eq "G_A"){
			    print OUT "$chr\t$start_bed\t$end_bed\tG_a\t$quality_score\t-\n";
			}
                    }
                }#end the loop                                                                                                    
            }

        }
    }
    close OUT;
    close FILE;
}






sub check_snp {
    my ($bases, $q1, $QUALITY_CUTOFF)=@_;
    my $result = 0;
    my $count_total =0; my $count=0;
    my @base=split (//,$bases);
    my @qualities =split //, $q1;
    for(my $i=0;$i<@base;$i++){
	my $quality_score = ord($qualities[$i]) - 33;
	if ($quality_score > $QUALITY_CUTOFF)
	{
	    my $ch=$base[$i];
	    if ($ch =~/[TtAaCcGg]/)
	    {
		$count++;
	    }
	    $count_total++;
	}
    }
    if ($count_total > 0)
    {
	my $percentage = 100*($count / $count_total);
	if ($count>5 && $percentage > 5)
	{
	    $result=1;
	}
    }
    return $result;
}



sub reverse_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCATGCA/;
        return $revcomp;
}

sub clean_bases{
    #remove unwanted sign (insertion / deletion (but retain the + and - and remove the first and last base information)             
    my ($pattern, $n)=@_; #for example : ,$,,                                                                                       
    #remove the ^ and the next character that correspond to the first base of the read.                                             
    $pattern =~ s/\^.//g;
    #remove the dollars sign from the last base of the read.                                                                       
    $pattern =~ s/\$//g;
    
    if ($pattern =~/\+(\d+)/ || $pattern=~/\-(\d+)/)
    {

        while ($pattern=~/\+(\d+)/)
        {
            my $size = $1;
            my $string_to_remove;
            for(my $i=0; $i<$size; $i++)
            {
                $string_to_remove = $string_to_remove.".";
            }
            $string_to_remove = '.\+\d+'.$string_to_remove;


            $pattern =~ s/$string_to_remove/\+/;

        }
        while ($pattern=~/\-(\d+)/)
        {
            my $size = $1;
            my $string_to_remove;
            for(my $i=0; $i<$size; $i++)
            {
                $string_to_remove = $string_to_remove.".";
            }
            $string_to_remove = '.\-\d+'.$string_to_remove;


            $pattern =~ s/$string_to_remove/\-/;
        }
    }
    return $pattern;
}
