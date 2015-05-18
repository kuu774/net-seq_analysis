#!/usr/bin/perl

##########################################
# parse_samtools_mpileup.pl
#
# hiroki.takahashi@chiba-u.jp
#
# Hiroki Takahashi
# Chiba University
#
##########################################

use strict;
use Getopt::Long;

my ($input_file, @line, $refID, $pos, $depth, $ref_base, $base_string, $base_quality);
my ($i, $same_as_ref, $A, $C, $G, $T, $N, $base_number);
my ($insertion, $insert_start, $insert_length, $insert_num, $digits, @insertions);
my ($deletion_start, $deletion_length, $deletion_num, $deletion, @deletions, $j, $d);

#$input_file = shift;

## add by HT
my $usage = <<USAGE;
Usage: perl parse_samtool_mpileup.pl -i <pileupfile> -p [0.7] -d [0]
        -i        input pileup file
        -p        5' tag probability, default=0.7
        -d        threshold of depth, default=0
USAGE

my ($input_file,$prob,$thre)= (undef,0.7,0);
GetOptions(
	   "i=s"=>\$input_file,
	   "p=f"=>\$prob,
	   "d=i"=>\$thre
);

unless ($input_file){
  print "Input file does not provide yet\n";
  print "\n$usage\n";
  exit(1);
}

if(! -e $input_file){
  print "Input file '$input_file' does not exists\n";
  print "\n$usage\n";
  exit(1);
}
## end HERE


open (IN, "$input_file");

print "refID\tposition\tref_base\tdepth\tA\tC\tG\tT\tN\t3prime_rna_end\t3prime_rna_end_minus\t3prime_rna_end_plus\tminus\tplus\tInsertions\tDeletions\n";

while (<IN>) {
    chomp;

    @line = split(/\s+|\t/, $_);
    $refID = $line[0];
    $pos = $line[1];
    $ref_base = $line[2];
    $depth = $line[3];
    $base_string = $line[4];
    $base_quality = $line[5];

    $same_as_ref = 0;
    $A = 0;
    $C = 0;
    $G = 0;
    $T = 0;
    $N = 0;
    $base_number = 0;
    @insertions = ();
    @deletions = ();
    my $fivePrime=0; ## add variable
    my ($forward, $reverse, $endFor, $endRev)=(0,0,0,0); ## add variable

    if($base_string=~/\^\^/){
      print "$base_string\n";
      print "This line should be checked!!\n";
      exit(1);
    }

    ## start HERE by HT
    for ($i=0; $i<= length$base_string; $i++) {
      my $tmp=substr($base_string, $i, 1);

      if($tmp eq '+' || $tmp eq '-'){ ## check INDEL
	my $indelNum=0; my $count=0;

	for($j=$i+1; $j<=length$base_string; $j++) { ## skip INDEl
	  $tmp=substr($base_string, $j, 1);
	  if($tmp=~/\d/){
	    $count++;
	  }
	  elsif($tmp !~ /\d/){
	    last;
	  }
	}
	my $size=substr($base_string, $i+1, $count);
	$i=$i+$count+$size;
      }

      elsif($tmp eq '$'){ ## get $ direction
	$tmp=substr($base_string, $i-1, 1); ## i+1 is bug, left of $ should be checked
	#$i+=1;        ## SKIP this line

	if($tmp =~ /\.|A|T|G|C|N/){ ## 3' end
	  #$forward++; ## SKIP this line
	}
	elsif($tmp =~ /,|a|t|g|c|n/){ ## 5' end, reverse direction
	  $endRev++;
	  #$reverse++; ## SKIP this line
	  $fivePrime++;
	}
      }

      elsif($tmp eq "^"){ ## get ^ direction ^K. K is mapping quality
	$tmp=substr($base_string, $i+2, 1);
	$i+=2;

	if($tmp =~ /\.|A|T|G|C|N/){ ## 5' end
	  $endFor++;
	  $forward++;
	  $fivePrime++;
	}
	elsif($tmp =~ /,|a|t|g|c|n/){  # 3' end
	  $reverse++;
	}
      }
      elsif($tmp =~ /\.|A|T|G|C|N/){
	$forward++;
      }
      elsif($tmp =~ /,|a|t|g|c|n/){
	$reverse++;
      }
    }
    ## until HERE by HT

    for ($i=0; $i<= length$base_string; $i++) {

	if ( substr($base_string, $i, 1) eq "+") {
	  $insert_start = $i;
	  $insert_length = "";
	  $digits = 0;

	   for ($j = $i+1; $j<= (length$base_string); $j++) {
	     if (substr($base_string, $j, 1) =~ /\d/) {       ## e.g. +23AC....

	       $insert_num = substr($base_string, $j, 1);
	       $insert_length = $insert_length.$insert_num;
	       $digits++;
	     }
	     elsif (substr($base_string, $j, 1) !~ /\d/) {
	       $insertion = substr($base_string, ($insert_start + $digits +1), $insert_length );
	       push @insertions, $insertion;
	       last;
	     }
	   }
	  $i=$i+$digits+$insert_length;
	} ## if ( substr($base_string, $i, 1) eq "+") {

      elsif ( (substr($base_string, $i, 1) eq "-") &&( substr($base_string, $i+1, 1) =~ /\d/)) {

	$deletion_start = $i;
	$deletion_length = "";
	$digits = 0;

	   for ($j = $i+1; $j<= (length$base_string); $j++) {

	     if (substr($base_string, $j, 1) =~ /\d/) {          ## e.g. +23AC....

	       #print substr($base_string, $j, 1), "$j ++\n"; ## comment out by HT

	       $deletion_num = substr($base_string, $j, 1);
	       $deletion_length = $deletion_length.$deletion_num;
	       $digits++;
	     }
	     elsif (substr($base_string, $j, 1) !~ /\d/) {
	       $deletion = substr($base_string, ($deletion_start + $digits +1), $deletion_length );
	       push @deletions, $deletion;
	       last;
	     }
	   }
	  $i=$i+$digits+$deletion_length;

	} ## if ( substr($base_string, $i, 1) eq "-") {

	else {
	    $same_as_ref++ if ((substr($base_string, $i, 1) eq ".") ||(substr($base_string, $i, 1) eq ","));
	    $A++ if ((substr($base_string, $i, 1) eq "A") ||(substr($base_string, $i, 1) eq "a"));
	    $C++ if ((substr($base_string, $i, 1) eq "C") ||(substr($base_string, $i, 1) eq "c"));
	    $G++ if ((substr($base_string, $i, 1) eq "G") ||(substr($base_string, $i, 1) eq "g"));
	    $T++ if ((substr($base_string, $i, 1) eq "T") ||(substr($base_string, $i, 1) eq "t"));
	    $N++ if ((substr($base_string, $i, 1) eq "N") ||(substr($base_string, $i, 1) eq "n"));
       }
    }
    @insertions = join ("; ", @insertions) if @insertions;
    @deletions = join ("; ", @deletions) if @deletions;


    if($depth>=$thre){
      if($depth==0){
	print "$refID\t$pos\t$ref_base\t$depth\t$A\t$C\t$G\t$T\t$N\t$fivePrime\t$endFor\t$endRev\t$forward\t$reverse\t@insertions\t@deletions\n";
      }
      elsif($fivePrime/$depth >= $prob){
	print "$refID\t$pos\t$ref_base\t$depth\t$A\t$C\t$G\t$T\t$N\t$fivePrime\t$endFor\t$endRev\t$forward\t$reverse\t@insertions\t@deletions\n";
      }
    }
}

