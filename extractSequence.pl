#!/usr/bin/perl

##########################################
# extractSequence.pl
#
# Hiroki Takahashi
# Mon 4th Nov, 2013 15:03:53 JST
#
##############################

use strict;
use Getopt::Long;

my $usage = <<USAGE;
Usage: perl extractSequence.pl -i <fasta> -peak <peak position file> -l [10] -comp [false]
        -i        input fasta file
        -peak     peak position file
        -l        width of peak, default=10
        -comp     need reverse complement sequence
        -h        print out this
USAGE

my ($fasta,$peakFile,$len,$comp,$help)= (undef,undef,10,"false",undef);
GetOptions(
	   "i=s"=>\$fasta,
	   "peak=s"=>\$peakFile,
	   "l=i"=>\$len,
	   "comp=s"=>\$comp,
	   "h"=>\$help
);

#$fasta="NC_007779.fasta";
#$peakFile="peakFile.txt";

unless ($fasta){
  print "Input file does not provide yet\n";
  print "\n$usage\n";
  exit(1);
}

if(! -e $fasta){
  print "Input file '$fasta' does not exists\n";
  print "\n$usage\n";
  exit(1);
}

unless ($peakFile){
  print "Input file does not provide yet\n";
  print "\n$usage\n";
  exit(1);
}

if(! -e $peakFile){
  print "Input file '$peakFile' does not exists\n";
  print "\n$usage\n";
  exit(1);
}

## check fasta file
my $count=0;
open(IN, $fasta);
while(my $line=<IN>){
  chomp($line);

  if($line=~/^>/){
    $count++;
  }
}
close(IN);

## check # of sequences
if($count>1){
  print "Your fasta file is multi fasta!!\n";
  print "Single sequence only!!\n";
}

## read fasta file
my $seq=''; my $seqLen=0; my $name='';
open(IN, $fasta);
while(my $line=<IN>){
  chomp($line);

  if($line=~/^>/){
    my @a=split(/ /, $line);
    $name=$a[0];
    $name=~s/^>//;
    next;
  }
  $seq.=$line;
}
close(IN);

$seqLen=length($seq);

## read peak file
my @pos=();
open(IN, $peakFile);
while(my $line=<IN>){
  chomp($line);
  push(@pos, $line);
}
close(IN);

my %count=();
@pos=grep{!$count{$_}++} @pos; ## remove duplicate
my @sort=sort{$a<=>$b} @pos;   ## sort position


my @comp=();
## print out
foreach my $p(@sort){
  my $start=$p-$len-1;
  my $end=$p+$len-1;

  if($p < 0 || $p >= $seqLen){
    print "peak '$p' is not located in the seqence!!\n";
    exit(1);
  }

  if($start<0){
    $start=0;
  }
  if($end>$seqLen){
    $end=$seqLen-1;
  }

  my $peakSeq=substr($seq, $start, ($end-$start+1));
  $peakSeq=~tr/a-z/A-Z/;

  print ">".$name."_".$p."\n";
  print $peakSeq."\n";

  if($comp eq "true"){
    my $compSeq=&getCompSeq($peakSeq);
    push(@comp,  ">".$name."_".$p."_complement\n".$compSeq."\n");
  }
}

## 20140108 edited by HT
## print complement sequence
if($comp eq "true"){
  print "\n# complement\n";
  foreach (@comp){
    print $_;
  }
}


## get reverse complement
##
sub getCompSeq($){
  my $orig=$_[0];

  my @a=split(//, $orig);
  my $comp=join("", reverse @a);
  #$comp=~tr/ACGTN//dc;
  $comp=~tr/AGCT/TCGA/;
  return $comp;
}
