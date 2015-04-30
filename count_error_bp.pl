#! /usr/bin/perl

##########################################
# count_error_bp.pl
#
# hiroki.takahashi@chiba-u.jp
#
# Hiroki Takahashi
# Chiba University
#
##########################################

use strict;
use warnings;
use Getopt::Long;

my $usage = <<USAGE;
Usage: perl count_error_bp.pl -i <sam file> -l <read length>
    -i     input sam file
    -l     read length (count only reads with more than equal 'l')
    -h     print out this
USAGE

my ($input,$len,$help) = (undef,21,undef);
GetOptions(
	   "i=s"=>\$input,
	   "l=i"=>\$len,
	   "h"=>\$help
);

#$input="test2.sam";
#$input="../samFiles/WT.sam";
#$input="../CThi10205.ACA.sam";

if($help){
  print $usage;
  exit(0);
}

unless ($input){
  print "Input file does not provide yet\n";
  print "\n$usage\n";
  exit(1);
}

if(! -e $input){
  print "Input file '$input' does not exists\n";
  print "\n$usage\n";
  exit(1);
}

## do parsing
my $readNum=0;    ## number of counted read with more than 'l'bp
my @count=();     ## count result
my %nucChange=();

## initialize two dimentional list
for(my $i=0; $i<$len; $i++){
  push(@count, 0);
}


open(IN, $input);
while(my $line=<IN>){
  chomp($line);

  if($line=~/^@/){ ## skip header
    next;
  }
  my @a=split(/\t/, $line);
  my ($flag, $cigar, $read)=($a[1], $a[5], $a[9]);

  if($flag!~/0|4|16/){
    print "Flag is $flag. Flag shoud be 0, 4, or 16.\n";
    die "Check Sam file\n";
  }

  if($flag==4){ ## skip unmapped reads
    next;
  }

  my $md='';
  foreach my $aa(@a){ ## get MD:Z tag
    if($aa=~/^MD:Z:/){
      $md=$'; ## 後方参照
    }
  }

  if(length($a[9])>=$len){ ## '=' or '>=' shold be determined
    my ($md2, $align)=&getMappingInfo($cigar, $md, $read);

    if($flag==16){ ## map in complement
      $md2=&getCompSeq($md2);
      $align=&getCompSeq($align);
    }

    ## start count match, mismatch, insertion, or deletion
    my $pos=0; $readNum++;

    for(my $i=0; $i<$len; $i++){
      if($align->[$i] eq '-'){
	next;
      }
      else{
	if($md2->[$i] eq 'M'){
	  $count[$pos]++;
	}
	elsif($md2->[$i] eq 'I'){
	  $nucChange{'I|'.$pos}++;
	}
	elsif($md2->[$i]=~/A|T|G|C/){
	  $nucChange{$md2->[$i].':'.($pos+1).'|'.$align->[$i]}++;
	}
	$pos++;
      }
    }
  }
}
print join(",", @count);
print "\n";
print "$readNum\n";


## output
## position and change pattern
## 1,2,3,...
## A->C
## A->G
## A->T
## C->A
## C->G
## C->T
## G->A
## G->C
## G->T
## T->A
## T->C
## T->G
my @change=('A:C','A:G','A:T','C:A','C:G','C:T','G:A','G:C','G:T','T:A','T:C','T:G',
	    'A:N','C:N','G:N','T:N');
my @row=();
push(@row, "change");
for(my $i=1; $i<=$len; $i++){
  push(@row, $i);
}

print join(",", @row);
print "\n";

for(my $j=0; $j<16; $j++){
  my @out=();
  push(@out, $change[$j]);
  for(my $i=1; $i<=$len; $i++){
    my @tmp=split(/:/, $change[$j]);
    my $key=$tmp[0].":".$i."|".$tmp[1];

    if(exists($nucChange{$key})){
      push(@out, $nucChange{$key});
    }
    else{
      push(@out, 0);
    }
  }

  print join(",", @out);
  print "\n";
}

#while((my $key, my $val)=each(%nucChange)){
#  print "$key\t$val\n";
#}


## get reverse complement
##
sub getCompSeq($){
  my $orig=$_[0];
  my @comp=reverse(@$orig);

  for(my $i=0; $i<=$#comp; $i++){
    $comp[$i]=~tr/ATGC/TACG/;
  }

  return \@comp;
}


## get alignment information
##
## CIGAR:25M1D3M
## MD:Z:25^T1A1
## read:    AACCCAGACGATGATGAAATGCCAA-ATC
## mapping: MMMMMMMMMMMMMMMMMMMMMMMMMTMAM
##
sub getMappingInfo($$$){
  my ($cigar, $md, $read)=@_;
  my @cig=&getCIGAR($cigar);
  my @MD=&getMDtagPlus($md);

  my @Ipos=(); my $judge=0; ## check insertions in CIGAR string
  for(my $i=0; $i<=$#cig; $i++){
    if($cig[$i] eq 'I'){
      push(@Ipos, $i);
      $judge=1;
    }
  }

  my @MD2=@MD;
  if($judge==1){ ## check whether there are insertions or not
    foreach my $pos(@Ipos){
      $MD2[$pos-1].='I';
      my $tmp=join("", @MD2);
      @MD2=split(//, $tmp);
    }
  }

  $judge=0; ## check insertions in CIGAR string
  my @Dpos=();
  for(my $i=0; $i<=$#cig; $i++){
    if($cig[$i] eq 'D'){
      push(@Dpos, $i);
      $judge=1;
    }
  }

  my @read=split(//, $read);
  if($judge==1){ ## check whether there are deletions or not
    foreach my $pos(@Dpos){
      $read[$pos-1].='-';
      my $tmp=join("", @read);
      @read=split(//, $tmp);
    }
  }


  if($cigar=~/^[0-9]+S/){ ## check 5' soft clipping
    my $tmp=$&;
    $tmp=~s/S//;
    $tmp+=0;
    for(my $i=0; $i<$tmp; $i++){
      unshift(@MD2, "S");
    }
  }

  if($cigar=~/[0-9]+S$/){ ## check 3' soft clipping
    my $tmp=$&;
    $tmp=~s/S//;
    $tmp+=0;
    for(my $i=0; $i<$tmp; $i++){
      push(@MD2, "S");
    }
  }
  return (\@MD2, \@read);
}


## parse MD tag into string, i.e., MMMMDMMMMNNN
## 'M': match
## 'D': deletion
## 'ATGC': mismatch
sub getMDtagPlus($){
  my @tag=(); my @str=();
  my $md=$_[0];

  my $md1=$md;
  #$md1=~s/\^[A-Z]+//g;

  while($md1=~/^0[A-Z]/){ ## check 5' end
    push(@tag, $&);
    $md1=$';
  }

  while($md1=~/^[A-Z]0|^[A-Z]|^\d+|^\^[A-Z]+/){
    push(@tag, $&);
    $md1=$'; ## 後方参照
  }

  my $num=0; my $char='';
  foreach (@tag){
    if(/^[A-Z]0$/ || /^[A-Z]$/ || /^0[A-Z]$/){ ## check mismatch
      $num=1;
      $_=~s/[0-9]//;
      $char=$_;
    }
    elsif(/\^/){ ## check deletion
      $_=~s/\^//;
      $num=1;
      $char=$_;
    }
    else{
      $num=$_;
      $num+=0;
      $char='M';
    }

    for(my $i=0; $i<$num; $i++){
      push(@str, $char);
    }
  }

  my @tmp=();
  foreach my $a(@str){
    my @a=split(//, $a);
    foreach my $b(@a){
      push(@tmp, $b);
    }
  }

  @str=@tmp;
  return @str;
}

## parse CIGAR
sub getCIGAR($){
  my @cig=(); my @str=();
  my $cigar=$_[0];
  while($cigar=~/([0-9]+[A-Z])/){
    push(@cig, $&);
    $cigar=$'; ## 後方参照
  }

  ## parse CIGAR string
  foreach (@cig){
    $_=~/[A-Z]/;
    my $char=$&;
    my $num=$`; ## 前方参照
    $num+=0;
    for(my $i=0; $i<$num; $i++){
      push(@str, $char);
    }
  }
  return @str;
}
