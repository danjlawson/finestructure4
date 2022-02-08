#!/usr/bin/perl
## 
use Switch;
use strict;
use POSIX;
use List::MoreUtils qw(uniq);

sub help {
  print("CONVERTS FAM (like) FILES INTO CHROMOPAINTER-STYLE ID FILES, OPTIONALLY USING AN EXTERNAL LIST OF INDIVIDUALS TO INCLUDE\n");
  print("FAM file must have at least the ID-column (specify with -fi), see below.\n");
  print("Usage:   idfileFromfam.pl <options> <in.fam> <out.ids>\n\n");
  print("OPTIONS\n");
  print("-l <list.txt>: List of IDs to retain, optionally with a population column. Default: <in.fam>\n");
  print("-i <x> : Use the x-th column of <list.txt> as the individual ID to match the fam file. Default: 0 (first column).\n");
  print("-p <x> : Use the x-th column of <list.txt> as a population ID. Default: 1 (second column). Can be set to the ind id column if pop ids are not available.\n");
  print("-fi <x> : Use the x-th column of the <in.fam> file as the individual ID. Default: 1 (second column).\n");
  print("-fp <x> : Use the x-th column of the <in.fam> file as the population ID. Default: 0 (first column).\n");
  print("-v     : Verbose mode\n");

die "\n";
}

my $listfile="";
my $indcol=0;
my $popcol=1;
my $famindcol=1;
my $fampopcol=0;
my $verbose=0;
my $argon=0;
my $unusedargs=0;
my $remove=0;

while($argon < scalar(@ARGV)){
  if(@ARGV[$argon] eq "-v"){
    $verbose=1;
    $argon++;
  }elsif(@ARGV[$argon] eq "-l"){
    $listfile=$ARGV[++$argon];
    $argon++;
    if($indcol<0) {$indcol=1;}
  }elsif(@ARGV[$argon] eq "-p"){
    $popcol=$ARGV[++$argon];
    $argon++;
  }elsif(@ARGV[$argon] eq "-i"){
    $indcol=$ARGV[++$argon];
    $argon++;
  }elsif(@ARGV[$argon] eq "-fp"){
    $popcol=$ARGV[++$argon];
    $argon++;
  }elsif(@ARGV[$argon] eq "-fi"){
    $famindcol=$ARGV[++$argon];
    $argon++;
  }else{
    $argon++;
    $unusedargs++;
  }
}

if($unusedargs != 2){
  help();
}

my $infile=$ARGV[scalar(@ARGV)-2];
my $outfile=$ARGV[scalar(@ARGV)-1];

if($verbose){
  print "Reading fam-like file $infile, creating id file $outfile\n";
  if($listfile ne "") {
    print "Including only individuals referenced in $listfile in column $indcol.\n";
    if($popcol>0) {print "Using population labels in $listfile from line $popcol\n";}
  }
}

my %indhash = ();


if($listfile eq "") {$listfile=$infile; }

open (LISTIN,"$listfile");
while(<LISTIN>)
  {
    chomp;
    my $line=$_;
    my @linearray=split(/\s+/,$line);
    if(scalar(@linearray)<$indcol || scalar(@linearray)<$popcol){
      my $numc = scalar(@linearray);
      die "Data has only $numc columns but we require max($indcol,$popcol). Are your -p or -i set right?";
    }
    #print "Reading ind $linearray[$indcol] and population $linearray[$popcol]\n";
    $indhash{$linearray[$indcol]} = $linearray[$popcol];
  }
close(LISTIN);

if($verbose){
  my $numinds=scalar(keys(%indhash));
  my $numpops=scalar(uniq(values(%indhash)));
  print "Read $numinds IDS and $numpops populations from $listfile\n";
}
  
open(IN,"$infile");
open(OUT,">$outfile");

while(<IN>)
{
  chomp;
  my $line=$_;
  my @linearray=split(/\s+/,$line);
  my $indid=$linearray[$famindcol];
  my $popid=$linearray[$fampopcol];
  my $keep=0;
  if(exists($indhash{$indid})){
    $popid=$indhash{$indid};
    $keep=1;
  }
  print OUT "$indid $popid $keep\n";
}
close(IN);
close(OUT);
