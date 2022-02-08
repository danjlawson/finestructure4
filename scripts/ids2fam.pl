#!/usr/bin/perl
## 
use Switch;
use strict;
use POSIX;

sub help {
print("CONVERTS CHROMOPAINTER-STYLE ID FILES TO PLINK STYLE FAM FILES\n");
print("usage:   ids2fam.pl <options> ids.txt out.fam\n\n");
print("OPTIONS\n");
print("-v     : Verbose mode\n");
print("-r     : Remove individuals that are supposed to be excluded in a chromoPainter run (third column 0)\n");

die "\n";
}

my $verbose=0;
my $argon=0;
my $remove=0;

for (my $i = 0; $i < scalar(@ARGV); ++$i){
  if(@ARGV[$i] eq "-v"){
    $verbose=1;
    $argon++;
  }
  if(@ARGV[$i] eq "-r"){
    $remove=1;
    $argon++;
  }
}

if(scalar(@ARGV) - $argon != 2){
  help();
}

my $infile=$ARGV[$argon++];
my $outfile=$ARGV[$argon++];

open(IN,"$infile");
open(OUT,">$outfile");

while(<IN>)
{
  chomp;
  my $line=$_;
  my @linearray=split(/\s+/,$line);
  if(scalar(@linearray)!=3) {
    die("Expect standard ChromoPainter id format with 3 columns: id pop inclusion\n");
  }
  if(!$remove || $linearray[2]!=0) {
    print OUT "$linearray[0] $linearray[1] 0 0 0 -9\n";
  }elsif ($verbose){
    print "Removing individual $linearray[0]\n";
  }
}
close(IN);
close(OUT);
