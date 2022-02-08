#!/usr/bin/perl
## 
use Switch;
use strict;
use POSIX;

sub help {
print("CONVERTS CHROMOPAINTER-STYLE ID FILES TO \n");
print("usage:   ids2ind.pl <options> in.ids out.ind\n\n");
print("OPTIONS\n");
print("-v     : Verbose mode\n");
print("-r     : Remove individuals (ids file third column 0) by setting the .ind population to \"Ignore\"\n");
print("-R     : Remove individuals (ids file third column 0) by omitting them completely.\n");

die "\n";
}

my $verbose=0;
my $argon=0;
my $remove=0;
my $removefull=0;

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
  if(!$removefull || $linearray[2]!=0) {
    my $pop=$linearray[1];
    if($linearray[2]==0 && $remove){
      $pop="Ignore";
    }
    print OUT "\t$linearray[0]\tU\t$pop\n";
  }elsif ($verbose){
    print "Removing individual $linearray[0]\n";
  }
}

close(IN);
close(OUT);
