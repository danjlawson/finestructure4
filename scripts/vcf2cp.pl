#!/usr/bin/perl
## CONVERTS PHASED IMPUTE2 OUTPUT TO CHROMOPAINTER-STYLE INPUT FILES


sub help {
print("CONVERTS PHASED SHAPEIT/IMPUTE2 OUTPUT TO CHROMOPAINTER-STYLE INPUT FILES\n");

print("usage:   perl vcf2cp.pl <options> input.vcf output_prefix\n");

print("where:\n");
print("        (i) input.vcf = phased vcf file with GT field.\n");
print("        (ii) output_prefix = filename prefix for chromopainter input file(s). The suffix \".phase\" is added\n\n");
print("The output is in CHROMOPAINTER v2 input format.\n");

print("<options>:\n");
print("-J:                 Jitter (add 1) snp locations if snps are not strictly ascending. By default, an error is produced.\n");
print("-m <val>:           Maximum individuals processed in a single pass. Setting this larger is faster but uses more memory. Default: 500.\n");
print("-g <val>:           Gap between chromosomes. By default, provides an error with multiple chromosome data, but you can set it to 100000000 or some other large value to make whole genome analysis possible.\n");
print("-p <val>:           Ploidy. Default: 2.\n");

print("NOTE: TO USE IN CHROMOPAINTER: You also need a recombination map. Create this with the \"convertrecfile.pl\" or \"makeuniformrecfile.pl\" scripts provided.\n\n");
print(" !!! WARNING:  THIS PROGRAM DOES NOT SUFFICIENTLY CHECK FOR MISSPECIFIED FILES. WE ARE NOT ACCOUNTABLE FOR THIS RUNNING INCORRECTLY !!!\n");
die "\n";
}

use Switch;
use VCF;
use strict;
use warnings;

my $chromosomegap=-1; ## If data from different chromosomes are provided, they are separated by this constant in BP.
my $numindsmaxsize=500;     ## ONLY READ IN THIS NUMBER OF VCF INDS AT A TIME (TO MINIMIZE RAM REQUIREMENTS -- MAKES PROGRAM RUN A FACTOR OF (N/$numindsmaxsize) SLOWER THAN IF YOU SET $numindsmaxsize=N, where N is the total number of individuals in your ".haps" IMPUTE2 output file, but also uses a factor of (N/$numindsmaxsize) less RAM
my $jitter=0; ## whether we jitter snp locations
my $ploidy=2; ## ploidy
my $in="";
my $outPRE="";


my $argon=0;
for (my $i = 0; $i < scalar(@ARGV); ++$i){
  print "$i\n";
  if($ARGV[$i] eq "-J"){
    $jitter=1;
  }elsif($ARGV[$i] eq "-m"){
    $i++;
    $numindsmaxsize = $ARGV[$i];
  }elsif($ARGV[$i] eq "-g"){
    $i++;
    $chromosomegap = $ARGV[$i];
  }elsif($ARGV[$i] eq "-p"){
    $i++;
    $ploidy = $ARGV[$i];
  }else{
    switch($argon){
      case 0 {$in="$ARGV[$i]";}
	case 1 {$outPRE="$ARGV[$i]";}
	else {
	  help();
	}
    }
    $argon++;
  }
}

if($outPRE eq "" || $argon != 2) {help();}
$outPRE =~ s/.phase$//;

print "Options in effect:\n";
print "  input file: $in\n";
print "  output ids file: $outPRE.ids\n";
print "  output phase file: $outPRE.phase\n";
print "  ploidy: $ploidy\n";
print "  jitter: $jitter\n";
print "  individuals per pass: $numindsmaxsize\n";
print "  chromosome gap (-ve for single chromosome analysis only): $chromosomegap\n";


#################################
## Process the header to obtain the individuals
print "Processing Header to identify individuals\n";
my $vcf = VCF->new(file=>$in);
$vcf->parse_header();

## IDs 
my @idnames = @{$vcf->{columns}};
@idnames = @idnames[ 9 .. $#idnames ];
my $ninds =(@idnames);
open(OUTI,">${outPRE}.ids");
for (my $n=0; $n<$ninds; $n+=1)
  {
    print OUTI "$idnames[$n]\n";
  }
print "Found $ninds Individuals\n";
print "Processed header, writing individual ids to ${outPRE}.ids .\n";
my $nhaps=$ninds * $ploidy;
close(OUTI);

#################################
## Initial pass of the chromosomes
print "Processing file to identify SNPs\n";
my @posvec=();
my $nsnps=0;
my $chrom="";
my $snpon=0;
while (my $line  = $vcf->next_line) {
  my @items = split(/\t/,$line);
  $nsnps += 1;
  ## CHECKING CHROMOSOMES
  if ("$chrom" eq ""){
    $chrom="$items[0]";
  }else{
    if("$chrom"!="$items[0]"){
      die("Require all SNPs to be on the same chromosome!\n");
    }
  }
  ## PROCESSING POSITIONS
  if(scalar(@posvec)>0){
    if(($items[1] <= $posvec[-1]) && $items[1]>=0){
      if(!$jitter){
	die("ERROR: SNPs are not strictly ascending, exiting. Rerun with -J to jitter the SNP locations, or remove multi-allelic SNPs.\n");
      }
      print("Duplication found: setting $items[1] to $posvec[-1]+1\n");
      $items[1]=$posvec[-1]+1;
    }
  }
  push(@posvec,$items[1]);
  $snpon += 1;
}
print "Found $nsnps SNPs.\n";
$vcf->close();

#######################


my @rsvec=();
my @genomat=();

## Handling splits
my $numsplits=int($ninds/$numindsmaxsize);
if (($numsplits*$numindsmaxsize)<$ninds)
{
    $numsplits=$numsplits+1;
}
print "Using $numsplits Split(s).\n";

open(OUTP,">${outPRE}.phase");

## Outer loop over the splits
for (my $a=0 ; $a < $numsplits ; $a+=1) {
  my $ta=$a+1;
  print "This is pass $ta of $numsplits.\n";
  @genomat=();
  $vcf = VCF->new(file=>$in);
  $vcf->parse_header();
  
  my $startIND=$a*$numindsmaxsize;
  my $endIND=($a+1)*$numindsmaxsize;
  if ($endIND > $ninds)
    {
      $endIND=$ninds;
    }
  $snpon=0;
  while (my $line  = $vcf->next_line) {
##     print "Processing SNP $snpon of $nsnps\n";
    chomp $line; 
    my @items = split(/\t/,$line);
    ## Loop over individuals
    for(my $n=$startIND;$n<$endIND;$n+=1){
      my @field=split(':',$items[9+$n]);
      my $gt=$field[0];
      my @gtv=split(/\|/, $gt);
      ### check length
      if( ($#gtv != 1)  && ($ploidy>1)){
	die ("Individual $idnames[$n] with field $gt at position index $snpon not phased or diploid?\n")
      }
      my $t1=($n-$startIND)*2;
      $genomat[(($n-$startIND)*$ploidy)][$snpon]=$gtv[0];
      if($ploidy>1){
	my $t2=($n-$startIND)*2+1;
	$genomat[(($n-$startIND)*$ploidy+1)][$snpon]=$gtv[1];
      }
    }
    $snpon += 1;
  }
  $vcf->close();

  print "Processed pass $ta, writing data to ${outPRE}.phase .\n";
  
## print out:	
  if ($a==0)
  {
    print OUTP "$nhaps\n";
    print OUTP "$nsnps\n";
    print OUTP "P @posvec\n";
  }
   for (my $i=0; $i < ($ploidy*($endIND-$startIND)); $i+=1)
    {
      for (my $j=0; $j < $nsnps; $j+=1)
	{
	  print OUTP "$genomat[$i][$j]";
	}
      print OUTP "\n";
    }
}

close(OUTP);
