#!/usr/bin/perl
## REMOVE SINGLETONS OR NON-SNPS FROM PHASE DATA
use strict;
use warnings;

sub help {
print("REMOVE SINGLETONS OR NON-SNPS FROM PHASE DATA\n");

print("usage:   perl phasesscreen.pl <phasefile> <outputphasefile>\n");
die "\n";
}

sub trim($){  # remove whitespace from beginning and end of the argument
        my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
}

###############################
## ARGUMENT PROCESSING

if(@ARGV != 2) {help();}

my $infile=$ARGV[0];
my $outfile=$ARGV[1];


####################################
## Define global variables

my $line;
my @linearray; # temporary array of the current line
my @genomat=();

####################################
## File IO

## Check we can read the input files

open IN,"$infile" or die "Could not open input file $infile\n";
open OUT, ">", $outfile or die "Could not create output file $outfile\n";

## Read the header
## read the PHASEFILE header
my $skip=1;
my @tmarr;
my @posvec;
my $usev2format=0;
my $ninds=0;
my $nsnps=0;
while ($skip) {
	my $tmp=trim(<IN>);
	my @tmpvals = split(/ /, $tmp);
	if($tmpvals[0] eq "P"){ # found the line with all the SNP locations
		@posvec= split(/ /, $tmp);
		shift @posvec;
		my $floc=tell(IN);
		$tmp=<IN>; # read the line of S's, if it exists
		if(substr($tmp, 0, 1) ne "S"){
		    seek(IN, $floc, 0); 
		    $usev2format=1;
		}
		$nsnps=trim(pop @tmarr);
		$ninds=trim(pop @tmarr);
		$skip=0;
	}else {
		push @tmarr, $tmpvals[0];
	}
}
if($usev2format==0){
    print "Detected Chromopainter v1 format\n";
    print "Detected $ninds individuals\n";
}else{
    print "Detected Chromopainter v2 format\n";
    print "Detected $ninds haplotypes\n";
}
print "And $nsnps SNPs\n";


################
## Read the data 
my $nhaps=0;
my @allcodes;
    while(<IN>)
    {
	$line=trim($_);
	@linearray=split(//,$line);
	
	for (my $i=0; $i < $nsnps; $i+=1)
	{
	  $genomat[$nhaps][$i]=$linearray[$i];
	   $allcodes[$i]{$linearray[$i]}=1;
	}
	$nhaps++;
    }
close(IN);


#################
## Count up the number of variable alleles
my @snpcounts;
for my $i ( 0 .. $#allcodes ) {
  for my $role ( keys %{ $allcodes[$i] } ) {
    if($role =~ m/[AGCTagct012]/){
      $snpcounts[$i]+=1;
    }else{
    }
  } 
}

################
## Remove invalid SNPs
my $snpon=0;
my $removed=0;
while($snpon<scalar(@{ $genomat[0] }) ){
  my $snpcount=$snpcounts[$snpon];
    if($snpcount<2 || scalar(@genomat) - $snpcount < 1 || ($snpon>0 && $posvec[$snpon-1]>=$posvec[$snpon]) ){ ## Exclude
      my $on=$snpon+$removed;
      if($snpcount<2)	{print "Removing $on as not enough variability\n";}
      if(($snpon>0 && $posvec[$snpon-1]>=$posvec[$snpon])) {print "Removing $on as duplicated position\n";}
      $nsnps--;
      $removed++;
      for (my $i=0; $i < scalar(@genomat); ++$i){
	splice(@{ $genomat[$i] },$snpon,1);
      }
      splice(@posvec,$snpon,1);
      splice(@snpcounts,$snpon,1);
    }else{
      $snpon++;
    }
}

#################
## Print PHASE format
print "Retained $nsnps SNPs\n";

if($usev2format==0) {print OUT "0\n";}
print OUT "$ninds\n";
print OUT "$nsnps\n";
print OUT "P @posvec\n";
if($usev2format==0) {
    for (my $j=0; $j < $nsnps; $j+=1)
    {
	print OUT "S";
    }
    print OUT "\n";
}

for (my $i=0; $i < $nhaps; $i+=1)
{
    for (my $j=0; $j < $nsnps; $j+=1)
    {
	print OUT "$genomat[$i][$j]";
    }
    print OUT "\n";
}
close(OUT);
