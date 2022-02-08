#!/usr/bin/perl
## CONVERTS PHASED BEAGLE OUTPUT TO CHROMOPAINTER-STYLE INPUT FILES

sub help {
print("CONVERTS PHASED BEAGLE OUTPUT TO CHROMOPAINTER-STYLE INPUT FILES\n");

print("usage:   perl beagle2chromopainter.pl <options> beagle_phased_output_file output_filename_prefix\n");

print("where:\n");
print("        (i) beagle_phased_output_file = filename of BEAGLE v3 or less (not vcf!) phased file (unzipped) that contains phased haplotypes\n");
print("        (ii) output_filename_prefix = filename prefix for chromopainter input file(s). The suffixes \".phase\" amd \".ids\" are added\n\n");
print("The output, by default, is in CHROMOPAINTER v2 input format. NOTE THAT ONLY BIALLELIC SNPS ARE RETAINED, i.e. we omit triallelic and non-polymorphic sites.\n");

print("<options>:\n");
print("-J:                 Jitter (add 1) to snp locations if snps are not strictly ascending. Otherwise an error is produced.\n");

print("<further options>   NOTE: YOU ONLY NEED THESE OPTIONS FOR BACKWARDS COMPATABILITY!\n");
print("-v1:                Produce output compatible with CHROMOPAINTER v1, i.e. include the line of \"S\" for each SNP. \n");
print("-f:                 By default, this script produces PHASE-style output, which differs from \n");
print("			   ChromoPainter input which requires an additional first line.  This option creates the correct\n");
print("		           first line for standard fineSTRUCTURE usage (i.e. the first line is \"0\", all other lines are appended)\n\n");

print(" !!! WARNING:  THIS PROGRAM DOES NOT SUFFICIENTLY CHECK FOR MISSPECIFIED FILES. WE ARE NOT ACCOUNTABLE FOR THIS RUNNING INCORRECTLY !!!\n");
print("NOTE: TO USE IN CHROMOPAINTER: You also need a recombination map. Create this with the \"convertrecfile.pl\" or \"makeuniformrecfile.pl\" scripts provided.\n\n");
die "\n";
}


sub trim($){  # remove whitespace from beginning and end of the argument
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

use Switch;
use strict;

###############################
## INPUT:

my $numindsmaxsize=10000000000;     ## ONLY READ IN THIS NUMBER OF BEAGLE INDS AT A TIME (TO MINIMIZE RAM REQUIREMENTS -- MAKES PROGRAM RUN A FACTOR OF (N/$numindsmaxsize) SLOWER THAN IF YOU SET $numindsmaxsize=N, where N is the total number of individuals in your ".haps" BEAGLE output file, but also uses a factor of (N/$numindsmaxsize) less RAM
## NOTE THAT IT ONLY WORKS IF THERE ARE NO OMITTED SNPS

###############################
## ARGUMENT PROCESSING

my $verbose=0;
my $v1=0; ## version 1 mode
my $fsmode=0; ## finestructure mode (i.e. start with an additional line containing 0)
my $jitter=0; ## whether we jitter snp locations

my $Mb = 1000000.0;
my $BEAGLEinfile="";
my $outfilePRE="";

my $argon=0;
for (my $i = 0; $i < scalar(@ARGV); ++$i){
	if(@ARGV[$i] eq "-f"){
	    $fsmode=1;
	}elsif(@ARGV[$i] eq "-v1"){
	    $v1=1;
	}elsif(@ARGV[$i] eq "-v"){
	    $verbose=1;
	}elsif(@ARGV[$i] eq "-J"){
	    $jitter=1;
	}else{
	    switch($argon){
		case 0 {$BEAGLEinfile="$ARGV[$i]";}
		case 1 {$outfilePRE="$ARGV[$i]";}
		else {
		    help();
		}
	    }
	    $argon++;
	}
}

if($outfilePRE eq "" || $argon != 2) {help();}
$outfilePRE =~ s/.phase$//;

##############################
## PROGRAM:

## (II) GET NUMBER OF SITES AND INDS: 
open(IN,"$BEAGLEinfile");
my $line=<IN>;
if(substr($line,0,4) eq "I id"){ # contains the extra spurious line that beagle adds, remove it
    if($verbose) {print("Found ID line\n");}
##    $line=<IN>;
}else{
    die("No header line?");
}
my @linearray=split(/\s+/,$line);
my $totalINDS=(@linearray-2)/2;
my $totalhaps=@linearray-2;
my $numsplits=int($totalINDS/$numindsmaxsize);
if (($numsplits*$numindsmaxsize)<$totalINDS)
{
    $numsplits=$numsplits+1;
}
my $nsites=0; # start at 0 since we already have the first line
while(<IN>)
{
    $line=$_;
    $nsites=$nsites+1;
}

if($verbose){
    print "Detected $nsites sites\n";
}
#print "$numsplits $nsitesFULL $totalINDS\n";
              ## (III) READ IN BEAGLE HAPLOTYPES AND MAKE CHROMOPAINTER HAPLOTYPE INPUT FILE:
open(OUT,">${outfilePRE}.phase");
if($fsmode==1) {
	print OUT "0\n";
}
for ($a=0; $a < $numsplits; $a+=1)
{
    my $startIND=$a*$numindsmaxsize;
    my $endIND=($a+1)*$numindsmaxsize;
    if ($endIND > $totalINDS)
    {
	$endIND=$totalINDS;
    }

                              ## read in:
    open(IN,"$BEAGLEinfile");
    my @rsvec=();
    my @posvec=();
    my @genomat=();
    my $snpcount=0;
    $line=<IN>;
#    if(substr($line,0,4) eq "I id"){ # contains the extra spurious line that beagle adds, remove it
#	$line=<IN>;
#    }
    trim($line);
    @linearray=split(/\s+/,$line);
    shift(@linearray);
    shift(@linearray);
    open(OUTIDS,">${outfilePRE}.ids");
    for(my $i=0;$i<scalar(@linearray);++$i){
	print OUTIDS "$linearray[$i]\n";
    }
    close(OUTIDS);

    while(<IN>)
    {
	$line=$_;
	@linearray=split(/\s+/,$line);
	push(@rsvec,$linearray[1]);
	my $safepos=$linearray[1];
	my @safepos2 =split(":",$safepos);
	$safepos=$safepos2[-1];
	$safepos =~ s/\D//g;
	
#	print("BEFORE: SNP location $linearray[2] lup $lastuniquesnp\n");
	if(scalar(@posvec)>0 &&($safepos <= $posvec[-1]) && $safepos>=0){
	    if(!$jitter){
		die("ERROR: SNPs are not strictly ascending, exiting. Rerun with -J to jitter the SNP locations.\n");
	    }
#	    print("Duplication found: setting $linearray[2] to $posvec[-1]+1\n");
	    $safepos=$posvec[-1]+1;
	}
#	print("AFTER:  SNP location $linearray[2] lup $lastuniquesnp\n");
	if(scalar(@posvec)>0 && ($safepos == $posvec[-1]) ){
	    die("Strange error due to jittering?\n");
	}
	shift(@linearray);
	shift(@linearray);
	my $oksnp=1;

	my %thash;
	my $tval=-1;
	for (my $i=0; $i < scalar(@linearray); ++$i) {
	    if( !exists($thash{$linearray[$i]} ) ){
		$thash{$linearray[$i]}=++$tval;
	    }
	}
	if($tval!=1){ $oksnp=0;
	}else{
	    my @arr=keys(%thash);
	    if(($arr[0]==1 ||$arr[1]==1)&&($arr[0]==0 ||$arr[1]==0)){## 0/1 format; retain it
		$thash{0}=0;
		$thash{1}=1;
	    }
	}
	if($oksnp){
	    push(@posvec,$safepos);
	    for (my $i=$startIND; $i < $endIND; $i+=1)
	    {
		$genomat[(($i-$startIND)*2)][$snpcount]=$thash{$linearray[(2*$i)]};
		$genomat[(($i-$startIND)*2+1)][$snpcount]=$thash{$linearray[(2*$i+1)]};
	    }
	    ++$snpcount;
	}else{
	    if($verbose) {print "Omitting SNP at position $safepos\n";}
	    $nsites--;
	}
    }

                                ## print out:	
    if ($a==0)
    {
	if($v1){
	    print OUT "$totalINDS\n";
	}else {
	    print OUT "$totalhaps\n";
	}
	print OUT "$nsites\n";
	print OUT "P @posvec\n";
	if($v1){
	    for (my $j=0; $j < $nsites; $j+=1)
	    {
		print OUT "S";
	    }
	    print OUT "\n";
	}
    }
    for (my $i=0; $i < (2*($endIND-$startIND)); $i+=1)
    {
	for (my $j=0; $j < $nsites; $j+=1)
	{
	    print OUT "$genomat[$i][$j]";
	}
	print OUT "\n";
    }
}
