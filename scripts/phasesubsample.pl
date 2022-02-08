#!/usr/bin/perl
## CONVERTS PHASE (CHROMOPAINTER) FORMAT TO BEAGLE FORMAT
use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

sub help {
print("EXTRACTS SNP RANGE FROM PHASE (CHROMOPAINTER) FORMAT\n");

print("usage:   perl phasesubsample.pl <options> <from> <to> <phasefile> <outputphasefile>\n");
print("Extract the SNPs [from to] inclusive, i.e. 1 2 extracts the 1st and 2nd SNPs.\n");
print("where:\n");
print("<from>:		First SNP to retain (1 is the first snp)\n");
print("<to>:		Final SNP to retain (L is the last snp)\n");
print("<phasefile>:		ChromoPainter/PHASE style SNP file, i.e. \n");
print("<outputphasefile>:       Output phase file\n\n");

print("<options>:\n");
print("-v: Verbose mode\n");
print("NB Compatible with chromopainter and chromopainterv2 phase formats. Updated 6th June 2017 to fix an out-by-one error.\n");
die "\n";
}

###############################
## ARGUMENT PROCESSING

my $from=0;
my $to=-1;
my $newnumsnps=-1;
my $phasefile="";
my $outfile="";
my $verbose=0;

GetOptions ('v|verbose' => \$verbose);
if(@ARGV != 4) {help();}

$from=$ARGV[0];
$to=$ARGV[1];
$phasefile=$ARGV[2];
$outfile=$ARGV[3];

if(!looks_like_number($from) || !looks_like_number($to) || $from<0 || $to<0 || $to<$from){
    die("Invalid arguments: from and to must be specified and must be in the range 1 to L\n");
}
$newnumsnps=$to-$from+1;

####################################
## Define global variables
my @snplocs; # location of the SNPS

my $numsnps=0; # number of SNPS defined in the file
my $numinds=0; # number of individuals defined in the file
my $numhaps=0; # number of haplotypes observed
my $ploidy=-1; # number of haps per ind

####################################
## File IO

## Check we can read the input files
open PHASEFILE, $phasefile or die $!;

## Create output files
open OUTFILE, ">", $outfile or die $!;

####################################
## Functions we need
sub trim($){  # remove whitespace from beginning and end of the argument
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

####################################
## Read the phasefile 

## read the PHASEFILE header
my $skip=1;
my @tmarr;
my $usev2format=0;
while ($skip) {
	my $tmp=<PHASEFILE>;
	my @tmpvals = split(/ /, $tmp);
	if($tmpvals[0] eq "P"){ # found the line with all the SNP locations
		@snplocs= split(/ /, trim($tmp));
		shift @snplocs;
		my $floc=tell(PHASEFILE);
		$tmp=<PHASEFILE>; # read the line of S's, if it exists
		if(substr($tmp, 0, 1) ne "S"){
		    seek(PHASEFILE, $floc, 0); 
		    $usev2format=1;
		}
		$numsnps=trim(pop @tmarr);
		$numinds=trim(pop @tmarr);
		$skip=0;
	}else {
		push @tmarr, $tmpvals[0];
	}
}
if($usev2format==0){
    print "Detected Chromopainter v1 format\n";
    print "Detected $numinds individuals\n";
}else{
    print "Detected Chromopainter v2 format\n";
    print "Detected $numinds haplotypes\n";
}
print "And $numsnps SNPs\n";

if($from>$numsnps || $to>$numsnps){
    die("from or to is larger than the number of SNPs detected ($numsnps) in the file\n");
}

## Print the new header
if($usev2format==0){
    print OUTFILE "0\n";
}
print OUTFILE "$numinds\n$newnumsnps\nP ";
for(my $i=$from-1;$i <$to;++$i){
    print OUTFILE "$snplocs[$i] ";
}
print OUTFILE "\n";

if($usev2format==0){
    for(my $i=$from-1;$i <$to;++$i){
	print OUTFILE "S";
    }
    print OUTFILE "\n";
}

# remaining lines are SNPs
while (my $tmp=<PHASEFILE>) {
    if($verbose){print "Reading haplotype $numhaps\n";}
    ++$numhaps;
    my @tarr=split(//,trim($tmp));
    if(scalar(@tarr)!=$numsnps){
	my $tmp=scalar(@tarr);
	die "Expected $numsnps SNPs on haplotype $numhaps, but received $tmp\n";
    }
    for(my $i=$from-1;$i <$to;++$i){
	print OUTFILE $tarr[$i];
    }
    print OUTFILE "\n";
}
close PHASEFILE;
close OUTFILE;
