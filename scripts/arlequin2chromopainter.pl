#!/usr/bin/perl
## CONVERTS PHASE (CHROMOPAINTER) FORMAT TO (CHROMOPAINTER)V2 
use strict;
use warnings;
use Getopt::Long;

sub help {
print("CONVERTS FROM ARLEQUIN (SNP) FORMAT TO CHROMOPAINTERv2\n");
print("Warning: only works if labels are in the form x_y\n");
print("usage: perl arlequin2chromopainter.pl <infile> <outputfileroot>\n");

print("with:\n");
print("<infile>:		Arlequin format .arp file, e.g. as created by fastsimcoal2\n");
print("<outputfileroot>:       Output phase & id file root (x.phase and x.ids will be created)\n\n");

print("<options>:\n");
print(":-d <val> : distance between chromosomes, in snps (default: 100000)");
print("-n <name>: predicate for individual label. (Default: IND if the first character is a number, \"\" otherwise)\n");
print("-v: Verbose mode\n");
die "\n";
}

###############################
## ARGUMENT PROCESSING
my $infile="";
my $outfile="";
my $idoutfile="";
my $idhapoutfile="";
my $verbose=0;
my $cdist=100000;
my $indname="default";

GetOptions ('v|verbose' => \$verbose,
	    'd|distance=s' => \$cdist,
	   'n|name=s' => \$indname);
if(@ARGV != 2) {help();}

$infile=$ARGV[0];
$outfile="$ARGV[1].phase";
$idoutfile="$ARGV[1].ids";
$idhapoutfile="$ARGV[1]_hap.ids";

####################################
## Define global variables
my @snplocs; # location of the SNPS

my $lastsnploc=0;
my $numsnps=0; # number of SNPS defined in the file
my $numhaps=0; # number of haplotypes observed

####################################
## File IO

## Check we can read the input files
open INFILE, $infile or die $!;


####################################
## Functions we need
sub trim($){  # remove whitespace from beginning and end of the argument
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

####################################
## Read the file 

# remaining lines are SNPs
my $curpos=0;
while (my $tmp=<INFILE>) {
  if(($curpos==0) && ($tmp =~ /^#[0-9]/ )){
    $tmp =~ s/^.//s;
    $tmp=trim($tmp);
    my @tarr=split(', ',$tmp);
    foreach (@tarr) {
      my $tval = $_+ $lastsnploc;
      push @snplocs,$tval;
    }
    $lastsnploc=$snplocs[-1]+$cdist;
  }elsif(($curpos==1)&&($tmp =~ /_/)) {
    ++$numhaps;
  }
  if($tmp =~/SampleData/){
    $numsnps=scalar(@snplocs);
    $curpos=1;
  }
}

if($verbose){print"Reading $numsnps SNPS for $numhaps haplotypes.\n";}

seek INFILE, 0, 0;
## Create output files
open OUTFILE, ">", $outfile or die $!;
open IDOUTFILE, ">", $idoutfile or die $!;
open IDHAPOUTFILE, ">", $idhapoutfile or die $!;

print OUTFILE "$numhaps\n$numsnps\nP";
for(my $i=0;$i <$numsnps;++$i){
    print OUTFILE " $snplocs[$i]";
}
print OUTFILE "\n";

# remaining lines are the haplotypes
my $hapon=0;
while (my $tmp=<INFILE>) {
  if($tmp =~ /_/){  
    $tmp=trim($tmp);
    my @tarr=split('\t',$tmp);
    $tarr[2]=trim($tarr[2]);
    if($verbose){print "Reading haplotype $tarr[0]\n";}
    if(length($tarr[2])!=$numsnps){
         my $tmp=length($tarr[2]);
         die "ERROR: Expected $numsnps SNPs (from positions) on haplotype $tarr[0], but received $tmp\n";
       }
    print OUTFILE "$tarr[2]\n";
    ++$hapon;

    if($indname eq "default"){
      my $tt=substr($tarr[0],0,1);
      if($tt =~/[0-9]/) {
	my @tids=split('_',$tarr[0]);
	if($tids[1]% 2==0){
	  my $tindid=($hapon/2);
	  print IDOUTFILE "IND$tindid\tPOP$tids[0]\t1\n";
	}
	print IDHAPOUTFILE "HAP$hapon\tPOP$tids[0]\t1\n";
      }else{
	if($hapon % 2 ==0) {print IDOUTFILE "$indname$tarr[0]\n";}
	print IDHAPOUTFILE "$indname$tarr[0]\n";
      }
    }else{
      if($hapon % 2 ==0) {print IDOUTFILE "$indname$tarr[0]\n";}
      print IDHAPOUTFILE "$indname$tarr[0]\n";
    }
  }
}
close INFILE;
close OUTFILE;
close IDOUTFILE;
close IDHAPOUTFILE;

