#!/usr/bin/perl
## TRANSPOSES FILES (WHICH MUST FIT IN MEMORY)

sub help {
print("TRANSPOSES FILES (WHICH MUST FIT IN MEMORY)\n");

print("usage:   perl transpose.pl <options> inputfile outputfile\n");

print("<options>:\n");
print("-d '<delim>'	:	Specify an input deliminator.  Defailt: ','.\n");
print("-D '<delim>'	:	Specify an output deliminator.  Default: ','.\n");
print("'\\t' is for tabs, ' ' for spaces, ',' for commas, '\w+' for all whitespace.\n");
print("-n		:	Specify that the first row is a header row without a value for the row headers.\n");
print("-w		:	Specify that whitespace is to NOT be stripped from the start and end of each value.\n");
print("EXAMPLE: perl transpose.pl -d '\\t' in.tab out.tab\n");
print("NB perl '\\R' is used for line endings, which should work under all platforms.\n");
die "\n";
}
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

use Switch;
use strict;

###############################
## INPUT:

my $remws=1;
my $norowheader=0;
my $delim=',';
my $delimout=',';
my $inputfile="";
my $outputfile="";

###############################
## ARGUMENT PROCESSING

my $argon=0;
my $rawargon=0;
for (my $i = 0; $i < scalar(@ARGV); ++$i){
	if(@ARGV[$i] eq "-w"){
		$remws=0;
	}elsif(@ARGV[$i] eq "-d"){
		my $tmparg="$ARGV[++$i]";
		my @tmpargarr=split (/'/,$tmparg);
		$delim=$tmpargarr[0];
		print "Using deliminator '$delim'\n";
	}elsif(@ARGV[$i] eq "-D"){
		my $tmparg="$ARGV[++$i]";
		my @tmpargarr=split (/'/,$tmparg);
		$delimout=$tmpargarr[0];
		print "Using deliminator '$delim'\n";
	}elsif(@ARGV[$i] eq "-n"){
		$norowheader=1;
	}else{
        switch($argon){
            case 0 {$inputfile="$ARGV[$i]";}
			case 1 {$outputfile="$ARGV[$i]";}
            else {
               help();
            }
        }
		$argon++;
	}
}
if($inputfile eq "" || $outputfile eq "") {
	help();
}
print "Reading from $inputfile, writing to $outputfile\n";
##############################
## PROGRAM:

              ## (I) GET RECOM-RATE INFO:
my @fulldata;
open(IN,"$inputfile");
my $linelen=0;
my $numlines=0;
while(<IN>)
{
    my $line=$_;
	$line =~ s/\R//g;
    my @linearray=split(/$delim/,$line);
#    my @linearray=split('$delim',$line);
	if($norowheader==1 && $numlines==0) {
		@linearray=("",@linearray);
	}
    if($linelen==0) {$linelen=scalar(@linearray);}
	$numlines++;
	push @fulldata, [ @linearray ];
}
print "Found $numlines lines of length $linelen\n";

## print:
open(OUT,">$outputfile");
for (my $j=0; $j < ($linelen); ++$j) {
	for (my $i=0; $i < ($numlines)-1; ++$i) {	
		if($remws==1){ $fulldata[$i][$j]=trim("$fulldata[$i][$j]");}
		print OUT "$fulldata[$i][$j]$delimout";
	}
	if($remws==1) {$fulldata[$numlines-1][$j]=trim("$fulldata[$numlines-1][$j]");}
	print OUT "$fulldata[$numlines-1][$j]\n";
}

