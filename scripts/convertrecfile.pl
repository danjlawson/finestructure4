#!/usr/bin/perl
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   http://www.gnu.org/licenses

package convertrecfile;
use Getopt::Std;
use Switch;
use strict;
use POSIX;


## help
sub showHelp {
print "-----convertrecfile.pl, create recombination maps for phase files from other maps.\n";
print "Copyright (C) 2014 Daniel Lawson (dan.lawson\@bristol.ac.uk) licenced under GPLv3\n";
print "This is free software with NO WARRANTY, you are free to distribute and modify; see http://www.gnu.org/licenses\n\n";

print "Usage: ./convertrecfile.pl <MAJOR MODE> <options> phasefile inrecfile outputrecfile\n";
print "phasefile is a valid chromopainter or chromopainter v2 inputfile ending in .phase\n";
print "inrecfile is a recombination file specified in one of the formats specified in <mode>\n";
print "outputrecfile will be a valid recombination file for use with ChromoPainter.\n";
print "MAJOR MODES: specified with -M. (Shortest unambiguous mode option will work)\n";
print " -M: <val>:      Specify the major mode. <val> can be:\n";
print "         hapmap: The hapmap format is specified as 4 columns: chromosome, Position(BP) Rate(cM/Mb) Map(cM)\n";
print "                 This uses columns 2 and 4 to reconstruct the map.\n";
print "         plink:  The more recent format is specified as 4 columns: chromosome, (unused), Map(cM) Position(BP)\n";
print "                 This uses columns 2 and 4 to reconstruct the map.\n";
print "         plain:  (default) Assumes that the data are specified in 2 columns, Position(BP) Rate(M/b)\n";
print "                 This is the mode assumed chromopainter (note: the rate is Morgans per base).\n";
print "Other important options:\n";
print " -v:	        Verbose mode.\n";
print " -h:             This help.\n";
print " -H:             Detailed help on the wide variety of different options, including different column\n";
print "                 separators, different units, reading of Culmulative vs non-culmulative distributions,\n";
print "                 and handling maps that do not cover the full range of the SNPs.\n";
print "EXAMPLE: ./convertrecfile.pl -M hap my_chr1.phase genetic_map_GRCh37_chr1.txt my_chr1.recombfile\n";

}

sub showHelpOptions {
print "\nOPTIONS: \n";
print " -I <val>:	Specify the range over which the recombination map is specified.  <val> can be:\n";
print "		snp:	(default) map is specified in bases and corresponds to the range in the phasefile.\n";
print "		unit:	map is specified over the range (0,1), over the range corresponding to the -N <val> parameter.\n";
print "			(NOTE: actually it is valid over any range (a,b) which is scaled into (N1,N2) for computation - see -N option).\n";
print " -N <val>,<val>:	Underlying range of bases (N1,N2) in the phasefile. Only used by \"-I unit\".\n		By default, assumed to be the observed range of the SNPs if not specified.\n";
print " -U <val>:	Unit of measure for the input map, which is output as MORGANS/BASE. Can be one of the following (specified as <u1>/<u2> or <u1> if a CDF or denominator is per base):\n";
print "		Morgans/Base:	(default) recombination is in morgans/base, and doesn't need rescaling.\n";
print "		centi-Morgans/Base:	recombination is in centi-morgans/base, and  will be multiplied by 1/100.\n";
print "		Morgans/MegaBase:	M/MB: Will be multiplied by 1/1000000.\n";
print "		centi-Morgans/MegaBase:	cM/MB: Will be multiplied by (1/100)*(1/1000000).\n";
print "		<numeric>:	specify a multiplier.\n";
print "	-c :            Output a culmulative map, as required by MOSAIC, instead of the default.\n";
print " -D <val>:	Distribution function; how the distribution is specified. <val> can be one of:\n";
print "		cdf:	Culmulative Distribution Function; map is a piecewise-linear (possibly unnormalised) culmulative function.\n";
print "		pdf:	(default) Probability Density Function: map is a piecewise-constant (possibly unnormalised) function.\n";
print " -T <val>:	Define how total recombination is to be calculated.\n		Note: you can compute this after the fact via EM in ChromoPainter (e.g. -in -iM -i 10). <val> can be:\n";
print "		norm:	(default) Use the normalisation of the distribution to calculate it. Note that if you don't have a normalised distribution, you can still use this mode, but then you must use ChromoPainter's EM estimation of Ne.\n";
print "		<num>:	The total recombination in Morgans is the number after x.  This overrides the -U option.\n";
print " -t <a,b>:       The fields to be used as the position, and the recombination rate, with the first column indexed as 1 (Default: 1,2).\n";
print " -d \"<val\":    Set the separator for fields in the (input) genetic map. The default is to use whitespace (\"\\\\s+\"); don't forget to escape backslashes.\n";
print " -i <val>:       Set the mode for the input SNP locations. By default this is \"phase\", but \"text\" is also valid, which expects one SNP location per line (instead of a phase file).\n";
print " -s:     Set zero recombination rate outside the range ofthe genetic map. The default is to use the mean rate over the whole map.\n";
print " -v:	Verbose mode.\n";
print " -h:	Show basic help message.\n";
print " -H:	Show this help message.\n";
print " -A:	Detailed help on how the assumptions behind the cdf and pdf specification.\n\n";
print "Example to reproduce the behaviour of hapmap mode:\n";
print "./convertrecfile.pl -U c -t 2,4 -D CDF my_chr1.phase genetic_map_GRCh37_chr1.txt my_chr1.recombfile\n";
print "Example to reproduce the behaviour of hapmap mode from the rates:\n";
print "./convertrecfile.pl -U c/M -t 2,3 -D PDF my_chr1.phase genetic_map_GRCh37_chr1.txt my_chr1.recombfile\n";
die "\n";
}

sub showHelpAssumptions {
print " *** Detailed assumptions ***\n";
print "ChromoPainter requires input of the form:\n(x_i\tdelta_i)\n";
print "where x_i is the SNP location in BP and delta_i is the recombination rate PER BASE to the NEXT SNP, i.e. SNP i+1 .\n";
print "All pdf's are assumed to be in this form.  Therefore the last value (for the final SNP L) is ignored, and delta_L=0 is typically assumed.\n";
print "All cdf's are assumed to be a piecewise-smooth linear function, i.e. to define the points (x_i\tR_i)\n";
print "Here R_i is the culmulative recombination rate function.  Therefore the FIRST value R_0=0, although again this is ignored in practice.\n";
print "The general formula used to convert between them is therefore:\n";
print "\tR_{i+1} = R_{i} + (x_{i+1} - x_{i}) * delta_{i}\n";
print "or:\tdelta_i = (R_{i+1} - R_{i})/(x_{i+1} - x_{i})\n";
print "\nInterpolation is performed on the basis of this CDF.\n";
print "\nRemember that values should be specified in Morgans.\n";
die "\n";
}

sub sethapmap { # sets things up for processing Hapmap style data, just as "-U c -t 2,4 -M CDF" would do 
    our $unit="centimorgans";
    our $cdvalmultiplier=1.0/100.0;
    our $locfield=2;
    our $recfield=4;
    our $req2fields=0;
    our $dmode="cdf";
}

sub setplink { # sets things up for processing Plink style data, just as "-U c -t 4,3 -M CDF" would do 
    our $unit="centimorgans";
    our $cdvalmultiplier=1.0/100.0;
    our $locfield=4;
    our $recfield=3;
    our $req2fields=0;
    our $dmode="cdf";
}

sub setplain { # sets things up for processing Hapmap style data, just as "-U c -t 2,4 -M CDF" would do 
    our $unit="morgans";
    our $cdvalmultiplier=1.0;
    our $locfield=1;
    our $recfield=2;
    our $req2fields=1;
    our $dmode="pdf";
}
#####################################
## Some utility functions
## Hacked way to do infinity
my $inf = 9**9**9;
my $neginf = -9**9**9;
my $nan = -sin(9**9**9);

sub isinf { $_[0]==$inf || $_[0]==$neginf }
sub isnan { ! defined( $_[0] <=> $inf ) }
sub signbit { substr( sprintf( '%g', $_[0] ), 0, 1 ) eq '-' } # useful for detecting negative zero

sub max ($$) { $_[$_[0] < $_[1]] }
sub min ($$) { $_[$_[0] > $_[1]] }
sub trim($){
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
sub isint{
    if (/\D/) { return(0); }
    return(1);
}

##################################################
## Reading in a phase file to obtain SNP locations
sub readPhase($){
## read the PHASEFILE header
    my $phasefile=shift;
    open PHASEFILE, $phasefile or die $!;
    my $skip=1;
    my @tmarr;
    my @snplocs;
    my $usev2format=0;
    my $numinds;
    my $numsnps;
    while ($skip) {
	my $tmp=<PHASEFILE>;
	my @tmpvals = split(/ /, $tmp);
	if($tmpvals[0] eq "P"){ # found the line with all the SNP locations
		@snplocs= split(/ /, $tmp);
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
    return(@snplocs);
}

sub readText($){
    my $phasefile=shift;
    print "Reading from text file $phasefile\n";
    open PHASEFILE, $phasefile or die $!;
    my @tarr;
    while (my $tmp=<PHASEFILE>) {
	$tmp=trim($tmp);
	print "SNP $tmp\n";
	if($tmp eq $tmp+0.0) {push @tarr, $tmp;}
    }
    close PHASEFILE;
    return(@tarr);
}

#####################################
## Define the variables 
my $maptype="snp";
my $baserangeA=$inf; # range conversion for -I unit and -N ..., minimum SNP location
my $baserangeB=$neginf;# range conversion for -I unit and -N ..., maximum observed SNP location
my $reclocA=$inf;# range conversion for -I unit, minimum observed recombination location
my $reclocB=$neginf;# range conversion for -I unit, maximum observed recombination location
our $unit="morgans";
our $dmode="pdf";
my $totrec="norm";
my $verbose=0;
our $cdvalmultiplier=1.0; # recombination map multiplier for converting from the native map scale to Morgans/base.
my $sep="\\s+"; # Separator for the input genetic map
our $locfield=1; # first field of the recombination file is assumed to be the location
our $recfield=2; # second field of the recombination file is assumed to be the recombination rate
our $req2fields=1; # by default, we require 2 fields per line. But if column numbers are specified we allow any number
our $inputmode="phase"; # by default we read phase files for the map locations
my $startfix=1; # if 1, we try to prevent a zero recombination rate up to the first SNP in the map
my $zerorecrate=0; # replace any zero recombination rates with this value (the mean if $startfix==1, 0 otherwise)
my $outputcum=0; # output a cumulative map instead of the default.
## Extract options
my %Options;
my $continue = getopts('I:N:U:M:T:t:d:D:i:hcsHAv', \%Options);


for my $key ( keys %Options ) {
        my $value = $Options{$key};
		switch ($key) {
		    case "M" {
			if($value =~ /^[hH]/){ sethapmap();
			}elsif($value =~ /^[pP]lain/){ setplain();
			}elsif($value =~ /^[pP]/){ setplink();
			}else{print "Invalid mode \"$value\": try \"hapmap\", \"plink\" or \"plain\"\n";}
		    }
			case "I" {
				if($value =~ /^[sS]/){ $maptype="snp";
				}elsif($value =~ /^[uU]/){ $maptype="unit";
				}else{
					print "ERROR: Unrecognised -I option: \"$value\"\n";
					$continue=0;
				}
			}
			case "N" {
				my @tmpvals = split(/,/, $value);
				if(scalar(@tmpvals)!=2){print "Invalid -N option! $value does not have two values separated by a comma.\n"; $continue=0;
				}else{
					$baserangeA=$tmpvals[0]+0;
					$baserangeB=$tmpvals[1]+0;
				}
			}
			case "U" {
			    my @vararr=split('/',$value);
			    if(scalar(@vararr)==1){ push(@vararr,"base");}
			    if ( ($value =~ /^[mM]/ || $value =~ /^[cC]/  || $value =~ /^[bB]/) && scalar(@vararr)!=2) {
				print("Specify recombination units (-U) either as a number, or \"a/b\" where \"a\" is either morgans/centi-morgans and \"b\" is base/megabase\n");
				$continue=0;
			    }else{
				if(scalar(@vararr)!=2){ 
					$unit="user";
					$cdvalmultiplier=$value+0; # convert to numeric
				}else{
				    if($vararr[0] =~ /^[mM]/){ $unit="morgans";
				    }elsif($vararr[0] =~ /^[cC]/) { 
					$unit="centimorgans"; 
					$cdvalmultiplier=1.0/100.0;
				    }
				    if($vararr[1] =~ /^[bB]/){ $unit="$unit/base";
				    }elsif($vararr[1] =~ /^[mM]/) { 
					$unit="$unit/Mbase"; 
					$cdvalmultiplier=$cdvalmultiplier / 1000000;
				    }
				}
				if($cdvalmultiplier+0.0 == 0.0){
				    $continue=0;
				    print "ERROR with -U flag: recombination multiplier is 0\n";
				}else{
				    print "Using recombination rate multiplier $cdvalmultiplier\n";
				}
			    }
			}
			case "D" {
				if($value =~ /^[pP]/){ $dmode="pdf";
				}elsif($value =~ /^[cC]/){ $dmode="cdf";
				}else{
					print "ERROR: Unrecognised -M option: \"$value\"\n";
					$continue=0;
				}
			}
			case "T" {
				if($value =~ /^[nN]/){ $totrec="norm";
				}else{
					$totrec=$value+0;
					if($totrec<=0) {
						print "ERROR: Invalid total recombination rate: must be >0.\n";
						$continue==0;						
					}
				}
			}
			case "t" {
			    my @arrval = split(/,/ , $value);
			    if(scalar(@arrval) != 2 || ((!isint($arrval[0]) || !isint($arrval[1])))) { 
				print "ERROR: Invalid specification of \"-t\": use a single comma to separate the column numbers (e.g. \"-t 3,4\").\n";
				$continue==0;						
			    }else{
				$locfield=$arrval[0];
				$recfield=$arrval[1];
				$req2fields=0;
			    }
 			}
		    case "i" { 
			if($value =~ /^[pP]/){ $inputmode="phase";
			}elsif($value =~ /^[tT]/){ $inputmode="text";
			}else{print "Invalid input mode -i $value\n"; $continue=0;}
		    }
			case "c" { $outputcum=1;}
			case "d" { $sep =$value;}
			case "s" { $startfix=0;}
			case "v" { $verbose=1;}
			case "h" { $continue=0;}
			case "H" { showHelpOptions();}
			case "A" { showHelpAssumptions();}
			else {print "ERROR: Invalid option \"$key\".\n";
			      $continue==0;}
		}
}

if($#ARGV+1 !=3) {
	$continue=0;
}
if($continue == 0) {
	showHelp();
	die "\n";
}

#########################################################################
## Set up the file names
my $phasefile = $ARGV[0];#'myres-1-1.phase';	
my $inputrecfile = $ARGV[1];#'input.recfile';	
my $outputfile = $ARGV[2];#'test.recombfile';

#if($verbose) {
my $options="OPTIONS: -I $maptype -U $unit -M $dmode -T $totrec -t $locfield,$recfield";
if(!isinf($baserangeA)){$options="-N $baserangeA:$baserangeB";}
if($startfix){ $options="$options -s";}
if($verbose){ $options="$options -v";}
print "$options\n";
print "Effective unit= $cdvalmultiplier\n";
print "FILES: PHASE FILE: $phasefile, REC FILE: $inputrecfile -> REC FILE: $outputfile\n";
#}

## Set up the outputfile
open OUTPUTFILE, ">", $outputfile or die $!;
print OUTPUTFILE "start.pos recom.rate.perbp\n";

#########################################################################
## Read in the SNPS
my @snps; 
switch($inputmode) {
    case "phase" {@snps=readPhase($phasefile);}
    case "text" {@snps=readText($phasefile);}
    else {die "Invalid input mode $inputmode.\n"; }
}

my $numsnps=scalar(@snps);
if ($verbose) {print "Found $numsnps SNPS in phase file\n";}
my $finalSNP=$snps[$numsnps-1];

if($verbose) {
    print "SNP range: ($snps[0],$snps[-1])\n";
}
#########################################################################
## decide on the range of the recombination map, if it is specified in a unit or undefined interval
if($maptype eq "unit"){
	if(isinf($baserangeA)) {
		$baserangeA=$snps[0];
		$baserangeB=trim($snps[$numsnps-1]);
		print "Setting SNP range to $baserangeA -- $baserangeB\n";
	}
	if($baserangeA>$snps[0] || $baserangeB < trim($snps[$numsnps-1])) {
		die "Invalid \"-N\" option specified: SNPs exist outside the range of the recombination map!\nrecombination range ($reclocA -- $reclocB), SNP range ($baserangeA -- $baserangeB)\n";
	}
	## GET THE RANGE OF THE INPUT FILE
	if(isinf($reclocA)){
		open INRECFILE, $inputrecfile or die $!;
		while (my $tmp=<INRECFILE>) {
			my @tmpvals = split(/ /, $tmp);
			my $tpos=trim($tmpvals[0]);
			if($tpos>$reclocB) {$reclocB = $tpos;}
			if($tpos<$reclocA) {$reclocA = $tpos;}
		}
		close INRECFILE;
		print "Observed range is $reclocA -- $reclocB\n";
	}
}

#########################################################################
## Process the recombination map by reading it in as a cdf.
open INRECFILE, $inputrecfile or die $!;
my @cplocs; #locations of changes in the value of the CDF gradient
my @cpvals; #value of the CDF at the locations in cplocs
my $cpcounter=0; # counter of the the line we're on (for validation)
my $lastval=0; # last recombination cdf value we observed
my $lastrec=0; # last recombination pdf value we observed
my $lastloc=-1; # last recombination change point location we observed
my $process=1; # only process valid lines where we can convert everything to numeric

while (my $tmp=<INRECFILE>) {
	$cpcounter++;
	$process=1;
	my @tmpvals = split(/$sep/, $tmp);
	my $tmpvalsv=scalar(@tmpvals);
	if(scalar(@tmpvals) < $recfield || scalar(@tmpvals) < $locfield) {$process=0;
	}elsif($req2fields==1 && scalar(@tmpvals)!=2){ $process=0;
	}elsif($tmpvals[$locfield-1] =~ /[a-df-zA-DF-Z]/){$process=0;
	}elsif($tmpvals[$recfield-1] =~ /[a-df-zA-DF-Z]/){$process=0;}
	if($process){
		my $recloc=trim($tmpvals[$locfield-1]);
		if("$maptype" eq "unit") {
			$recloc=$baserangeA + ($baserangeB-$baserangeA)*($recloc-$reclocA)/($reclocB-$reclocA);
		}
		if($dmode eq "cdf"){
			push @cplocs, $recloc;
			push @cpvals, $cdvalmultiplier * trim($tmpvals[$recfield-1]);
		}elsif($dmode eq "pdf"){# convert back to cdf
		    if($lastloc!=-1){
			push @cpvals, $lastval;
			push @cplocs, $lastloc;
			my $val=$lastrec*($recloc-$lastloc) * $cdvalmultiplier;
			$lastval+= $val;
#				print "$lastloc, $lastval -- $tmpvals[$locfield-1], $tmpvals[$recfield-1]\n";
		    }
		    $lastrec=trim($tmpvals[$recfield-1]);
		    $lastloc=$recloc;
		}
	}else{
		print "WARNING: NOT PROCESSING LINE $cpcounter: $tmp\n";
	}
}
if($dmode eq "pdf" && $process){ # add te final line
    push @cpvals, $lastval;
    push @cplocs, $lastloc;
}
close INRECFILE;
if(scalar(@cplocs)==0) {
    die "ERROR: No valid lines found in the location of SNPs. Something is wrong.\n";
}

#########################################################################
## Normalisation and tweaking to the recombination rates
if($totrec ne "norm"){
	$cdvalmultiplier=$totrec/$cpvals[-1]; # maximum value of the CDF ovserved
	print "Setting total recombination to $totrec MORGANS PER BASE using the multiplier $cdvalmultiplier\n";
	for(my $i = 0; $i < scalar(@cplocs); ++$i){
		$cpvals[$i]=$cpvals[$i]*$cdvalmultiplier;
	}
}
if($startfix){ # try to set the first value to the rate of the second value
    $zerorecrate = $cpvals[scalar(@cplocs)-1]/($cplocs[scalar(@cplocs)-1]-$cplocs[0]);
    my $t1=$cpvals[scalar(@cplocs)-1];
    my $t2=$cplocs[scalar(@cplocs)-1];
    my $t3=$cplocs[0];
    if($verbose){print "Using mean recombination rate of $zerorecrate for any SNPs with zero recombination rate.\n";}
}

#########################################################################
## We now have a working valid recombination file
if ($verbose) {
	print "Map range: ($cplocs[0],$cplocs[-1])\n";
	print "Recombination CDF calculated as follows:\n";
	print "Location, CDF\n";
		for(my $i = 0; $i < scalar(@cplocs); ++$i){
		print "$cplocs[$i], $cpvals[$i]\n";
	}
	print "\n";
}


#########################################################################
## Loop over all SNPs to calculate their recombination position
my $recinterval=-1; # interval of the SNP cdf we are currently in
my $cdfdiff=0; # current CDF difference from last SNP
#my $recval=($cpvals[1]-$cpvals[0])/($cplocs[1] -$cplocs[0]); # current recombination rate per base pair value
my $recval=0; # current recombination rate per base pair value
my $tmpcdfdiff=0;
my $outputgap=(floor($numsnps/100));

my $snploc=$snps[0]; # last snp we saw
## Get to the start of the SNPs in the map
while ($snploc > $cplocs[$recinterval+1] && $recinterval+1<scalar(@cplocs)) {
    ++$recinterval;
    $recval=($cpvals[$recinterval+1]-$cpvals[$recinterval])/($cplocs[$recinterval+1] -$cplocs[$recinterval]); # Assign the first recombination rate when we encounter the map
}
 
## Iterate through the snps
for (my $i = 0; $i < $numsnps-1; ++$i){ # snp 0 is a fake snp to be ignored
  if($verbose && $numsnps>100 && ceil($i/$outputgap)*$outputgap==$i && $i>0 ){
      my $perc=100*$i/(100*floor($numsnps/100)); 
      print " $perc\%: Processing SNP $i of $numsnps\n";
  }
  $snploc=$snps[$i]; # last snp we saw
  my $nextsnploc=$snps[$i+1]; # we output the distance to the NEXT SNP
  my $lastloc=$snploc; # last rec change point OR snp we saw

  if($recinterval<0 && $snploc >= $cplocs[0]){ 
      $recinterval=0;
      $recval=($cpvals[1]-$cpvals[0])/($cplocs[1] -$cplocs[0]);
  } # Assign the first recombination rate when we encounter the map
  $cdfdiff=0;
  my $tmplen=scalar(@cplocs);
  my $tmplen2=scalar(@snps);
  while ($nextsnploc > $cplocs[$recinterval+1] && $recinterval+1<scalar(@cplocs)) {
 	my $cdjump=$recval*($cplocs[$recinterval+1]-$lastloc);
	my $tmp=$cplocs[$recinterval+1];
#      print "RECLOC: cdjump=$cdjump, lastloc=$lastloc, nextrecloc=$tmp, reciterval=$recinterval (snploc=$snploc,nextsnploc=$nextsnploc, rv=$recval)\n";
	$cdfdiff+= $cdjump;
	$recinterval++;
	$lastloc=$cplocs[$recinterval];
	if($recinterval+1<scalar(@cplocs)) {
	    $recval=($cpvals[$recinterval+1]-$cpvals[$recinterval])/($cplocs[$recinterval+1] -$cplocs[$recinterval]);
	}else {
	    $recval=$zerorecrate;
	}
  }
  if($nextsnploc > $cplocs[-1]) {
      $recval=0; # rate is zero outside the map
  }
  if($nextsnploc >= $lastloc) {
	my $cdjump=$recval*($nextsnploc -$lastloc);
 #     print "SNPLOC: cdjump=$cdjump, lastloc=$lastloc, nextsnploc=$nextsnploc  (snploc=$snploc,nextsnploc=$nextsnploc, rv=$recval)\n";
	$cdfdiff +=$cdjump;
  }
  my $overallrecrate = $cdfdiff/($nextsnploc-$snploc);
  if(($snploc < $cplocs[0] || $snploc> $cplocs[-1] )&& $overallrecrate==0) {
      print "WARNING: SNP difference $snploc-$nextsnploc is outside the map range so is being given recombination rate $zerorecrate Morgans/base.\n";
      $overallrecrate = $zerorecrate;
  }
  $tmpcdfdiff+=$cdfdiff;
  if($outputcum == 1) {
    printf OUTPUTFILE "%i %.20f\n", $snploc, $tmpcdfdiff;
  } else {
    printf OUTPUTFILE "%i %.20f\n", $snploc, $overallrecrate;
  }
#  print "Current: SNP=$snploc RATE=$overallrecrate CDF=$tmpcdfdiff\n";
  $lastloc=$snploc;
}
printf OUTPUTFILE "%i %.20f\n", $snps[-1],0.0;

close OUTPUTFILE;


