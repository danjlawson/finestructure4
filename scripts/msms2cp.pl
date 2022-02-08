#!/usr/bin/perl
use Switch;
use strict;
use POSIX;

sub help {
print("CONVERTS MSMS/SCRM/MS/MSPRIME OUTPUT TO CHROMOPAINTER-STYLE INPUT FILES\n");
print("usage:   perl msms2cp.pl <options> msmsoutput.txt output_filename_prefix\n\n");
print("OPTIONS\n");
print("-c1    : Output chromopainter version1 format\n");
print("-n <x> : Multiplier for the SNP locations (default: 1000000)\n");
print("-p <x> : Specify the ploidy (default:2 for diploid; needed only for CP version 1)\n");
print("-ms <x>: Specify ms mode, and give the number of *haplotypes* in it (because ms doesn't include this in the header)\n");
print("-v     : Verbose mode\n");


die "\n";
}

my $verbose=0;
my $MSMSinfile="";
my $outfilePRE="";
my $argon=0;
my $cpversion1=0;
my $snpfactor=1000000;
my $msver="unknown";
my $nhaps=0;
my $ploidy=2;

for (my $i = 0; $i < scalar(@ARGV); ++$i){
	if(@ARGV[$i] eq "-c1"){
	    $cpversion1=1;
	}elsif(@ARGV[$i] eq "-p"){
	    $ploidy=@ARGV[++$i];
	}elsif(@ARGV[$i] eq "-n"){
	    $snpfactor=@ARGV[++$i];
	}elsif(@ARGV[$i] eq "-ms"){
	    $msver="ms";
	    $nhaps=@ARGV[++$i];
	}elsif(@ARGV[$i] eq "-msprime"){
	    $msver="msprime";
	    $nhaps=@ARGV[++$i];
	}elsif(@ARGV[$i] eq "-v"){
	    $verbose=1;
	}else{
	    switch($argon){
		case 0 {$MSMSinfile="$ARGV[$i]";}
		case 1 {$outfilePRE="$ARGV[$i]";}
		else {
		    help();
		}
	    }
	    $argon++;
	}
}


if($outfilePRE eq "" || $argon != 2) {help();}


my $nsnps=0;
my $totalINDS=0;
my @snplocs;
my @popsizes;
my $npops=0;
## PROGRAM:
## (II) GET NUMBER OF SITES AND INDS: 
open(IN,"$MSMSinfile");
open(OUT,">${outfilePRE}.phase");
my $starteddata=0;

my $lineno=0;
while(<IN>)
{
    chomp;
    my $line=$_;
    ++$lineno;
    my @linearray;
    my $start=substr($line,0,2);
    if(($lineno>1) & ($msver eq "ms" | $msver eq "unknown")){
	if($lineno==2){
	    if($msver ne "ms") {
		die "Detected MS output, but this is not specified. You must use \"-ms <n>\" to specify the number of haplotypes for ms output\n";
	    }
	    if($verbose) {print "Detected MS output in MS mode\n";}
	    $start="se";
	    $line="segsites: $line";
	    $totalINDS=$nhaps/$ploidy;
	    if($verbose){print "Informed of $nhaps haplotypes ($totalINDS inds for CPv1 only)\n";}
	}elsif($lineno==3){
	    $start="po";
	    $starteddata=1;
	    $line="positions: $line";
	}
      }
    
    if($start eq "ms" || $start eq "sc" || index($line, "macs")!=-1 || index($line, "mspms") != -1){
      if($msver eq "ms"){
	$msver=substr($line,0,4);
	die "Detected $msver output, but you have specified MS. Rerun omitting the \"-ms\" flag.";
      }
      if(index($line, "mspms") != -1){
	$msver="mspms";
      }elsif(index($line, "macs") != -1){
	$msver="macs";
      }else{
	$msver=substr($line,0,4);
      }
      if($verbose) {print "Detected $msver output\n";}
      @linearray=split(/\s+/,$line);
      my $I=-1;
      for(my $i=0;$i<scalar(@linearray);++$i){
	if($linearray[$i] eq "-I"){ $I=$i;}
      }
      for(my $i=$I+2;$i<$I+2+$linearray[$I+1];++$i){
	$nhaps+=$linearray[$i];
	push (@popsizes,$linearray[$i]);
      }
      $totalINDS=$nhaps/$ploidy;
      if($verbose){print "Detected $nhaps haplotypes ($totalINDS inds for CPv1 only)\n";}
    }elsif($start eq "se"){
      @linearray=split(/\s+/,$line);
      $nsnps=$linearray[1];
      if($verbose) {print "Detected $nsnps SNPs\n";}
    }elsif($start eq "po"){
      @snplocs=split(/\s+/,$line);
      shift(@snplocs);
      for(my $i=0;$i<scalar(@snplocs);++$i){
	$snplocs[$i]=floor($snplocs[$i]*$snpfactor);
	if($i>0 && $snplocs[$i] <= $snplocs[$i-1]){
	  $snplocs[$i]=$snplocs[$i-1]+1;
	}
	#	    print "SNP $i location $snplocs[$i]\n";
      }
      if($cpversion1){
	print OUT "0\n";
	print OUT "$totalINDS\n";
      }else{
	print OUT "$nhaps\n";
      }
      print OUT "$nsnps\n";
      print OUT "P @snplocs\n";
      if($cpversion1){
	for (my $j=0; $j < $nsnps; $j+=1)
	  {
	    print OUT "S";
	  }
	print OUT "\n";
      }	
      $starteddata=1;
    }elsif(($start eq "//")|| ($start eq "")){
	# Do nothing
	if($verbose) {print "Ignoring line:$line\n";}
    }elsif($starteddata){
	print OUT "$line\n";
	
    }
}
close(IN);
close(OUT);

$npops=scalar(@popsizes);
print "Making HAPLOID ID file ${outfilePRE}.hap.ids with $npops populations\n";

open(OUT,">${outfilePRE}.hap.ids");
my $dippossible=1;
my $hapid=1;
for (my $i=0; $i < $npops;++$i)
  {
    if($popsizes[$i] % 2 !=0) {$dippossible=0;}
    my $popid=$i+1;
    for (my $j=0; $j < $popsizes[$i]; $j+=1)
      {
	print OUT "HAP$hapid POP$popid 1\n";
	++$hapid;
      }
  }
close(OUT);

if($dippossible==1){
  print "Making DIPLOID ID file ${outfilePRE}.ids with $npops populations\n";
  open(OUT,">${outfilePRE}.ids");
  my $indid=1;
  for (my $i=0; $i < $npops;++$i)
    {
      my $popid=$i+1;
      for (my $j=0; $j < $popsizes[$i]/2; $j+=1)
	{
	  print OUT "IND$indid POP$popid 1\n";
	  ++$indid;
	}
    }
  close(OUT);
}else{
  print "Data are not compatible with a diploid population, so omitting creating this ID file.\n"
}


if($msver eq "ms"){
    my $obshaps=$lineno-3;
    if($nhaps != $obshaps){
	print "WARNING! $nhaps haplotypes were specified, but the data appears to contain $obshaps haplotypes. You can fix this manually in the output file, or you can rerun specifying the right number of haplotypes!\n";
    }else{
	print "Information: The number of haplotypes for ms ($nhaps) appears to have been correctly specified.\n";
    }
}

