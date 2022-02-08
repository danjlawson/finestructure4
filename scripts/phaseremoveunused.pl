#!/usr/bin/perl
## 
use Switch;
use strict;
use POSIX;

sub help {
print("REMOVE UNUSED INDIVIDUALS FROM A PHASE FILE, MAKING A NEW VALID PHASE AND ID FILE\n");
print("usage:   ids2ind.pl <options> in.ids out.ids in.phase out.phase\n\n");
print("OPTIONS\n");
print("-v     : Verbose mode\n");

die "\n";
}

sub trim($){  # remove whitespace from beginning and end of the argument
        my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
}

my $verbose=0;
my $argon=0;

for (my $i = 0; $i < scalar(@ARGV); ++$i){
  if(@ARGV[$i] eq "-v"){
    $verbose=1;
    $argon++;
  }
}

if(scalar(@ARGV) - $argon != 4){
  help();
}

my $idfile=$ARGV[$argon++];
my $idfileout=$ARGV[$argon++];
my $infile=$ARGV[$argon++];
my $outfile=$ARGV[$argon++];

open(ID,"$idfile");
open(IDOUT,"$idfileout");

my @retain;
my $nids=0;

while(<ID>)
{
  chomp;
  my $line=$_;
  my @linearray=split(/\s+/,$line);
  if(scalar(@linearray)!=3) {
    die("Expect standard ChromoPainter id format with 3 columns: id pop inclusion\n");
  }
  push @retain, $linearray[2];
  if ($linearray[2] ==0 && $verbose){
    print "Removing individual $linearray[0]\n";
  }else{
    print IDOUT "$line\n";
    ++$nids;
  }
}
close(ID);
close(IDOUT);
if($verbose){ print "Wrote $nids ID lines to file $idfileout\n"; }

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

my @hretain;
my $ploidy=1;
if($ninds==scalar(@retain)){
  print "Detected haploid data\n";
}elsif($ninds==2*scalar(@retain)){
  print "Detected diploid data\n";
  $ploidy=2;
}else{
  $ploidy=$ninds/scalar(@retain);
  die "INVALID PLOIDY: $ploidy\n";
}
$nids*=$ploidy;

for(my $i=0;$i<scalar(@retain);++$i){
  for(my $j=0;$j<$ploidy;++$j){
    push @hretain, $retain[$i];
  }
}

print OUT "$nids\n";
print OUT "$nsnps\n";
print OUT "P @posvec\n";

my $nhaps=0;
while(<IN>)
  {
    my $line=trim($_);
    if($hretain[$nhaps]==1){
      print OUT "$line\n";
    }
    $nhaps++;
  }
close(IN);

close(OUT);
