#!/bin/bash
help(){
    echo "Convert HGDP fastphase to phase. Usage: fastphase2phase.sh <hapguess_switch.out.gz file> <.final.gz file> <outputfile>";
}

function ceiling() {
    float_in=$1
    ceil_val=${float_in/.*}
    ceil_val=$((ceil_val+1))
}

if [ "$#" -le 2 ] ; then
    help
    exit 0
fi
if [ "$1" == "-h" ] ; then
    help
    exit 0
fi
file="$1"
locfile="$2"
out="$3"
echo Processing $locfile
locline=`zcat $locfile  | cut -f4 | tail -n +2 | tr '\n' ' '`
nsnpslocline=`echo $locline | wc -w` 
echo "Found $nsnpslocline snps in the snp locations file"
echo Converting $file to $out
tmp=`mktemp`
zcat $file | sed '/^[^0-9]/d' | sed 's/ //g' | sed '/^\s*$/d' | head -n 1876 >  $tmp
nline=`wc $tmp | tr -s ' '`
#echo "wc finds $nline"
nind=`echo $nline | cut -f 1 -d ' '`
echo "We found $nind haplotypes"
nsnpind=`echo $nline | cut -f 3 -d ' '`
echo "And a matrix size is $nsnpind"
nsnp=`echo "$nsnpind / $nind - 1" | bc`
nsnpl=`echo "$nsnpind / $nind - 1 " | bc -l`
nsnplc=${nsnpl%.*} #`ceiling $nsnpl`
#echo "0 $nsnpl A $nsnplc B $nsnpslocline"
if [ $nsnplc -ne $nsnpslocline ] ; then
    echo "ERROR: Number of snps in the locations file ($nsnpslocline) does not match that in the phased data file ($nsnps)"
    exit 1
fi
echo "Confirming $nsnp snps"
echo "$nind" > $out
echo "$nsnp" >> $out
echo "P $locline" >> $out
cat $tmp >> $out
rm $tmp
echo "Successfully created $out"