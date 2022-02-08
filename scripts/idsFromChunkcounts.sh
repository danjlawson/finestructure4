#!/bin/sh
if [ $# -ne 2 ] ; then
    echo "Usage: idsFromChunkcounts.sh chunkcounts.out ids.txt"
    echo "where chunkcounts.out is a chromopainter output file (with or without the #Cfactor line) and ids.txt will contain the ids from that file, usable as an input file with -t."
    exit 0
fi

in=$1
out=$2

cut -f1 -d' ' $in | grep -v "#Cfactor" | tail -n +2 > $out
