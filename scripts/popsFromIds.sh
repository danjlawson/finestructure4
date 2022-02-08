#!/bin/sh
if [ $# -ne 2 ] ; then
    echo "Usage: popsFromIds.sh ids.out ids_with_pops.txt"
    echo "where ids.txt is only the ID labels, and ids_with_pops.txt will cintain guessed populations (by removing numbers) as well as the inclusion column, for use in advanced mode in chromopainter."
    exit 0
fi

in=$1
out=$2

cut -f1 -d' ' $in | awk '{f=$0;gsub(/[0-9]/,"",f);printf("%s\t%s\t1\n",$1,f )}' > $out
