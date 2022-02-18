#!/bin/bash
fs="fs"
awk="awk"
sloc="../scripts/"
eloc="../examples/"
fmt="fold -w 120 -s"
cmds="info actions output parameters input stages tools s1indfrac"
tools="cp combine fs"
scripts="finestructuregreedy.sh makeuniformrecfile.pl vcf2cp.pl convertrecfile.pl phasescreen.pl phasesubsample.pl beagle2chromopainter.pl impute2chromopainter.pl msms2cp.pl"
examples="example1/example1.sh example2/example2.sh example3/example_hgdp.sh"
#transpose.pl neaverage.pl chromopainterindivrename.pl phase2beagle.pl  ped2ippca.pl  chromopainter2impute2.pl

$fs -h | $fmt > fshelp.txt

for cmd in $cmds; do
    $fs -h $cmd | $fmt  > fshelp$cmd.txt
done

for tool in $tools; do
    $fs $tool -h | grep -A 9999999 "Usage" | $fmt  > fshelp$tool.txt
done

for script in $scripts; do
    $sloc$script | $fmt  > fsscripts$script.txt
done

for example in $examples; do
    cat $eloc$example | $fmt  > fs`basename $example`.txt
done

fs -V | head -n 1 | cut -f2 -d' ' > fsversion.txt
pdflatex manual.tex && pdflatex manual.tex

### HTML version
htlatex manual.tex "xhtml,2"
cat default.css >> manual.css

mkdir -p backupdir
sidebarfile="sidebar.html"
files=`ls manual*.html`
for file in $files ; do
    lineno=`grep -n '</head><body' $file | $awk -F ':' '{print $1}'`
    lineno2=$(( $lineno - 1 ))
    echo "Processing $file (keeping up to line $lineno)"
    head -n $lineno2 $file > $file.tmp.html
    cat $sidebarfile >> $file.tmp.html
    tail -n +$lineno $file | sed 's|</body></html>|</td></tr></table></body></html>|' >> $file.tmp.html
    cp $file backupdir/$file
    mv $file.tmp.html $file
done
