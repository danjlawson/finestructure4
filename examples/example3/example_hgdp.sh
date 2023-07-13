#!/bin/sh
## This example is intended to be run on a HPC machine. In particular, it will only work verbatim on a torque-based qsub system.  Naturally, you can adapt it to run on other systems. The only part to change is "qsub_run.sh"



####### Download the sample ID info
wget -nc http://hgdp.uchicago.edu/Phased_data/H938_Clumps.clst.txt_check.Continent.format.order.NR.gz
zcat H938_Clumps.clst.txt_check.Continent.format.order.NR.gz  | cut -f2 -d' ' | awk 'NR%2==0{printf("%s%i %s 1\n",$1,NR/2,$1);}' > hgdp.ids

####### Download the genetic map
wget -nc http://hapmap.ncbi.nlm.nih.gov/downloads/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz
tar -xzvf genetic_map_HapMapII_GRCh37.tar.gz

####### Create scripts that download the raw data and convert it into chromopainter input format
mkdir -p getdata
cmdfile="getdata_raw.sh"
rm -f $cmdfile
for chr in `seq 21 22`; do
    cmdf="getdata/getdata_chr$chr.sh"
    echo "#!/bin/sh" > $cmdf
    echo "wget -nc http://hgdp.uchicago.edu/Phased_data/chrom${chr}_hapguess_switch.out.gz" >> $cmdf
    echo "wget -nc http://hgdp.uchicago.edu/Phased_data/chrom${chr}.final.gz" >> $cmdf
    echo "fastphase2cp.sh chrom${chr}_hapguess_switch.out.gz chrom${chr}.final.gz chrom${chr}.phase" >> $cmdf
    echo "convertrecfile.pl -M hapmap chrom${chr}.phase genetic_map_GRCh37_chr${chr}.txt chrom${chr}.recombfile" >> $cmdf
    chmod +x $cmdf
    echo "$cmdf" >> $cmdfile
done

## If you want the list of RSIDs
zcat chrom*.final.gz | tail -n +2 | cut -f1 > hgdprsids.txt

######## Submit these scripts to the qsub system
qsub_run.sh -n 8 -f getdata_raw.sh -w 1
## (nb: you could use cat getdata_raw.sh | parallel if you are using a single machine with many cores)

## Now we start the finestucture run proper.  First we do the painting.
fs hgdp_2chr.cp -n -phasefiles chrom21.phase chrom22.phase -recombfiles chrom21.recombfile chrom22.recombfile -idfile hgdp.ids -hpc 1 -s1indfrac 0.1 -go # We are just using chr21-22.
qsub_run.sh -f hgdp_2chr/commandfiles/commandfile1.txt -n 8 -m 16 # choose -m to get the number of jobs we want. There are 186 commands to run, so this gives us 12 jobs. Don't go above a couple of hundred.
fs hgdp_2chr.cp -indsperproc 100 -go # Set indsperproc again to control the amount of jobs. This is a short-term hack because the default (1) produces too many files - 22 * 938 * 9 = 185,724 which breaks filesystem lookups
qsub_run.sh -f hgdp_2chr/commandfiles/commandfile2.txt -n 8 # We're already doing 100 inds per node requested. So there are only 20 commands this time, and we don't need -m to run them serially
 ## Finestructure section, no more painting.  So things scale as N^3 with no dependence on L.
fs hgdp_2chr.cp -s3iters 1000000 -go
qsub_run.sh -f hgdp_2chr/commandfiles/commandfile3.txt -n 2
## But we received a the following messages:
#WARNING: Failed Gelman-Rubin MCMC diagnostic check with maximum potential scale reduction factor 1.39013 (threshold 1.3)
#WARNING: Stage 3 convergence criterion failed! Use "-ignoreGR" to continue regardless. (Set the parameter "-threshGR:-1" to ignore the GR statistic in future.)  Re-running MCMC for longer: this is attempt 1 of 5.
## So we rerun, with automatically doubled duration:
fs hgdp_2chr.cp -go
qsub_run.sh -f hgdp_2chr/commandfiles/commandfile3.txt -n 2
## ... And again...
fs hgdp_2chr.cp -go
qsub_run.sh -f hgdp_2chr/commandfiles/commandfile3.txt -n 2
## ... And again...
fs hgdp_2chr.cp -go
qsub_run.sh -f hgdp_2chr/commandfiles/commandfile3.txt -n 2
### Finally, convergence is achieved. We can run the tree. This is quicker
fs hgdp_2chr.cp -go
qsub_run.sh -f hgdp_2chr/commandfiles/commandfile4.txt -n 2
## And we are done:
fs hgdp_2chr.cp -go

