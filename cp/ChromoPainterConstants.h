#ifdef __cplusplus
extern "C" {
#endif

#ifndef CHROMOPAINTERCONSTANTS_H
#define CHROMOPAINTERCONSTANTS_H

  static const char *helpfilestring="to run: use './ChromoPainter' with following options:\n \
       -g <geno.filein>  (REQUIRED; no default)\n			\
       -r <recommap.filein>  (REQUIRED for LINKED MODE; no default. Specify \"-u\" to run in UNLINKED mode.)\n\
       -t <name file> (REQUIRED; no default.) File containing the names of individuals, and optionally, whether they are to be included in the analysis.\n\
       -f <donorlist.filein>  file listing breakdown of donor haps by population (required unless using -a switch)\n\
       -i <int>  number of EM iterations for estimating parameters (default=0)\n\
       -in  maximize over recombination scaling constant (N_e) using E-M\n\
       -ip  maximize over copying proportions using E-M\n\
       -im  maximize over mutation (emission) probabilities using E-M\n\
       -iM  maximize over global mutation (emission) probability using E-M\n\
       -s <int>  number of samples per recipient haplotype (default=0)\n\
       -n <double>  recombination scaling constant start-value (N_e; default=400000 divided by total number of haplotypes in <geno.filein>)\n\n\
       -p  specify to use prior copying probabilities in donor list file\n\
       -m <double>  specify to use mutation (emission) probabilities in donor list file (and provide self-copying mutation rate\n\
       -M <double>  global mutation (emission) probability (default=Li & Stephen's (2003) fixed estimate)\n\
       -k <double>  specify number of expected chunks to define a 'region' (default=100)\n\
       -j  specify that individuals are haploid\n\
       -u  specify that data are unlinked\n\
       -a <a_1> <a_2>  condition individuals a_1 through a_2 on every other individual (use '-a 0 0' to do all inds)\n\
       -d create the file '.transitionprobs.out.gz' containing the probability that there was a recombination event between each locus\n\
       -b  print-out zipped file with suffix '.copyprobsperlocus.out' containing prob each recipient copies each donor at every SNP (note: file can be quite large)\n\
       -y  do NOT print individual count numbers next to population labels in output files (only relevant if '-a' switch is used)\n\
       -o <outfile-prefix>  (default = 'geno.filein')\n\
       -J jitter SNP locations if they are invalid.  SNPs that have the same location are placed 1 basepair after the previous SNP.\n\
       -e <nloci> extract out a random block of data to run EM parameter estimation on, of length nloci.\n\
       -l <start> <end> process only loci in the range [start end), treating the rest of the data as if it did not exist. (Note that <start> begins at 0, so <end> should be the number of loci or less (zero to keep everything from start locus)\n\
       -S <seed> Set the RNG seed. Currently effects: a) the random choice of block being used by -e. If you set the same seed, you will get the same block! b) samples from the algorithm (-s option)\n\
       -Rg <geno.filein> Recipient genotype file, structured as in -g. If provided, this adds a recipient population consisting of all individuals in this file. If in all-vs-all mode (-a) ONLY INDIVIDUALS IN THIS FILE ARE RECIPIENTS and ALL INDIVIDUALS IN THE ORIGINAL GENOTYPE FILE are DONORS (so no longer all-vs-all). If in donor mode, the DONORS and RECIPIENTS in the original genotype file are respected, and these are added to the end. You do not need to specify any recipients in the donor file, this file counts as a single recipient population.\n\
       -Rt <name file> Recipient id file, structured as in -t. Only use if -Rg is given; recipients will be named (and appropriately omitted, if the inclusion column is zero) from this file. If omitted, recipients are called IND1-N.\n\
";

#endif

#ifdef __cplusplus
}
#endif 
