#!/bin/sh
# This example is on simulated data, indicative of the split between africa, europe and asia
# It was generated using the following commands, using the tool "msms" available from :

#> msms -N 10000 -t 5000 -r 4400 100000000 -ms 70 1 -I 3 30 20 20 -en 0.2 3 0.1 -ej 0.25 3 2 -en 0.25 2 0.1 -em 0.25 1 2 50 -ej 0.4 2 1 > example.txt
#> msms2cp.pl -n 10000000 example.txt example_cp 
#> makeuniformrecfile.pl example_cp.phase example_cp.recombfile

# We have 15 African individuals (IND1-15) 10 Europeans (Inds16-25) and 10 Asians (Inds 26-35)

## This is how we process this in finestructure:
fs example_cp.cp -n -phasefiles example_cp.phase -recombfiles example_cp.recombfile -idfile example_cp.ids -s1minsnps 5000 -s3iters 10000 -s4iters 10000 -go
# Takes about 4 mins on a 2nd gen Intel i5 dual core laptop.

# See the help for those commands for details:
#> fs -h s1minsnps
#Help for Parameter s1minsnps : Minimum number of SNPs for EM estimation (for chromopainter -e, default: 10000)
#> fs -h s3iters
#Help for Parameter s3iters : Number of TOTAL iterations to use for MCMC. By default we assign half to burnin and half to sampling. (default: 100000)
#$ fs -h s4iters
#Help for Parameter s4iters : Number of maximization steps when finding the best state from which the tree is built. (default: 10000)

# There isn't a lot to look at - the model is certain about the clustering, and it is right.  We'll explore more in example2.

