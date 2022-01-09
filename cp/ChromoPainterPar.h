#ifdef __cplusplus
extern "C" {
#endif

  // Start the header
#ifndef CHROMOPAINTERPAR_H
#define CHROMOPAINTERPAR_H

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <time.h>
#include <string.h>

#include "ChromoPainterError.h"
  
  struct param_t {
    jmp_buf err;
    FILE *out;
    int geno_find;
    int recom_find;
    int donorlist_find;
    int outfile_find;
    int idfile_find;
    int EMiter_find;
    int geno2_find;
    int idfile2_find;
    int numsamples_find;
    int ne_find;
    int mut_find;
    int region_size_find;
    int copy_prop_em_find;
    int recom_em_find;
    int mutation_em_find;
    int mutationALL_em_find;
    int allvsall;
    int haploid_ind;
    int unlinked_ind;
    int prior_donor_probs_ind;
    int mutation_rate_ind;
    int print_file9_ind;
    int indcount_suppress_ind;
    int jitter_locations;

    int start_val; // -a start individual
    int end_val; // -a end individual
    int errormode; // How errors are reported; for the GUI

    int emloci;// Number of loci to randomly pick inside a data block 

    int startlocus;// Initial locus of data to retain
    int endlocus;// Final locus of data 
    
    int rseed; // random number seed

    int EMruns;
    int samplesTOT;
    
    int verbose; // verbose mode
    int vverbose; // very verbose mode (not recommended for general use)

    int printnorecprobs; // whether we print the probability of no recombination

    // REAL PARAMETERS START HERE
    double small_recom_val; // Smallest recombination rate allowed
    double small_copy_val; // (!!!) copy props per hap not allowed to go below this value, even if E-M wants to make them lower (!!!)
    double region_size;    // number of chunks per region -- used to look at variability in copied chunks across regions in order to estimate "c" in Dan Lawson's fineSTRUCTURE
    
    double mut_rate_self;// mutation rate for the self population (if -m)
    double GlobalMutRate;// mutation rate for all populations (if -M)
    double N_e; // Initial value of Ne
    double readN_e; // Initial value of Ne, read from cmd line


  };
    
  struct param_t *DefaultParam(); // Create default parameters
  //  struct param_t *CreateFiles(); // Create default parameters

  void setNe(int nhaps,struct param_t *Par);/// calculate N_e based on the number of haplotypes
  void parameterCheck(struct param_t *Par);/// Check the parameters

  void DestroyParam(struct param_t *Par); //Destroy the parameters



#endif
  // End the header

#ifdef __cplusplus
}
#endif
