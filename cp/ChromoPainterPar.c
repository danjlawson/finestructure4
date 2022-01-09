#include "ChromoPainterPar.h"


struct param_t *DefaultParam() {

  struct param_t *Par;
  Par = malloc(sizeof(struct param_t));
  Par->out=stdout;  
  Par->geno_find=0;
  Par->recom_find=0;
  Par->donorlist_find=0;
  Par->outfile_find=0;
  Par->idfile_find=0;
  Par->EMiter_find=0;
  Par->idfile2_find=0;
  Par->geno2_find=0;
  Par->numsamples_find=0;
  Par->ne_find=0;
  Par->mut_find=0;
  Par->region_size_find=0;
  Par->copy_prop_em_find=0;
  Par->recom_em_find=0;
  Par->mutation_em_find=0;
  Par->mutationALL_em_find=0;
  Par->allvsall=0;
  Par->haploid_ind=0;
  Par->unlinked_ind=0;
  Par->prior_donor_probs_ind=0;
  Par->mutation_rate_ind=0;
  Par->print_file9_ind=0;
  Par->indcount_suppress_ind=0;
  Par->jitter_locations=0;
  Par->mut_rate_self=-1;

  Par->start_val=0;
  Par->end_val=0;
  Par->errormode=0;
  Par->rseed=-1; // Random number seed

  Par->emloci=-1; // Number of loci to randomly place in a block
  Par->startlocus=0; // Initial locus
  Par->endlocus=0; // End locus (negative to mean "all")

  Par->verbose=0; // verbose mode 
  Par->vverbose=0;// very verbose mode (debugging only)

  Par->printnorecprobs=0; // whether we print the probability of no recombination

  Par->EMruns=0;
  Par->samplesTOT=10;

  Par->GlobalMutRate=-9;// default mutation rate is negative, meaning a default value is used
  Par->small_recom_val=0.000000000000001;    // lower limit for small genmap rates
  Par->small_copy_val=0.000000000000001; // (!!!) copy props per hap not allowed to go below this value, even if E-M wants to make them lower (!!!)
  Par->region_size = 100;    // number of chunks per region -- used to look a

  Par->N_e = 400000;    // Effective population size (which will be divided by nhaps)
  Par->readN_e = Par->N_e;    // Effective population size specified 
  return(Par);
}

void DestroyParam(struct param_t *Par){
  if(Par==NULL) return;
  if(Par->out!=stdout) fclose(Par->out);
  free(Par);
}


void parameterCheck(struct param_t *Par){

  // CHECK EVERYTHING: this should be moved to a validation function
  if (((Par->geno_find==0) || (Par->recom_find==0)) && (Par->unlinked_ind==0)) { printf("Error with command line (Each of -g and -r MUST be specified if data are linked). Exiting...\n"); stop_on_error(1,Par->errormode,Par->err);}
  if ((Par->geno_find==0) && (Par->unlinked_ind==1)) { printf("Error with command line (-g MUST be specified). Exiting...\n"); stop_on_error(1,Par->errormode,Par->err);}
  if ((Par->recom_find==1) && (Par->unlinked_ind==1)) { printf("Data specified as containing unlinked sites (-u). Ignoring supplied recombination rate file....\n");}
  if ((Par->mutation_em_find==1) && (Par->mutationALL_em_find==1))
    {
      printf("You have specified to estimate a global mutation (emission) rate and population-specific mutation (emission) rates. Please choose only one of the '-im' and '-iM' switches. Exiting...\n");
      stop_on_error(1,Par->errormode,Par->err);
    }
  if ((Par->mutation_rate_ind==1) && (Par->mut_find==1))
    {
      printf("You have provided values for both a global mutation (emission) rate ('-M') and population-specific mutation (emission) rates ('-m'). Please choose only one of the '-m' and '-M' switches. Exiting...\n");
      stop_on_error(1,Par->errormode,Par->err);
    }
  if ((Par->mutation_rate_ind==1) && (Par->mutationALL_em_find==1))
    {
      printf("You have specified to estimate a global mutation (emission) rate; will ignore population-specific mutation (emission) rates. If you wish to use donor-specific mutation rates, use the '-im' switch.\n");
    }
  if((Par->startlocus>0 || Par->endlocus>0)&&Par->emloci>0){ 
    printf("You have requested both a specific start region (-l) and a random start region (-e). Exiting due to ambiguity...\n"); stop_on_error(1,Par->errormode,Par->err);
  }
}

void setNe(int nhaps,struct param_t *Par) {
  if (Par->ne_find==0)  Par->N_e=Par->readN_e/nhaps;
}
