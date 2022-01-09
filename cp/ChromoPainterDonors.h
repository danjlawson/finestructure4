#ifdef __cplusplus
extern "C" {
#endif

#ifndef CHROMOPAINTERDONORS_H
#define CHROMOPAINTERDONORS_H

#include <stdlib.h>
#include <stdio.h>

#include "ChromoPainterData.h"
#include "ChromoPainterPar.h"
#include "ChromoPainterInfiles.h"
#include <stdbool.h>
  
  struct donor_t {
    // Population level information
    int ndonorpops; // number of donor populations  (in a given run)
    int nrecippops; // number of recipient populations (= number of recipient inds in all-vs-all cases)
    int ndonorvecsize; // number of donors that may be present (i.e. nrecipinds for all-vs-all mode, otherwise number of donor pops)
    char ** donorlabels; // label of each donor population
    char ** recippoplabels; // label of each recipient population
    int *recippops;/// recipient population for each recipient individual
    
    // individual level recipient information, for recipients...
    int nrecipinds; /// Number of recipient individuals in total
    char ** reciplabels;
    // individual level recipient information, and for donors...
    int ndonorhaps; // number of donors haploptypes in total
    int * ndonorhaps_vec; // number of donor hpalotypes in each populatiom
    // only populated when NOT in allvsall mode, and..
    double * ndonorprobs; /// Par->prior_donor_probs_inds==1
    double * ndonormutrates;/// Par->mutation_rate_ind==1
  };

  struct donor_t *createDonors(struct infiles_t *Infiles,struct ids_t *Ids,struct ids_t *Ids2,struct data_t *Data,struct param_t *Par);
  int getNpops(struct infiles_t *Infiles,struct param_t *Par,bool countdonors);
  void readDonorPops(struct donor_t *Donors,struct ids_t *Ids,struct ids_t *Ids2,struct infiles_t *Infiles,struct data_t *Data,struct param_t *Par);
  int validateDonorCount(int nhaps,int ndonorpops,struct param_t *Par);
  void populateDonors(struct donor_t *Donors,struct ids_t *Ids,struct ids_t *Ids2,struct data_t *Data,struct param_t *Par);
  void reallocateDataHaps(struct data_t *Data,struct donor_t *Donors, struct ids_t *Ids, struct param_t *Par);
  //char *** p_reciplabels,char *** donorlabels, int **ndonorhaps,

  //  void initializeCopyVecs(int ** p_pop_vec,double ** p_MutProb_vec,double ** p_copy_prob,double ** p_copy_probSTART, int ndonorpops, int *ndonorhaps,double *ndonorprobs, double *ndonormutrates,struct data_t *Data, struct param_t * Par);

  void clearDonors(struct donor_t *Donors, struct param_t * Par);

#endif

#ifdef __cplusplus
}
#endif
