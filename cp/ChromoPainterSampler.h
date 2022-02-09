#ifdef __cplusplus
extern "C" {
#endif

#ifndef CHROMOPAINTERSAMPLER_H
#define CHROMOPAINTERSAMPLER_H


  
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <time.h>
#include <string.h>

#include "zlib.h"
#include "ChromoPainterPar.h"
#include "ChromoPainterData.h"
#include "ChromoPainterDonors.h"
#include "ChromoPainterOutfiles.h"
#include "ChromoPainterError.h"
#include "ChromoPainterRecmap.h"

  /*  struct fb_t {
    int * newh;
    int ** existing_h;
    double ** Alphamat;
    double * MutProb_vec;
    int *p_Nhaps;
    int *p_Nloci;
    double * copy_prob;
    double * copy_probSTART;
    double * TransProb;
  };  */
 ///* A data type in which all the forward backward details are stored

  ///
  double thetaLiStephens(double * MutProb_vec, int *p_nchr,int *p_Nhaps);// calculate theta, as given in Li and Stephens. Additionally, set the mutation vector to this value if mutprob<0

  double deltaLiStephens(double * TransProb, double * pos, double p_rhobar, double * lambda, int *p_Nloci,struct param_t *Par);/// Calculate the effective recombination rate between loci

  /// FORWARD ALGORITHM
  //  double InitialiseForward(struct fb_t *Fb); // Initialise the forward algorithm
  //  double forwardAlgorithm(struct fb_t *Fb, struct param_t *Par);  // Perform the forward algorithm inference

  double InitialiseForward(int * newh, int ** existing_h, double ** Alphamat, double * MutProb_vec, int *p_Nhaps,int *p_Nloci, double * copy_probSTART, double * TransProb); // Initialise the forward algorithm

  double forwardAlgorithm(int * newh, int ** existing_h, double ** Alphamat, double * MutProb_vec, int *p_Nhaps,int *p_Nloci,double * copy_prob, double * copy_probSTART, double * TransProb, struct param_t *Par);  // Perform the forward algorithm inference

  void  backwardAlgorithm(int finalrun,int ndonorpops,int ind_val,double Alphasum,double p_rhobar, double * N_e_new,int * newh, int ** existing_h, double ** Alphamat, double * lambda, double delta,double * MutProb_vec, int *p_Nhaps,int *p_Nloci,double * copy_prob,double * copy_prob_new,double * copy_prob_newSTART, double *corrected_chunk_count, double *expected_chunk_length, double * expected_differences,double *regional_chunk_count_sum_final,double *regional_chunk_count_sum_squared_final, int *num_regions, double * copy_probSTART, double * TransProb,int * pop_vec,double *pos, double * snp_info_measure, struct files_t *Outfiles, struct param_t *Par);///Perform the backward algorithm, including printing some stuff

  /// Additional things
  double ** sampler(double ** copy_prob_new_mat, int * newh, int ** existing_h, int *p_Nloci, int *p_Nhaps, int *p_nchr,  double p_rhobar, double * MutProb_vec, int * allelic_type_count_vec, double * lambda, double * pos, double * copy_prob, double * copy_probSTART, int * pop_vec,int * cond_mat_haplotypes, int ndonorpops, int run_num, int ind_val, struct files_t *Outfiles,struct param_t *Par); // compute the likilihood in its entirity

  int loglik(struct copyvec_t *Copyvec, struct donor_t *Donors, struct data_t *Data, struct ids_t *Ids,struct infiles_t *Infiles, struct files_t *Outfiles, struct param_t * Par);
  //  int loglik(int nhaps_startpop, int *p_nloci, int p_nhaps, double N_e_start, double * recom_map, double * MutProb_vec, double * copy_prob, double * copy_probSTART, int * pop_vec, char *filename, struct files_t *Outfiles, struct param_t * Par); // do em,  and everything else too


#endif

#ifdef __cplusplus
}
#endif



