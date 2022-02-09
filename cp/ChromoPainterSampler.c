#include "ChromoPainterSampler.h"

#include <time.h>

#define MIN_THETA 1e-8
#define MIN_NE 1e-8
#define SMALL_NUM 1e-20

///////////////////////////////////////
///////////////////////////////////////
// Sampler

double thetaLiStephens(double * MutProb_vec, int *p_nchr,int *p_Nhaps){
  int i;
  double sum = 0;
  for(i = 1; i < *p_nchr; i++){
    sum = sum + 1.0/i;
  }
  double Theta = 1.0 / sum;
  for (i=0; i < *p_Nhaps; i++)
    {
      if (MutProb_vec[i]<0) MutProb_vec[i]=0.5 * Theta/(*p_Nhaps + Theta);
      //if (MutProb_vec[i]<rounding_val) MutProb_vec[i]=rounding_val;
    }
  return(Theta);
}

double deltaLiStephens(double * TransProb, double * pos, double p_rhobar, double * lambda, int *p_Nloci,struct param_t *Par){
  double delta = 1.0;
  int locus;

  if (Par->unlinked_ind==0 && lambda[0] >= 0) TransProb[0] = 1 - exp(-1 * (pos[1]-pos[0]) * delta * p_rhobar*lambda[0]);
  if (Par->unlinked_ind==1 || lambda[0] < 0) TransProb[0] = 1.0;

  for(locus = 1; locus < *p_Nloci - 1; locus++)
    {
      delta = 1.0;
      if (Par->unlinked_ind==0 && lambda[locus] >= 0) TransProb[locus] = 1 - exp(-1 * (pos[locus+1]-pos[locus]) * delta * p_rhobar*lambda[locus]);
      if (Par->unlinked_ind==1 || lambda[locus] < 0) TransProb[locus] = 1.0;
    }
  return(delta);
}

double InitialiseForward(int * newh, int ** existing_h, double ** Alphamat, double * MutProb_vec, int *p_Nhaps,int *p_Nloci, double * copy_probSTART, double * TransProb) {
//double InitialiseForward(struct_t *Fb, double * TransProb){
  int i;
  double Alphasum=0;
  double ObsStateProb;

  // Initial values
  for (i=0; i < *p_Nhaps; i++)
    {
      if(newh[0]==9) {
	ObsStateProb=1.0;
      }else if(newh[0]==8) {
	ObsStateProb = (1-SMALL_NUM) * (newh[0] == existing_h[i][0]) + SMALL_NUM * (newh[0] != existing_h[i][0]);
      }else{
	ObsStateProb = (1-MutProb_vec[i]) * (newh[0] == existing_h[i][0]) + MutProb_vec[i] * (newh[0] != existing_h[i][0]);
      }
      
      Alphamat[i][0] = log(copy_probSTART[i]*ObsStateProb);
      Alphasum = Alphasum + exp(Alphamat[i][0])*TransProb[0];
    }
  Alphasum=log(Alphasum);
  return(Alphasum);
}

double forwardAlgorithm(int * newh, int ** existing_h, double ** Alphamat, double * MutProb_vec, int *p_Nhaps,int *p_Nloci,double * copy_prob, double * copy_probSTART, double * TransProb, struct param_t *Par) {
//double forwardAlgorithm(struct_t *Fb, struct param_t *Par){
  // Perform the forward step
  // return alphasum, the sum of the forward weightings (logged)
  if(Par->vverbose) fprintf(Par->out,"        sampler: forwards algorithm\n");
      /* FORWARDS ALGORITHM: (Rabiner 1989, p.262) */
      /* INITIALIZATION: */
  int locus;  // loci index
  int i; // haplotype index
  double Alphasumnew; // alpha for the next locus
  double ObsStateProb;
  double large_num;

  if(Par->vverbose) fprintf(Par->out,"        forward Algorithm (initialising) \n");
  double Alphasum = InitialiseForward(newh, existing_h, Alphamat, MutProb_vec,p_Nhaps,p_Nloci,copy_probSTART,TransProb);

  if(Par->vverbose) fprintf(Par->out,"        forward Algorithm (computing) \n");

  // Perform the forward pass
  for (locus=1; locus < *p_Nloci; locus++)
    {
      Alphasumnew = 0.0;
      large_num = -1.0*Alphasum;
      for (i=0; i < *p_Nhaps; i++)
	{
	  if(newh[locus]==9) {
	    ObsStateProb=1.0;
	  }else if(newh[locus]==8) {
	    ObsStateProb = (1-SMALL_NUM) * (newh[locus] == existing_h[i][locus]) + SMALL_NUM * (newh[locus] != existing_h[i][locus]);
	  }else{
	    ObsStateProb = (1-MutProb_vec[i]) * (newh[locus] == existing_h[i][locus]) + MutProb_vec[i] * (newh[locus] != existing_h[i][locus]);
	  }

	  Alphamat[i][locus] = log(ObsStateProb*copy_prob[i]*exp(Alphasum+large_num) + ObsStateProb*(1-TransProb[(locus-1)])*exp(Alphamat[i][(locus-1)]+large_num)) - large_num;
	  if (locus < (*p_Nloci - 1)) Alphasumnew = Alphasumnew + exp(Alphamat[i][locus]+large_num)*TransProb[locus];
	  if (locus == (*p_Nloci - 1)) Alphasumnew = Alphasumnew + exp(Alphamat[i][locus]+large_num);
	}
      Alphasum = log(Alphasumnew)-large_num;
    }

  // Check that all is well

  if (isnan(Alphasum))
    {
      fprintf(Par->out,"Sampler::forwardAlgorithm error: Negative or NaN likelihood. Could be because emission or transition probabilities are too low??...Exiting...\n");
      /* fprintf(Par->out,"Alpha matrix for debugging\n"); */
      /* for (locus=1; locus < *p_Nloci; locus++) */
      /* 	{ */
      /* 	  fprintf(Par->out,"LOCUS %i ",locus); */
      /* 	  for (i=0; i < *p_Nhaps; i++) */
      /* 	    { */
      /* 	      fprintf(Par->out," %f",Alphamat[i][locus]); */
      /* 	    } */
      /* 	  fprintf(Par->out,"\n"); */
      /* 	} */
      //      fprintf(Par->out,"HAP COPYPROBSTART COPYPROB log(COPYPROB)\n");
      //for (i=0; i < *p_Nhaps; i++) fprintf(Par->out,"%d %lf %lf %lf\n",i,copy_probSTART[i],copy_prob[i],log(copy_prob[i]));
      stop_on_error(1,Par->errormode,Par->err);
    }
  if(Par->vverbose) fprintf(Par->out,"        forward Algorithm (complete) \n");
  return(Alphasum);
}

///////////////////////////////////////////////
// Backwards algorithm

void  backwardAlgorithm(int finalrun,int ndonorpops,int ind_val,double Alphasum,double p_rhobar, double * N_e_new,int * newh, int ** existing_h, double ** Alphamat, double * lambda, double delta,double * MutProb_vec, int *p_Nhaps,int *p_Nloci,double * copy_prob,double * copy_prob_new,double * copy_prob_newSTART, double *corrected_chunk_count, double *expected_chunk_length, double * expected_differences,double *regional_chunk_count_sum_final,double *regional_chunk_count_sum_squared_final, int *num_regions, double * copy_probSTART, double * TransProb,int * pop_vec,double *pos, double * snp_info_measure, struct files_t *Outfiles, struct param_t *Par){

  double total_regional_chunk_count,total_gen_dist;
  double Betasum, Betasumnew;
  double large_num;
  int i,locus;
  double total_prob,total_prob_from_i_to_i,total_prob_to_i_exclude_i;
  double total_prob_from_i_exclude_i,total_prob_from_any_to_any_exclude_i;
  double constant_exclude_i,constant_from_i_to_i,constant_exclude_i_both_sides;
  
  double ObsStateProb, ObsStateProbPREV;

  double total_ind_sum;
  double expected_chunk_length_sum;
  double sum_prob;
  double * exp_copy_pop=malloc(ndonorpops * sizeof(double));
  double * BetavecPREV = malloc(*p_Nhaps * sizeof(double));
  double * BetavecCURRENT = malloc(*p_Nhaps * sizeof(double));
  double * expected_transition_prob = malloc((*p_Nloci-1)*sizeof(double));
  double * ind_snp_sum_vec = malloc(ndonorpops * sizeof(double));

  double * regional_chunk_count = malloc(*p_Nhaps * sizeof(double));
  double * regional_chunk_count_sum = malloc(ndonorpops * sizeof(double));

  double rounding_val = 1.0/10000000.0;  // for regional_counts; c is a bit lame
  
  for (i=0; i < ndonorpops; i++)
    {
      regional_chunk_count_sum[i] = 0.0;
      ind_snp_sum_vec[i]=0.0;
    }
  for (i = 0; i < *p_Nhaps; i++)
    {
      copy_prob_new[i] = 0.0;
      corrected_chunk_count[i] = 0.0;
      expected_chunk_length[i] = 0.0;
      expected_differences[i] = 0.0;
      regional_chunk_count[i] = 0.0;
    }
  total_regional_chunk_count=0.0;
  *num_regions=0;

  /* BACKWARDS ALGORITHM: (Rabiner 1989, p.263) */
  /* INITIALIZATION: */
  Betasum = 0.0;
  if (finalrun)
    {
      for (i=0; i < ndonorpops; i++)
	exp_copy_pop[i]=0.0;
    }

  for(i=0; i < *p_Nhaps; i++)
    {
      
	  if(newh[locus]==9) {
	    ObsStateProb=1.0;
	  }else if(newh[locus]==8) {
	    ObsStateProb = (1-SMALL_NUM) * (newh[(*p_Nloci-1)] == existing_h[i][(*p_Nloci-1)]) + SMALL_NUM * (newh[(*p_Nloci-1)] != existing_h[i][(*p_Nloci-1)]);
	  }else{
	    ObsStateProb = (1-MutProb_vec[i]) * (newh[(*p_Nloci-1)] == existing_h[i][(*p_Nloci-1)]) + MutProb_vec[i] * (newh[(*p_Nloci-1)] != existing_h[i][(*p_Nloci-1)]);
	  }

      BetavecPREV[i] = 0.0;
      Betasum = Betasum + TransProb[(*p_Nloci-2)]*copy_prob[i]*ObsStateProb*exp(BetavecPREV[i]);
      if (finalrun) exp_copy_pop[pop_vec[i]]=exp_copy_pop[pop_vec[i]]+exp(BetavecPREV[i]+Alphamat[i][(*p_Nloci-1)]-Alphasum);

      // for estimating new mutation rates:
      expected_differences[i]=expected_differences[i]+exp(Alphamat[i][(*p_Nloci-1)]-Alphasum)*(newh[(*p_Nloci-1)] != existing_h[i][(*p_Nloci-1)]);
    }
  if (finalrun)  printCopyProbs(exp_copy_pop,ind_val,pos[*p_Nloci-1],Outfiles,Par);


           /* INDUCTION: */
  if(Par->vverbose) fprintf(Par->out,"        sampler: induction\n");
  Betasum = log(Betasum);
  /* CALCULATE EXPECTED NUMBER OF TIMES OF COPYING TO EACH DONOR POP (Rabiner 1989, p.263,265 or Scheet/Stephens 2006 Appendix C): */
  for (locus = (*p_Nloci-2); locus >= 0; locus--)
    {
      Betasumnew = 0.0;
      large_num = -1.0*Betasum;
      total_prob=0.0;
      
      constant_exclude_i = 0.5;
      constant_from_i_to_i = 1.0;
      constant_exclude_i_both_sides = 0.0;
      expected_chunk_length_sum=0.0;
      sum_prob=0.0;
      if (finalrun)
	{
	  for (i=0; i < ndonorpops; i++)
	    exp_copy_pop[i]=0.0;
	}
      for (i = 0; i < *p_Nhaps; i++)
	{
	  if(newh[locus]==9) {
	    ObsStateProb=1.0;
	  }else if(newh[locus]==8) {
	    ObsStateProb = (1-SMALL_NUM) * (newh[locus] == existing_h[i][locus]) + SMALL_NUM * (newh[locus] != existing_h[i][locus]);
	  }else{
	    ObsStateProb = (1-MutProb_vec[i]) * (newh[locus] == existing_h[i][locus]) + MutProb_vec[i] * (newh[locus] != existing_h[i][locus]);
	  }

	  if(newh[locus+1]==9) {
	    ObsStateProbPREV=1.0;
	  }else if(newh[locus+1]==8) {
	    ObsStateProbPREV = (1-SMALL_NUM) * (newh[(locus+1)] == existing_h[i][(locus+1)]) + SMALL_NUM * (newh[(locus+1)] != existing_h[i][(locus+1)]);
	  }else {
	    ObsStateProbPREV = (1-MutProb_vec[i]) * (newh[(locus+1)] == existing_h[i][(locus+1)]) + MutProb_vec[i] * (newh[(locus+1)] != existing_h[i][(locus+1)]);
	  }
	  
	  BetavecCURRENT[i] = log(exp(Betasum+large_num) + (1-TransProb[locus]) * ObsStateProbPREV*exp(BetavecPREV[i] + large_num)) - large_num;
	  if (locus > 0) Betasumnew = Betasumnew + TransProb[(locus-1)]*copy_prob[i]*ObsStateProb*exp(BetavecCURRENT[i] + large_num);
	  if (locus == 0) copy_prob_newSTART[i] = exp(Alphamat[i][0] + BetavecCURRENT[i] - Alphasum);
	  total_prob = total_prob + exp(Alphamat[i][(locus+1)]+BetavecPREV[i]-Alphasum)-exp(Alphamat[i][locus]+BetavecPREV[i]- Alphasum)*ObsStateProbPREV*(1-TransProb[locus]);

	  copy_prob_new[i] = copy_prob_new[i] + exp(Alphamat[i][(locus+1)]+BetavecPREV[i]-Alphasum)-exp(Alphamat[i][locus]+BetavecPREV[i]- Alphasum)*ObsStateProbPREV*(1-TransProb[locus]);

	  total_prob_from_i_to_i = exp(Alphamat[i][locus]+BetavecPREV[i]-Alphasum)*ObsStateProbPREV*(1-TransProb[locus]+TransProb[locus]*copy_prob[i]);
	  total_prob_to_i_exclude_i = exp(Alphamat[i][(locus+1)]+BetavecPREV[i]-Alphasum)-exp(Alphamat[i][locus]+BetavecPREV[i]- Alphasum)*ObsStateProbPREV*(1-TransProb[locus]+TransProb[locus]*copy_prob[i]);
	  total_prob_from_i_exclude_i = exp(Alphamat[i][locus]+BetavecCURRENT[i]-Alphasum) - exp(Alphamat[i][locus]+BetavecPREV[i]-Alphasum)*ObsStateProbPREV*(1-TransProb[locus]+TransProb[locus]*copy_prob[i]);
	  total_prob_from_any_to_any_exclude_i = 1.0-exp(Alphamat[i][locus]+BetavecCURRENT[i]-Alphasum)-exp(Alphamat[i][(locus+1)]+BetavecPREV[i]-Alphasum)+exp(Alphamat[i][locus]+BetavecPREV[i]-Alphasum)*ObsStateProbPREV*(1-TransProb[locus]+TransProb[locus]*copy_prob[i]);
	  
	  regional_chunk_count[i]=regional_chunk_count[i]+(exp(Alphamat[i][(locus+1)]+BetavecPREV[i]-Alphasum)-exp(Alphamat[i][locus]+BetavecPREV[i]- Alphasum)*ObsStateProbPREV*(1-TransProb[locus]));
	  total_regional_chunk_count=total_regional_chunk_count+(exp(Alphamat[i][(locus+1)]+BetavecPREV[i]-Alphasum)-exp(Alphamat[i][locus]+BetavecPREV[i]- Alphasum)*ObsStateProbPREV*(1-TransProb[locus]));
	  ind_snp_sum_vec[pop_vec[i]]=ind_snp_sum_vec[pop_vec[i]]+(exp(Alphamat[i][(locus+1)]+BetavecPREV[i]-Alphasum)-exp(Alphamat[i][locus]+BetavecPREV[i]- Alphasum)*ObsStateProbPREV*(1-TransProb[locus]));
	  
	  corrected_chunk_count[i]=corrected_chunk_count[i]+(exp(Alphamat[i][(locus+1)]+BetavecPREV[i]-Alphasum)-exp(Alphamat[i][locus]+BetavecPREV[i]- Alphasum)*ObsStateProbPREV*(1-TransProb[locus]));
	  if (Par->unlinked_ind==0 && lambda[locus]>=0) expected_chunk_length[i]=expected_chunk_length[i]+100*(pos[locus+1]-pos[locus])*delta*lambda[locus]*(constant_from_i_to_i*total_prob_from_i_to_i+constant_exclude_i*(total_prob_to_i_exclude_i+total_prob_from_i_exclude_i)+constant_exclude_i_both_sides*total_prob_from_any_to_any_exclude_i);  // multiply by 100 to get cM
	  expected_chunk_length_sum=expected_chunk_length_sum+constant_from_i_to_i*total_prob_from_i_to_i+constant_exclude_i*(total_prob_to_i_exclude_i+total_prob_from_i_exclude_i)+constant_exclude_i_both_sides*total_prob_from_any_to_any_exclude_i;

	  // for estimating new mutation rates:
	  expected_differences[i]=expected_differences[i]+exp(Alphamat[i][locus]+BetavecCURRENT[i]-Alphasum)*(newh[locus] != existing_h[i][locus]);
	  BetavecPREV[i] = BetavecCURRENT[i];

	  if (finalrun) exp_copy_pop[pop_vec[i]]=exp_copy_pop[pop_vec[i]]+exp(BetavecCURRENT[i]+Alphamat[i][locus]-Alphasum);

	  sum_prob=sum_prob+total_prob_from_i_to_i+total_prob_to_i_exclude_i+total_prob_from_i_exclude_i;

	  // for calculating the total probability of a recombination between locus and i
	  //fprintf(Par->out,"TPFL[%i]=%f\n",i,total_prob_from_i_to_i);
	  //	  probfromlocus[i]<-PREVprobfromlocus[i]*total_prob_from_i_to_i;
	} // end for every haplotype
      
      if (finalrun) printCopyProbs(exp_copy_pop,ind_val,pos[locus],Outfiles,Par);

      expected_transition_prob[locus]=total_prob;


      if (locus > 0) Betasum = log(Betasumnew) - large_num;

      if ((total_regional_chunk_count+rounding_val) >= Par->region_size)
	{// a new region
	  for (i = 0; i < *p_Nhaps; i++)
	    {
	      //fprintf(Par->out,"%d %d %lf %lf\n",i,pop_vec[i],regional_chunk_count[i],regional_chunk_count_sum[pop_vec[i]]);
	      regional_chunk_count_sum[pop_vec[i]]=regional_chunk_count_sum[pop_vec[i]]+regional_chunk_count[i];
	      regional_chunk_count[i]=0.0;
	    }
	  for (i = 0; i < ndonorpops; i++)
	    {
	      regional_chunk_count_sum_final[i]=regional_chunk_count_sum_final[i]+regional_chunk_count_sum[i];
	      regional_chunk_count_sum_squared_final[i]=regional_chunk_count_sum_squared_final[i]+pow(regional_chunk_count_sum[i],2.0);
	      regional_chunk_count_sum[i]=0.0;
	    }
	  total_regional_chunk_count=0.0;
	  *num_regions=*num_regions+1;
	} // end  if end of a region

      total_ind_sum=0.0;
      for (i = 0; i < ndonorpops; i++)
	total_ind_sum=total_ind_sum+ind_snp_sum_vec[i];
      for (i = 0; i < ndonorpops; i++)
	{
	  snp_info_measure[i]=snp_info_measure[i]+pow((ind_snp_sum_vec[i]/total_ind_sum),2.0);
	  ind_snp_sum_vec[i]=0.0;
	}
    } // end loop over loci

  // print the total recombination rates
  if(finalrun) printTransitionProb(expected_transition_prob,ind_val,(int)(*p_Nloci),Outfiles);

  for (i=0; i < ndonorpops; i++)
    snp_info_measure[i]=snp_info_measure[i]/(*p_Nloci);

  if(Par->vverbose) fprintf(Par->out,"        sampler: compute expectations\n");

  /* CALCULATE EXPECTED NUMBER OF TOTAL TRANSITIONS, IN ORDER TO ESTIMATE N_e (Scheet/Stephens 2006 Appendix C (C3)): */
  total_prob=0.0;
  total_gen_dist=0.0;
  for (locus = 0; locus < (*p_Nloci-1); locus++)
    {
      if (Par->unlinked_ind==0 && lambda[locus] >= 0) total_gen_dist=total_gen_dist+(pos[(locus+1)]-pos[locus])*delta*lambda[locus];
      if (Par->unlinked_ind==0 && lambda[locus] >= 0) total_prob=total_prob+((p_rhobar*(pos[(locus+1)]-pos[locus])*delta*lambda[locus])/(1.0-exp(-1.0*p_rhobar*(pos[(locus+1)]-pos[locus])*delta*lambda[locus])))*expected_transition_prob[locus];
    }
  if (Par->unlinked_ind==0){
    *N_e_new = total_prob/total_gen_dist;
    if(*N_e_new<MIN_NE) *N_e_new=MIN_NE;
  }
  if (Par->unlinked_ind==1) *N_e_new = 0.0;
  
  /* CALCULATE SOMETHING ANALAGOUS TO EXPECTED NUMBER OF TIMES EACH HAP i IS VISITED, CONDITIONAL ON THE OBSERVED DATA (I.E  (27) AND PARAGRAPH UNDER (38) IN RABINER 1989, Proceedings of the IEEE 77(2):257-286), BUT -- AS WE'RE ONLY COUNTING CHUNKS -- SUBTRACT OUT TIMES YOU DO NOT SWITCH */
  for (i=0; i < *p_Nhaps; i++) corrected_chunk_count[i]=corrected_chunk_count[i]+copy_prob_newSTART[i];

   free(exp_copy_pop);
   free(BetavecPREV);
   free(BetavecCURRENT);

   free(ind_snp_sum_vec);
   free(expected_transition_prob);
   free(regional_chunk_count);
   free(regional_chunk_count_sum);

}

///////////////////////////////////////////////
double ** sampler(double ** copy_prob_new_mat, int * newh, int ** existing_h, int *p_Nloci, int *p_Nhaps, int *p_nchr, double p_rhobar, double * MutProb_vec, int * allelic_type_count_vec, double * lambda, double * pos, double * copy_prob, double * copy_probSTART, int * pop_vec, int * cond_mat_haplotypes,int ndonorpops, int run_num, int ind_val, struct files_t *Outfiles, struct param_t *Par)
{


  int i, j, locus;
  double sum;
  double prob, total_prob;
  double total_prob_from_i_to_i,total_prob_to_i_exclude_i,total_prob_from_i_exclude_i,total_prob_from_any_to_any_exclude_i;
  int num_regions;
  double * TransProb = malloc( ((*p_Nloci)-1) * sizeof(double));
  double N_e_new;
  int * sample_state = malloc(*p_Nloci * sizeof(int));
             //correction to PAC-A rho_est
  double random_unif, random_unifSWITCH;
  double no_switch_prob;
  double large_num;
  double ** Alphamat = malloc(*p_Nhaps * sizeof(double *));
  double * copy_prob_new = malloc(*p_Nhaps * sizeof(double));
  double * copy_prob_newSTART = malloc(*p_Nhaps * sizeof(double));
  double * Alphasumvec = malloc(*p_Nloci * sizeof(double));
  double * corrected_chunk_count = malloc(*p_Nhaps * sizeof(double));
  double * expected_chunk_length = malloc(*p_Nhaps * sizeof(double));
  double * expected_differences = malloc(*p_Nhaps * sizeof(double));
  double * regional_chunk_count_sum_final = malloc(ndonorpops * sizeof(double));
  double * regional_chunk_count_sum_squared_final = malloc(ndonorpops * sizeof(double));
  double * snp_info_measure=malloc(ndonorpops * sizeof(double));

  if(Par->vverbose) fprintf(Par->out,"        sampler: initializing\n");

  for(i=0 ; i< *p_Nhaps ; i++)
    {
      Alphamat[i] = malloc(*p_Nloci * sizeof(double));
    }
  for (i=0; i < ndonorpops; i++)
    {
      regional_chunk_count_sum_final[i] = 0.0;
      regional_chunk_count_sum_squared_final[i] = 0.0;
      snp_info_measure[i]=0.0;
    }

  if(Par->vverbose) fprintf(Par->out,"        sampler: constructing theta\n");
				// Theta as given in Li and Stephens
  double Theta=thetaLiStephens(MutProb_vec,p_nchr,p_Nhaps);
  if(Par->vverbose) fprintf(Par->out,"(Using theta=%f)\n",Theta);
    // TransProb[i] is probability of copying mechanism "jumping" between
  //   loci i and i+1
  if(Par->vverbose) fprintf(Par->out,"        sampler: constructing delta\n");
  double delta=deltaLiStephens(TransProb,pos,p_rhobar,lambda,p_Nloci,Par);
  if(Par->vverbose) fprintf(Par->out,"(Using delta=%f)\n",delta);

  if(Par->vverbose) fprintf(Par->out,"        sampler: forwards algorithm\n");
      /* FORWARDS ALGORITHM: (Rabiner 1989, p.262) */
  double Alphasum = forwardAlgorithm(newh, existing_h, Alphamat, MutProb_vec, p_Nhaps,p_Nloci,copy_prob, copy_probSTART, TransProb,Par);

  if(Outfiles->usingFile[2]) fprintf(Outfiles->fout3," %.10lf",Alphasum);

  if(Par->vverbose) fprintf(Par->out,"        sampler: backwards algorithm\n");
  int finalrun= (run_num == (Par->EMruns-1));
  if(run_num <= (Par->EMruns-1)){
    backwardAlgorithm(finalrun,ndonorpops,ind_val,Alphasum,p_rhobar,&N_e_new,newh,existing_h,Alphamat,lambda,delta,MutProb_vec,p_Nhaps,p_Nloci,copy_prob,copy_prob_new,copy_prob_newSTART, corrected_chunk_count, expected_chunk_length, expected_differences,regional_chunk_count_sum_final,regional_chunk_count_sum_squared_final, &num_regions, copy_probSTART, TransProb, pop_vec,pos,snp_info_measure,Outfiles,Par);
  }

  //////////////////////////////// 
      /* print-out samples if we've done enough iterations: */
   if (finalrun)
     {
       if(Par->vverbose) fprintf(Par->out,"        sampler: printing.\n");

       N_e_new = p_rhobar;

       for (i=0; i < *p_Nhaps; i++)
	 {
	   copy_prob_new[i] = copy_prob[i];
	   copy_prob_newSTART[i] = copy_probSTART[i];
	 }

           /* calculate Alphasums (for efficient sampling): */
       for (locus=0; locus < *p_Nloci; locus++)
	 {
	   Alphasumvec[locus] = 0.0;
	   large_num = Alphamat[0][locus];
	   for (i = 1; i < *p_Nhaps; i++)
	     {
	       if (Alphamat[i][locus] > large_num)
		 large_num = Alphamat[i][locus];
	     }
	   large_num = -1.0*large_num;
	   for (i = 0; i < *p_Nhaps; i++)
	     Alphasumvec[locus] = Alphasumvec[locus] + exp(Alphamat[i][locus]+large_num);
	   Alphasumvec[locus] = log(Alphasumvec[locus]) - large_num;
	 }

              /* SAMPLING ALGORITHM: (from Falush, Stephens, & Pritchard (2003) Genetics 164:1567-1587) */
       for (j = 0; j < Par->samplesTOT; j++)
	 {
	   //fprintf(Par->out,"sample %d\n",j);
	      /* sample last position: */
	   total_prob = 0.0;
	   large_num = Alphamat[0][(*p_Nloci-1)];
	   for (i = 1; i < *p_Nhaps; i++)
	     {
	       if (Alphamat[i][(*p_Nloci-1)] > large_num)
		 large_num = Alphamat[i][(*p_Nloci-1)];
	     }
	   large_num = -1.0*large_num;
	   random_unif = (double) rand()/RAND_MAX;
	   total_prob = Alphasumvec[(*p_Nloci-1)];
	   prob = 0.0;
	   for (i = 0; i < *p_Nhaps; i++)
	     {
	       prob = prob + exp(Alphamat[i][(*p_Nloci-1)]+large_num);
	       if (random_unif <= exp(log(prob)-large_num-total_prob))
		 {
		   sample_state[(*p_Nloci-1)] = i;
		   break;
		 }
	     }

              /* sample remaining positions: */
	   for (locus = (*p_Nloci-2); locus >= 0; locus--)
	     {
	        // first sample prob you switch and see if you need to
                   // if you do need to switch, you need to go through the below loop to figure out where to switch to
               large_num = -1.0 * Alphasumvec[locus];
	       total_prob = log(exp(Alphasumvec[locus]+large_num)*TransProb[locus]*copy_prob[sample_state[(locus+1)]] + exp(Alphamat[sample_state[(locus+1)]][locus]+large_num)*(1.0-TransProb[locus]))-large_num;
	       no_switch_prob = exp(log(exp(Alphamat[sample_state[(locus+1)]][locus]+large_num)*(1.0-TransProb[locus])) - large_num - total_prob);
               random_unifSWITCH = (double) rand()/RAND_MAX;
	       if (random_unifSWITCH <= no_switch_prob) sample_state[locus] = sample_state[(locus+1)];

	       //if (j ==0 && locus > 9500) fprintf(Par->out,"%d %d %lf %lf %lf %lf %lf\n",locus,sample_state[(locus+1)],no_switch_prob,Alphamat[sample_state[(locus+1)]][locus],large_num,total_prob,1.0-TransProb[locus]);
               if (random_unifSWITCH > no_switch_prob)
		 {
		   total_prob = 0.0;
		   large_num = Alphamat[0][locus];
		   for (i = 1; i < *p_Nhaps; i++)
		     {
		       if (Alphamat[i][locus] > large_num)
			 large_num = Alphamat[i][locus];
		     }
		   large_num = -1.0*large_num;

		   random_unif = (double) rand()/RAND_MAX;
		   total_prob = log(exp(Alphasumvec[locus]+large_num)*TransProb[locus]*copy_prob[sample_state[(locus+1)]]) - large_num;
		   prob = 0.0;
		   for (i = 0; i < *p_Nhaps; i++)
		     {
		       prob = prob + exp(Alphamat[i][locus]+large_num)*TransProb[locus]*copy_prob[sample_state[(locus+1)]];
		       if (random_unif <= exp(log(prob)-large_num-total_prob))
			 {
			   sample_state[locus] = i;
			   break;
			 }
		     }
		 }
	     }

	   if(Outfiles->usingFile[0]) {
	     gzprintf(*Outfiles->fout1,"%d",j+1);
	     for (i = 0; i < *p_Nloci; i++)
	       {
		 gzprintf(*Outfiles->fout1," %d",cond_mat_haplotypes[sample_state[i]]+1);
	       }
	     gzprintf(*Outfiles->fout1,"\n");
	   }
	 }
     }
   if(Par->vverbose) fprintf(Par->out,"        sampler: organising memory.\n");

  for (i = 0; i < *p_Nhaps; i++)
     copy_prob_new_mat[0][i] = copy_prob_new[i];
   for (i = 0; i < *p_Nhaps; i++)
     copy_prob_new_mat[1][i] = copy_prob_newSTART[i];
   for (i = 0; i < *p_Nhaps; i++)
     copy_prob_new_mat[2][i] = corrected_chunk_count[i];
   for (i = 0; i < *p_Nhaps; i++)
     copy_prob_new_mat[3][i] = expected_chunk_length[i];
   for (i = 0; i < *p_Nhaps; i++)
     copy_prob_new_mat[4][i] = expected_differences[i];
   for (i = 0; i < ndonorpops; i++)
     copy_prob_new_mat[5][i] = regional_chunk_count_sum_final[i];
   for (i = 0; i < ndonorpops; i++)
     copy_prob_new_mat[6][i] = regional_chunk_count_sum_squared_final[i];
   for (i = 0; i < ndonorpops; i++)
     copy_prob_new_mat[7][i] = snp_info_measure[i];
   copy_prob_new_mat[0][(*p_Nhaps)] = N_e_new;
   copy_prob_new_mat[1][(*p_Nhaps)] = num_regions;

   if(Par->vverbose) fprintf(Par->out,"        sampler: freeing memory.\n");
   for (i=0; i < *p_Nhaps; i++)
     {
       free(Alphamat[i]);
     }
   free(Alphamat);
   free(TransProb);
   free(sample_state);
   free(Alphasumvec);
   free(copy_prob_new);
   free(copy_prob_newSTART);
   free(corrected_chunk_count);
   free(expected_chunk_length);
   free(expected_differences);
   free(regional_chunk_count_sum_final);
   free(regional_chunk_count_sum_squared_final);
   free(snp_info_measure);

   return(copy_prob_new_mat);
}

///////////////////////////////////////
///////////////////////////////////////
// Log likelihood

int loglik(struct copyvec_t *Copyvec, struct donor_t *Donors, struct data_t *Data,struct ids_t *Ids, struct infiles_t *Infiles, struct files_t *Outfiles, struct param_t * Par)
  //int loglik(int nhaps_startpop, int *p_nloci, int p_nhaps, double N_e_start, double * recom_map, double * MutProb_vec, double * copy_prob, double * copy_probSTART, int * pop_vec, char *filename, struct files_t *Outfiles, struct param_t * Par)
{
  int nhaps_condpop = Data->condhaps;
  int nind_condpop  = Data->condhaps /Data->hapsperind;

  char *step;
  char line[2047];
  char waste[2047];
  int i, j, m, n, count, r, h;
  int nhaps, num_regions_tot;
  int ndonorpops=Outfiles->Donors->ndonorpops;
  double sum_total_diff;

  double * total_back_prob = malloc((ndonorpops) * sizeof(double));
  double * total_back_probSTART = malloc((ndonorpops) * sizeof(double));
  double * total_counts = malloc((ndonorpops) * sizeof(double));
  double * total_lengths = malloc((ndonorpops) * sizeof(double));
  double * total_differences = malloc((ndonorpops) * sizeof(double));
  double * total_region_counts = malloc((ndonorpops) * sizeof(double));
  double * total_squared_region_counts = malloc((ndonorpops) * sizeof(double));
  double * snp_info_measure_final = malloc((ndonorpops) * sizeof(double));
  double N_e_new, N_e;
  double total_prob, total_probSTART;
  int * copy_pop_vec = malloc(Copyvec->ndonors * sizeof(int));
  double * copy_prob_new = malloc(Copyvec->ndonors * sizeof(double));
  double * copy_prob_newSTART = malloc(Copyvec->ndonors * sizeof(double));
  double * MutProb_vec_new = malloc(Copyvec->ndonors * sizeof(double));
  double ** back_prob;// Allocated in sampler // = malloc(8 * sizeof(double *));
  double ** copy_prob_pop = malloc(2 * sizeof(double *));
  int * ndonorhaps_vec=malloc((ndonorpops)*sizeof(int));

  //  for (i=0; i < 8; i++)
  //    back_prob[i] = malloc((Data->condhaps+1) * sizeof(double));
  for (i=0; i < 2; i++)
    copy_prob_pop[i] = malloc((ndonorpops) * sizeof(double));

  if(Par->vverbose) fprintf(Par->out,"sampler::loglik allocated memory, using %d donor pops\n",ndonorpops);

  for (i=0; i < (ndonorpops); i++)
    {
      total_back_prob[i] = 0.0;
      total_back_probSTART[i] = 0.0;
    }

  if(Par->verbose) fprintf(Par->out,"Constructed basic likelihood quantities; processing from ind %d to %d\n",Par->start_val+1,Par->end_val);

  total_prob = 0.0;
  total_probSTART = 0.0;
  for (m = Par->start_val; m < Par->end_val; m ++)
    {
      setIndAsRecipient(m,Data,Ids,Par);// set up data to paint the m-th individual
      nhaps = Data->current_donor_nind*Data->hapsperind;
      int nhapseff=Donors->ndonorhaps + Data->reciphaps;

      int indon=0;
      for (i=0; i < ndonorpops; i++) {
	ndonorhaps_vec[i]=Outfiles->Donors->ndonorhaps_vec[i];
	
	if(!Par->geno2_find &&
	   Donors->recippops[m]<Donors->nrecippops &&
	   Donors->recippops[m]>=0 &&
	   strcmp(Donors->donorlabels[i],Donors->recippoplabels[Donors->recippops[m]])==0) {
	  if(Par->verbose) fprintf(Par->out,"Painting individual from population %s, removing them from it as a donor\n",Donors->donorlabels[i]);
	  ndonorhaps_vec[i]-=Data->hapsperind;
	  //nhaps-=Data->hapsperind;
	  nhapseff=Data->condhaps;
	}
	for (j=0; j <ndonorhaps_vec[i]; j++){
	  copy_pop_vec[indon++]=i;
	}
      }

      int * allelic_type_count_vec=getallelic_type_count_vec(Data);

      fprintf(Par->out,"Processing Recipient number %d, named %s\n",m+1,Donors->reciplabels[m]);
      if(Outfiles->usingFile[2]) fprintf(Outfiles->fout3,"EMPAR for %s\n",Donors->reciplabels[m]);

      indon=0;
      for (i=0; i < ndonorpops; i++)
	{
	  for (j=0; j < ndonorhaps_vec[i]; j++)
	    {
	      if (Par->prior_donor_probs_ind==0){
		copy_prob_new[indon] = 1.0/nhaps;
	      }else{
		copy_prob_new[indon] = Donors->ndonorprobs[i]/ndonorhaps_vec[i];
	      }
		     indon++;
	    }
	}
      
      for (i=0; i < Copyvec->ndonors; i++)
	{
	  copy_prob_newSTART[i] = copy_prob_new[i];
	  MutProb_vec_new[i] = Copyvec->MutProb_vec[i];
	  //	  if(Par->vverbose) fprintf(Par->out,"Mutation rate %f for pop %i\n",MutProb_vec_new[i],i);
	  if(MutProb_vec_new[i]<MIN_THETA && MutProb_vec_new[i]>=0) {
	    MutProb_vec_new[i]=MIN_THETA;
	    if(Par->vverbose) fprintf(Par->out,"Mutation rate %f for ind %i below threshold, set to %f\n",MutProb_vec_new[i],i,MIN_THETA);
	  }
	}

      int finalrun=0;

      for (r=0; r < Par->EMruns; r++)
	{
	  if(r==Par->EMruns-1) finalrun=1;
	  if(Par->verbose) fprintf(Par->out,"EM iteration %i (of %i) for individual %i\n",r+1,Par->EMruns,m+1);
	  if(Outfiles->usingFile[2]) fprintf(Outfiles->fout3,"%d",r);

	  total_prob = 0.0;
	  total_probSTART = 0.0;
	  //	  for (i=0; i < (ndonorpops); i++)
	  for (i=0; i < ndonorpops; i++)
	    {
	      total_back_prob[i] = 0.0;
	      total_back_probSTART[i] = 0.0;
	      total_counts[i]=0.0;
	      total_lengths[i]=0.0;
	      total_differences[i]=0.0;
	      total_region_counts[i]=0.0;
	      total_squared_region_counts[i]=0.0;
	      snp_info_measure_final[i]=0.0;
	    }

	  N_e_new=0.0;
	  num_regions_tot=0;
	  for(h=0; h < Data->hapsperind; h++)
	    {
	      if(Par->vverbose) {
		fprintf(Par->out,"     ... processing haplotype %i of %i : ",h+1,(2-Par->haploid_ind));
		int t;
	  	for(t=0;t<10;++t)fprintf(Par->out,"%i",
					 Data->ind_chromosomes[h][t]);
		fprintf(Par->out," ...\n");
	      }
	      if (r==0) { // recompute Ne, in case this recipient has different number of donors than others
		if (Par->ne_find==0){ N_e=Par->readN_e/nhaps;
		}else N_e=Par->N_e;
	      }


	      // write the headers
	      if (r == (Par->EMruns-1)){
		if(Outfiles->usingFile[0] && Par->samplesTOT > 0) {
		  gzprintf(*Outfiles->fout1,"HAP %d %s\n",h+1, Donors->reciplabels[m]);
		}
		if (Outfiles->usingFile[8] && Par->print_file9_ind==1) gzprintf(*Outfiles->fout9,"HAP %d %s\n",h+1, Donors->reciplabels[m]);
	      }
	      /* SAMPLE FROM PAC CONDITIONAL ON COPY-PROBS: */
	      if(Par->verbose) fprintf(Par->out,"     ... performing sample\n");
	      if(Par->vverbose) {
		fprintf(Par->out,"   using %d SNPs, %d donor haps, %d total haps, Ne=%f, Mutation prob=%f\n",
			Data->nsnps,nhaps,
			Data->nhapstotal,N_e,MutProb_vec_new[0]);
	      }


	      double ** back_prob = malloc(8 * sizeof(double *));
	      for (i=0; i < 8; i++)
		back_prob[i] = malloc((nhaps+1) * sizeof(double));

	      back_prob = sampler(back_prob,
				  Data->ind_chromosomes[h], 
				  Data->cond_chromosomes, 
				  &Data->nsnps,
				  &nhaps, 
				  &nhapseff,
				  N_e, MutProb_vec_new, allelic_type_count_vec, Copyvec->recom_map, Data->positions, copy_prob_new, copy_prob_newSTART, copy_pop_vec, Data->current_donor_haps, Donors->ndonorpops, r, m, Outfiles,Par);

    
	      N_e_new=N_e_new+back_prob[0][nhaps];
	      num_regions_tot=num_regions_tot+back_prob[1][nhaps];

	      if(Par->vverbose) fprintf(Par->out,"     ... processing sample\n");

                 /* GET NEW COPY-PROBS BASED ON PAC SAMPLES: */
	      count = 0;
	      for (i = 0; i < (Donors->ndonorpops); i++)
		{
		  for (j = 0; j < ndonorhaps_vec[i]; j++)
		    {
		      total_back_prob[i] = total_back_prob[i] + back_prob[0][count];
		      total_prob = total_prob + back_prob[0][count];

		      total_back_probSTART[i] = total_back_probSTART[i] + back_prob[1][count];
		      total_probSTART = total_probSTART + back_prob[1][count];

		      total_counts[i] = total_counts[i] + back_prob[2][count];
		      total_lengths[i] = total_lengths[i] + back_prob[3][count];
		      total_differences[i] = total_differences[i] + back_prob[4][count];
		      //if (i < 2) fprintf(Par->out,"%d %d %d %d %d %d %lf %lf %lf %lf\n",m,r,h,i,j,count,back_prob[2][count],back_prob[3][count],total_counts[i],total_lengths[i]);

		      count = count + 1;
		    }
		}
	      for (i = 0; i < (Donors->ndonorpops); i++)
		{
		  total_region_counts[i] = total_region_counts[i] + back_prob[5][i];
		  total_squared_region_counts[i] = total_squared_region_counts[i] + back_prob[6][i];
		  snp_info_measure_final[i] = snp_info_measure_final[i] + back_prob[7][i];
		}

	      for (i = 0; i < 8; i++)
		free(back_prob[i]);
	      free(back_prob);
	    }
	  if(Outfiles->usingFile[2]) fprintf(Outfiles->fout3," %.10lf %.10lf\n",N_e,MutProb_vec_new[0]);
	  //if (estimate_mutationALL_ind==0) fprintf(fout3,"\n");
	  //if (estimate_mutationALL_ind==1) fprintf(fout3," %lf\n",MutProb_vec_new[0]);
	  if (Par->recom_em_find==1) N_e=N_e_new/Data->hapsperind;

	  for (i = 0; i < (Donors->ndonorpops); i++)
	    {
	      copy_prob_pop[0][i] = total_back_prob[i]/total_prob;
	      copy_prob_pop[1][i] = total_back_probSTART[i]/total_probSTART;
	    }

	  if(Par->vverbose) fprintf(Par->out,"creating copy probs and mutation probs\n");
	  /* RESET COPY-PROBS and MUTATION-PROBS: */
	  // (first check for probabilities of 0:)
	  for (i=0; i < (Donors->ndonorpops); i++)
	    {
	      if (copy_prob_pop[0][i] <= 0)
		copy_prob_pop[0][i] = Par->small_copy_val*ndonorhaps_vec[i];

	      if (copy_prob_pop[1][i] <= 0)
		copy_prob_pop[1][i] = Par->small_copy_val*ndonorhaps_vec[i];
	    }
	  total_prob = 0.0;
	  total_probSTART = 0.0;
	  for (j=0; j < (Donors->ndonorpops); j++)
	    {
	      total_prob = total_prob + copy_prob_pop[0][j];
	      total_probSTART = total_probSTART + copy_prob_pop[1][j];
	    }
	  for (j=0; j < (Donors->ndonorpops); j++)
	    {
	      copy_prob_pop[0][j] = copy_prob_pop[0][j]/total_prob;
	      copy_prob_pop[1][j] = copy_prob_pop[1][j]/total_probSTART;
	    }

	  if (Par->copy_prop_em_find==1)
	    {
	      if(Par->vverbose) fprintf(Par->out,"creating copy probs\n");
	      count = 0;
	      for (i=0; i < (ndonorpops); i++)
		{
		  for (j=0; j < ndonorhaps_vec[i]; j++)
		    {
		      copy_prob_new[count] = copy_prob_pop[0][i]/ndonorhaps_vec[i];
		      copy_prob_newSTART[count] = copy_prob_pop[1][i]/ndonorhaps_vec[i];
		      count = count + 1;
		    }
		}
	    }

	  if (Par->mutation_em_find==1)
	    {
	      if(Par->vverbose) fprintf(Par->out,"creating mutation probs\n");
	      count = 0;
	      for (i=0; i < (ndonorpops); i++)
		{
		  for (j=0; j < ndonorhaps_vec[i]; j++)
		    {
		      MutProb_vec_new[count] = total_differences[i]/(Data->nsnps*(2-Par->haploid_ind));
		      count = count + 1;
		    }
		}
	    }

	  if (Par->mutationALL_em_find==1)
	    {
	      sum_total_diff=0.0;
	      for (i=0; i < (Donors->ndonorpops); i++)
		sum_total_diff=sum_total_diff+total_differences[i]/(Data->nsnps*(2-Par->haploid_ind));
	      count = 0;
	      for (i=0; i < (Donors->ndonorpops); i++)
		{
		  for (j=0; j < ndonorhaps_vec[i]; j++)
		    {
		      MutProb_vec_new[count] = sum_total_diff;
		      count = count + 1;
		    }
		}
	    }

	  /* print props, lengths, counts, and differences: */
 	  if (finalrun) printSummary(m,num_regions_tot,copy_prob_pop,total_counts,total_lengths,total_differences,total_region_counts,total_squared_region_counts,Outfiles,Par);
	}
      free(allelic_type_count_vec);
    }

  for(i=0;i<2;i++) free(copy_prob_pop[i]);
  free(copy_prob_pop);
  free(copy_pop_vec);
  free(copy_prob_new);
  free(copy_prob_newSTART);
  free(MutProb_vec_new);
  //  for (i = 0; i < 8; i++)
  //    free(back_prob[i]);
  //  free(back_prob);
  free(total_back_prob);
  free(total_back_probSTART);
  free(ndonorhaps_vec);
  free(total_counts);
  free(total_lengths);
  free(total_differences);
  free(total_region_counts);
  free(total_squared_region_counts);
  free(snp_info_measure_final);

  return(1);
}
