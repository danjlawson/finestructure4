#include "ChromoPainterRecmap.h"

//int ** p_pop_vec,double ** p_MutProb_vec,double ** p_copy_prob,double ** p_copy_probSTART, 
// int ndonorpops, int *ndonorhaps,double *ndonorprobs, double *ndonormutrates

struct copyvec_t * initializeCopyVecs(struct donor_t *Donors,struct data_t *Data, struct param_t * Par) {

  struct copyvec_t *Copyvec=malloc(sizeof(struct copyvec_t));

  int ndonors=0;
  int i,j,count;
  // compute the number of donor haplotypes
  for(i=0;i<Donors->ndonorpops;i++){
    ndonors+=Donors->ndonorhaps_vec[i];
  }
  if (ndonors != Data->condhaps)
    {
      fprintf(Par->out,"Found %d donor haplotypes, but %d haplotypes in the phase file. These must match if used. Exiting....\n",ndonors,Data->condhaps);
      //      stop_on_error(1,Par->errormode,Par->err);
    }

  // allocate memory
  Copyvec->ndonors=ndonors;
  Copyvec->pop_vec= malloc(ndonors * sizeof(int));
  Copyvec->MutProb_vec = malloc(ndonors * sizeof(double));
  Copyvec->copy_prob = malloc(ndonors * sizeof(double));
  Copyvec->copy_probSTART = malloc(ndonors * sizeof(double));

  // (0) INITIALIZE copy_prob, MutProb_vec, and pop_vec:
  if (Par->prior_donor_probs_ind==0) { // if we don't have prior info for the probs, its simple
    for (i=0; i < Copyvec->ndonors; i++) Copyvec->copy_prob[i] = 1.0/Copyvec->ndonors;
    if(Par->vverbose) fprintf(Par->out,"Assigning copy probability %f to each of the %d donors\n",1.0/Copyvec->ndonors,Copyvec->ndonors);
  }

  if (Par->prior_donor_probs_ind==1)
  {// account for the prior probabilities, assigning uniformly within populations
    if(Par->vverbose) fprintf(Par->out,"Assigning assigned donor probabilities uniformally within each population\n");
    count=0;
    for (i=0; i < Donors->ndonorpops; i++)
      {
	for (j=0; j < Donors->ndonorhaps_vec[i]; j++)
	  {
	    Copyvec->copy_prob[count] = Donors->ndonorprobs[i]/Donors->ndonorhaps_vec[i];
	    count = count + 1;
	  }
      }
  }
  
  ///////////////////////////////
  // mutation probs
  if ((Par->mutation_rate_ind==0) || (Par->mutationALL_em_find==1))
    if(Par->vverbose) fprintf(Par->out,"Assigning global mutation rate %f to the %d donors\n",Par->GlobalMutRate,ndonors);
   {// global mutation rate is being used
     for (i=0; i < ndonors; i++) {
       Copyvec->MutProb_vec[i]=Par->GlobalMutRate;
     }
   }
 
  if ((Par->mutation_rate_ind==1) && (Par->mutationALL_em_find==0))
    { // local mutation rate
    if(Par->vverbose) fprintf(Par->out,"Assigning local mutation rates\n");
    count=0;
    for (i=0; i < Donors->ndonorpops; i++)
      {
	for (j=0; j < Donors->ndonorhaps_vec[i]; j++)
	  {
	    Copyvec->MutProb_vec[count] = Donors->ndonormutrates[i];
	    count = count + 1;
	  }
      }
    }

  ///////////////////////////////
  // Copy probability
  if(Par->vverbose) fprintf(Par->out,"Creating initial copying probability vectors\n");
  for (i=0; i < ndonors; i++)
    Copyvec->copy_probSTART[i] = Copyvec->copy_prob[i];

  count = 0;
  for (i=0; i < Donors->ndonorpops; i++)
    {
      for (j=0; j < Donors->ndonorhaps_vec[i]; j++)
	{
	  Copyvec->pop_vec[count] = i;
	  count = count + 1;
	}
    } 
  return Copyvec;
}


void assignRecMap(struct copyvec_t *Copyvec, struct infiles_t *Infiles,struct data_t *Data, struct param_t *Par) {
  ///////////////////////////////////
  // READING IN THE RECOMBINATION MAP
  char *step;
  char line[2047];
  double bpval;
  int i=0, j=0;

  if(Par->vverbose) fprintf(Par->out, "Assigning recombination map for %d SNPs\n",Data->nsnps);
  Copyvec->recom_map_size = Data->nsnps-1;
  Copyvec->recom_map = malloc(Copyvec->recom_map_size * sizeof(double));
  if ((Par->unlinked_ind==1) && (Par->recom_find==0))
    {
      for (j=0; j < (Copyvec->recom_map_size); j++) Copyvec->recom_map[j]=-9.0;
    }
  if (Par->recom_find==1)
    {
      if (!openRecmap(Infiles,Par->out)) { fprintf(Par->out,"error opening recom map input file: %s\n",Infiles->recmap); stop_on_error(1,Par->errormode,Par->err);}
      strcpy(line,"");
 if(fgets(line,2047,Infiles->frecmap)==NULL){stop_on_error(1,Par->errormode,Par->err);};   // header
      //
      int jj=-1;
      j=-1;
      while(jj < Copyvec->recom_map_size-1)
	{
	  if(fgets(line,2047,Infiles->frecmap)==NULL){stop_on_error(1,Par->errormode,Par->err);};
	  j++;
	  if(j<Par->startlocus || j>=Par->endlocus) continue;
	  jj++;
	  step=line;
	  reading(&step,"%lf",&bpval);    // basepair position

	  if (bpval != Data->positions[jj])
	    {
	      if(Par->jitter_locations)
		fprintf(Par->out,"Warning: genetic map position difference at basepair %d (%lf vs %lf). This will occur if jittering occurred or if using an invalid map.\n",jj+1,bpval,Data->positions[jj]);
	      else {
		fprintf(Par->out,"basepair positions do not match between %s and phase file at basepair %d (%lf vs %lf). Exiting....\n",Infiles->recmap,jj+1,bpval,Data->positions[jj]);
		stop_on_error(1,Par->errormode,Par->err);
	      }
	    }
	  reading(&step,"%lf",&Copyvec->recom_map[jj]);
	  if (Copyvec->recom_map[jj] >= 0 && Copyvec->recom_map[jj] <= Par->small_recom_val)
	    {
	      fprintf(Par->out,"Warning: recom rate very low at basepair %lf (%lf). Assuming recomb rate between this snp and next one is %lf....\n",Data->positions[jj],Copyvec->recom_map[jj],Par->small_recom_val);
	      Copyvec->recom_map[jj]=Par->small_recom_val;
	    }
	  if (Copyvec->recom_map[jj]<0)
	    {
	      fprintf(Par->out,"recom rate < 0 at basepair %lf. Assuming recomb rate of infinity between this snp and next one....\n",Data->positions[jj]);
	      //fprintf(Par->out,"recom rate must be > 0 (basepair %lf)!! Exiting....\n",Data->positions[j]);
	      //stop_on_error(1,Par->errormode,Par->err);
	    }
	}
      if(fgets(line,2047,Infiles->frecmap)==NULL){stop_on_error(1,Par->errormode,Par->err);};
      closeRecmap(Infiles);
    }

  // check ordering of snps (only allowed to be less than previous position if recom_map<0 at position -- i.e. suggesting new chromosome):
  for (i=1; i < Copyvec->recom_map_size; i++)
    {
      if (Data->positions[i]<Data->positions[(i-1)] && (Copyvec->recom_map[(i-1)]>=0))
	{
	  if(Par->unlinked_ind==1){
	    fprintf(Par->out,"WARNING: positions in phase file are not listed in increasing order between SNPs %i-%i (basepairs %lf and %lf). Continuing as loci are unlinked....\n",i-1,i,Data->positions[(i-1)],Data->positions[i]);
	  }else{
	    fprintf(Par->out,"positions in phase file are not listed in increasing order at basepairs %lf and %lf. Exiting....\n",Data->positions[(i-1)],Data->positions[i]);
	    stop_on_error(1,Par->errormode,Par->err);
	  }
	}
    }
}


void clearCopyvec(struct copyvec_t *Copyvec) 
{
  if(Copyvec==NULL) return;
  if(Copyvec->recom_map_size>=0) free(Copyvec->recom_map);
  if(Copyvec->ndonors>=0) {
    free(Copyvec->MutProb_vec);
    free(Copyvec->pop_vec);
    free(Copyvec->copy_prob);
    free(Copyvec->copy_probSTART);
  }
  free(Copyvec);
}
