#include "ChromoPainterDonors.h"

void reallocateDataHaps(struct data_t *Data,struct donor_t *Donors, struct ids_t *Ids, struct param_t *Par) {
  free(Data->cond_chromosomes);
  free(Data->condnums);
  free(Data->recipnums);
  // quantities that we now know after reading the donor file
  Data->condhaps  = Donors->ndonorhaps;
  Data->reciphaps = Donors->nrecipinds*Data->hapsperind;
  Data->recipinds = Donors->nrecipinds;
  if(Par->vverbose) fprintf(Par->out,"reallocateDataHaps: Using %d cond haps, %d recipient individuals\n",Donors->ndonorhaps,Donors->nrecipinds);
  
  //allocate memory
  Data->cond_chromosomes=malloc(Data->condhaps*sizeof(int *));
  Data->condnums  = malloc((Data->condhaps)/Data->hapsperind*sizeof(int));
  Data->recipnums = malloc((Data->reciphaps)/Data->hapsperind*sizeof(int));


  
  // assign conditioning haplotypes
  int i,k,ii=0,iii=0;
  for(i=0;i<Donors->ndonorpops;++i){ // i is the donor pop, k is the raw ID index. ii is the index within the donors and iii the haplotype index...
    for(k=0;k<Ids->nind_tot;++k){
      while(k<Ids->nind_tot && Ids->include_ind_vec[k]==0) ++k;
      if(k>=Ids->nind_tot) break;
      char *refid=Ids->popid[k];
      if(Par->allvsall && Par->geno2_find) refid=Ids->id[k];
      if(strcmp(refid,Donors->donorlabels[i])==0){
	Data->condnums[ii]=k;
	if(Par->vverbose)fprintf(Par->out,"Assigning donor chromosome %i to raw data chromosome %i\n",iii,k*Data->hapsperind);
	Data->cond_chromosomes[iii]=Data->all_chromosomes[k*Data->hapsperind];
	if(Data->hapsperind==2) {
	  Data->cond_chromosomes[iii+1]= Data->all_chromosomes[k*Data->hapsperind+1];
	}
	++ii;
	iii+=Data->hapsperind;
	if(Par->vverbose){
	  fprintf(Par->out,"Donors::reallocateDataHaps Assigned data structures for donor individual %i (called %s) in population %s : ",k,Ids->id[k],Ids->popid[k]);
	  int t;
	  for(t=0;t<10;++t)fprintf(Par->out,"%i",Data->all_chromosomes[k*Data->hapsperind][t]);
	  fprintf(Par->out,"\n");
	}
      }
    }
  }
  if(iii!=Data->condhaps){
    fprintf(Par->out,"Donors::reallocateDataHaps Logic error: Expected to count %i donor inds but counted %i instead\n",Data->condhaps/Data->hapsperind,ii);
    stop_on_error(1,Par->errormode,Par->err);
  }

  // assign recipient haplotypes
  ii=0;
  for(i=0;i<Donors->nrecipinds;++i){
      for(k=0;k<Ids->nind_tot;++k){
	while(k<Ids->nind_tot && Ids->include_ind_vec[k]==0) ++k;
	if(k>=Ids->nind_tot) break;
	char *refid=Ids->id[k];
	if(Par->geno2_find) refid=Ids->id[k]; 
	if(strcmp(refid,Donors->reciplabels[i])==0){
	  Data->recipnums[ii]=k;
	  ++ii;
	  if(Par->vverbose){
	    fprintf(Par->out,"Donors::reallocateDataHaps Assigned data structures for recipient individual %i, location %i in the ID file (called %s) :",ii,k,Ids->id[k]);
	    int t;
	    for(t=0;t<10;++t)fprintf(Par->out,"%i",Data->all_chromosomes[k*Data->hapsperind][t]);
	    fprintf(Par->out,"\n");

	  }
	}
      }
  }
  if(ii!=Data->recipinds){
    fprintf(Par->out,"Donors::reallocateDataHaps Logic error: Expected to count %i recipient inds but counted %i instead\n",Data->recipinds,ii);
    stop_on_error(1,Par->errormode,Par->err);
  }
  if(Par->vverbose) fprintf(Par->out,"reallocateDataHaps: Using %d cond haps, %d recipient population(s)\n",Data->condhaps,Donors->nrecippops);
  
  // Populate the data 
  setIndAsRecipient(0,Data,Ids,Par);
  if(Par->vverbose){
     fprintf(Par->out,"Donors::reallocateDataHaps Assigned %i condition haplotypes\n",Data->condhaps);
  }

}

struct donor_t *createDonors(struct infiles_t *Infiles,struct ids_t *Ids,struct ids_t *Ids2,struct data_t *Data,struct param_t *Par) {
  struct donor_t *Donors=malloc(sizeof(struct donor_t)); 
 
  if(Par->donorlist_find){// donor mode
    readDonorPops(Donors,Ids,Ids2,Infiles,Data,Par);
    if(Par->vverbose) fprintf(Par->out,"Found %i donor populations in donor file\n",Donors->ndonorpops);
  }else{// all-vs-all mode
    if(Par->vverbose) fprintf(Par->out,"Creating Donors in all-vs-all mode\n");
    populateDonors(Donors,Ids,Ids2,Data,Par);
  }
  
  if (Par->end_val==0) Par->end_val=Data->recipinds;
  if(Par->end_val-1 > Data->nhapstotal/Data->hapsperind) {
    fprintf(Par->out,"Donors::createDonors error: Requested a final individual (%d, via -a?) greater than the number of individuals (%d) provided!\n",Par->end_val-1, Data->nhapstotal/Data->hapsperind);
    stop_on_error(1,Par->errormode,Par->err);
  }
  return(Donors);
}

int getNpops(struct infiles_t *Infiles,struct param_t * Par, bool countdonors) {
  ////////////////////////////////////////////////
  // open third file (to get information on donor population hap numbers)
  // returns -1 if the file is not specified
  char line[2047];
  char *step;
  char *tpopname=malloc(2047*sizeof(char));
  char *tpopstatus=malloc(2047*sizeof(char));

  int npops=0;
  if (Par->donorlist_find) {
    if(Par->verbose) fprintf(Par->out,"Donors::getNpops Opening Donor list %s\n",Infiles->donorlist);
    if(!openDonorlist(Infiles,Par->out)){ fprintf(Par->out,"Donors::getNpops error opening %s\n",Infiles->donorlist); stop_on_error(1,Par->errormode,Par->err);}
    if(fgets(line,2047,Infiles->fdonorlist)==NULL &&!feof(Infiles->fdonorlist)){printf("Donors::getNpops Error reading %s",Infiles->donorlist);stop_on_error(1,Par->errormode,Par->err);};
    while(!feof(Infiles->fdonorlist) || strlen(line)>0)
      {
	step = line;
	reading(&step,"%s",tpopname);
	reading(&step,"%s",tpopstatus);
	if(tpopstatus[0]=='\n') {printf("Donors::getNpops Error reading %s",Infiles->donorlist);stop_on_error(1,Par->errormode,Par->err);}
	if(tpopstatus[0]=='D' && countdonors){
	  ++npops;
	}else if(tpopstatus[0]!='D' && !countdonors){
	  ++npops;
	}
	if(feof(Infiles->fdonorlist)) break;
	strcpy(line,"");
	if(fgets(line,2047,Infiles->fdonorlist)==NULL &&!feof(Infiles->fdonorlist)){printf("Donors::getNpops Error reading %s",Infiles->donorlist);stop_on_error(1,Par->errormode,Par->err);};
      }
    closeDonorlist(Infiles);
    if(Par->vverbose && countdonors) fprintf(Par->out, "Donors::getNpops Counted %i donor population(s)\n",npops);
    if(Par->vverbose && !countdonors) fprintf(Par->out,"Donors::getNpops Counted %i recipient population(s)\n",npops);
  }
  free(tpopname);
  free(tpopstatus);
  return(npops);
}
  /////////////////////////////////////////////
  // 

//&Donors->reciplabels,&Donors->donorlabels,&Donors->ndonorhaps,&Donors->ndonorprobs,&Donors->ndonormutrates
//char *** p_reciplabels,char *** p_donorlabels, int ** p_ndonorhaps,double ** p_ndonorprobs, double ** p_ndonormutrates
void readDonorPops(struct donor_t *Donors,struct ids_t *Ids,struct ids_t *Ids2,struct infiles_t *Infiles,struct data_t *Data,struct param_t *Par) {
  double totaldonorprobs=0;
  int i=0,j=0;
  char line[2047];
  char *step;

  if (!Par->donorlist_find) {
    fprintf(Par->out,"Donors::readDonorPops Logic error: we should not read donor pops when not using the donor pops file!\n"); stop_on_error(1,Par->errormode,Par->err);
  }
  if(Par->allvsall ==1) {
    // The number of donor pops is determined by all-vs-all mode
    fprintf(Par->out,"Donors::readDonorPops Logic error: we should not read donor pops when in all-vs-all mode!\n"); stop_on_error(1,Par->errormode,Par->err);
  }
  
  /// Get the number of donor populations, and create all variables needed
  Donors->nrecipinds=(int) ceil(Data->nind);
  Donors->ndonorpops=getNpops(Infiles,Par,true);
  Donors->nrecippops=getNpops(Infiles,Par,false);
  int usingdonorfile=1;
  if(Donors->nrecippops==0) {
    Donors->nrecippops=1;
    usingdonorfile=0;
  }
  Donors->ndonorhaps=0;/// Will be calculated later

  // allocate memory
  Donors->ndonorhaps_vec = malloc(Donors->ndonorpops * sizeof(int));
  Donors->donorlabels =malloc(Donors->ndonorpops * sizeof(char *));
  Donors->recippoplabels =malloc(Donors->nrecippops * sizeof(char *));
  for(i=0;i<Donors->ndonorpops;i++)
    Donors->donorlabels[i]=malloc(2047*sizeof(char));
  for(i=0;i<Donors->nrecippops;i++)
    Donors->recippoplabels[i]=malloc(2047*sizeof(char));
  
  if (Par->prior_donor_probs_ind==1) {
    Donors->ndonorprobs = malloc(Donors->ndonorpops * sizeof(double));
  }
  if (Par->mutation_rate_ind==1) {
    Donors->ndonormutrates = malloc(Donors->ndonorpops * sizeof(double));
  }

  if(Par->vverbose) {
    fprintf(Par->out,"Donors::readDonorPops created %i recipient population(s), %i donor population(s)\n",Donors->nrecippops,Donors->ndonorpops);
  }

  // Now process them
  if(!openDonorlist(Infiles,Par->out)) { fprintf(Par->out,"Donors::readDonorPops error opening %s\n",Infiles->donorlist); stop_on_error(1,Par->errormode,Par->err);}
  int ndonors = 0; // number of assigned donor haplotypes
  int nrecips = 0; // number of assigned recipient populations
  //  int nlines=Donors->ndonorpops;
  //for (i=0; i < nlines; i++)
  i=-1;
  while(!feof(Infiles->fdonorlist) || strlen(line)>0)
    {
      strcpy(line,"");
      if(fgets(line,2047,Infiles->fdonorlist)==NULL){continue;}
      if(line[0]=='\n') continue;
      i++;
      if(i>Donors->ndonorpops+Donors->nrecippops*usingdonorfile){
	printf("Donors::readDonorPops error: More populations in file than expected!\n");
	stop_on_error(1,Par->errormode,Par->err);
      }

      step=line;
      if(Par->vverbose) {
	printf("Donors::readDonorPops line %i:%s",i,line);
      }
      char *tpopname=malloc(2047*sizeof(char));
      char *tpopstatus=malloc(2047*sizeof(char));
      
      reading(&step,"%s",tpopname);
      reading(&step,"%s",tpopstatus);
      if(tpopstatus[0]=='D'){
	strcpy(Donors->donorlabels[ndonors],tpopname);
	if (Par->prior_donor_probs_ind==1) // there should be another column of donors
	  {
	    reading(&step,"%lf",&Donors->ndonorprobs[ndonors]);
	    if (Donors->ndonorprobs[ndonors]<=0.0000000000000001)
	      {
		fprintf(Par->out,"Donors::readDonorPops Donor copying probabilities must be > 0. Exiting...\n");
		stop_on_error(1,Par->errormode,Par->err);
	      }
	    totaldonorprobs=totaldonorprobs+Donors->ndonorprobs[ndonors];
	  }
	if (Par->mutation_rate_ind==1) // there should be another column of mutations
	  {
	    reading(&step,"%lf",&Donors->ndonormutrates[ndonors]);
	    if (Donors->ndonormutrates[ndonors]>1)
	      {
		fprintf(Par->out,"Donors::readDonorPops Donor mutation (emission) probabilities must be <= 1 (use a negative number to specify default). Exiting...\n");
		stop_on_error(1,Par->errormode,Par->err);
	      }
	  }
      
	if(Par->vverbose)
	  {
	    fprintf(Par->out,"Donors::readDonorPops Donor Population %i label %s",ndonors, Donors->donorlabels[ndonors]);
	    if(Par->prior_donor_probs_ind) printf(" probability %f",Donors->ndonorprobs[ndonors]);
	    if(Par->mutation_rate_ind) printf(" mutation rate %f",Donors->ndonormutrates[ndonors]);
	    printf("\n");
	}
	++ndonors;
      }else{ // a recipient population
	strcpy(Donors->recippoplabels[nrecips],tpopname);
	if(Par->vverbose)
	{
	  fprintf(Par->out,"Donors::readDonorPops Recipient population %i label %s\n",nrecips, Donors->recippoplabels[nrecips]);
	}
	++nrecips;
      }
      if(feof(Infiles->fdonorlist)) break;
    }
  if(usingdonorfile==0){
	strcpy(Donors->recippoplabels[0],"RecipFile");
  }
  
  // remove unused individuals
  for(j=0;j<Ids->nind_tot;++j){
    int rem=1;
    if(Ids->include_ind_vec[j]==0) rem=0;
    for(i=0;i<Donors->ndonorpops;++i) if(strcmp(Donors->donorlabels[i],Ids->popid[j])==0) rem=0;
    for(i=0;i<Donors->nrecippops;++i) if(strcmp(Donors->recippoplabels[i],Ids->popid[j])==0) rem=0;
    if(usingdonorfile==0 && j>=Ids->nind_tot-Ids2->nind_tot) rem=0; // Don't remove individuals from the target population 
    // If we get here, this individual is supposed to be retained but isn't wanted in the analysis
    if(rem){
      if(Par->verbose) fprintf(Par->out,"Donors::readDonorPops omitting individual %i (called %s) as its population %s is not used as either a donor or a recipient.\n",j,Ids->id[j],Ids->popid[j]);
      Ids->include_ind_vec[j]=0;
      Ids->nind_inc--;
    }
  }

  // assign haplotypes to donors
  for(i=0;i<Donors->ndonorpops;++i) {
    Donors->ndonorhaps_vec[i]=0;
    for(j=0;j<Ids->nind_tot;++j){
      if(Ids->include_ind_vec[j] && strcmp(Donors->donorlabels[i],Ids->popid[j])==0){
	Donors->ndonorhaps_vec[i]+=2-Par->haploid_ind;
	Donors->ndonorhaps+=2-Par->haploid_ind;
	if(Par->vverbose)  fprintf(Par->out,"Donors::readDonorPops Assigning ind %s to donor population %s\n",Ids->id[j], Donors->donorlabels[i]);
      }
    }
  }

  // assign haplotypes to recipients.
  // First count the number of recipients
  Donors->nrecipinds=0;
  if(usingdonorfile==0){ // we are doing things from Idfile2
    Donors->nrecipinds=Ids2->nind_inc;
  }else{
    for(i=0;i<Donors->nrecippops;++i){
      for(j=0;j<Ids->nind_tot;++j){
	if(Ids->include_ind_vec[j] && strcmp(Donors->recippoplabels[i],Ids->popid[j])==0){
	  ++Donors->nrecipinds;
	}
      }
    }
  }
  if(Par->vverbose)  fprintf(Par->out,"Donors::readDonorPops Counted %i recipient individuals\n",Donors->nrecipinds);
  
  // Assign memory for them
  Donors->reciplabels =malloc(Donors->nrecipinds * sizeof(char *));
  Donors->recippops = malloc(Donors->nrecipinds * sizeof(char *));
  for(i=0;i<Donors->nrecipinds;i++)
      Donors->reciplabels[i]=malloc(2047*sizeof(char));

  // And then actually assign them 
  int assignedrecips=0;
  for(i=0;i<Donors->nrecipinds;i++){ // in all-vs-all mode, every ind is a population
    Donors->recippops[i]= -1;
  } // No population for undefined recipient individuals
  if(usingdonorfile==0){
    int jj=0;
    for(j=0;j<Ids2->nind_tot;++j){
      if(Ids2->include_ind_vec[j]){
	if(Par->vverbose)   fprintf(Par->out,"Donors::readDonorPops Assigning ind %s to recipient population %s\n",Ids2->id[j], Donors->recippoplabels[0]);
	strcpy(Donors->reciplabels[assignedrecips],Ids2->id[j]);
	Donors->recippops[jj] = 0;
	jj++;
	++assignedrecips;
      }// end if ids included
    }// end for over j
  }else{
    int jj=0;
    for(i=0;i<Donors->nrecippops;++i){
      for(j=0;j<Ids->nind_tot;++j){
	if(Ids->include_ind_vec[j]){
	  if(strcmp(Donors->recippoplabels[i],Ids->popid[j])==0) {
	    if(Par->vverbose)   fprintf(Par->out,"Donors::readDonorPops Assigning ind %s to recipient population %s\n",Ids->id[j], Donors->recippoplabels[i]);
	    strcpy(Donors->reciplabels[assignedrecips],Ids->id[j]);
	    Donors->recippops[jj] = i;
	    jj++;
	    ++assignedrecips;
	  }
	}// end if ids included
      }// end for over j
    }// end for over i
  }// end else
  //  Data->nhaps_startpop=assignedrecips*Data->hapsperind;
  
  closeDonorlist(Infiles);
  /// Finished reading file
  //  if(Par->start_val==0) Par->start_val=1;
  if(Par->end_val==0) Par->end_val=Donors->nrecipinds;

  reallocateDataHaps(Data,Donors,Ids,Par);// Create the correct number of condition an recipient haplotypes, now we know the donor structure
  

  // Now we sanity check and tidy up the populations
  if (Par->prior_donor_probs_ind==1) {// we have donor probs; do they check out?
    if (((totaldonorprobs > 1.00001) || (totaldonorprobs < 0.99999)) )
    {
      fprintf(Par->out,"Donors::readDonorPops Probabilities across all donors in %s does not sum to 1.0 (instead sums to %lf)! Exiting....\n",Infiles->donorlist,totaldonorprobs);
      stop_on_error(1,Par->errormode,Par->err);
    }
  }// end if we hve donor probs

  if(Par->vverbose) {
    fprintf(Par->out,"Donors::readDonorPops Created %i donor populations containing %i haplotypes:\n",Donors->ndonorpops,Donors->ndonorhaps);
  }
  for(i=0;i<Donors->ndonorpops;i++){
    fprintf(Par->out,"  %i: %s, %i\n",i+1,Donors->donorlabels[i],Donors->ndonorhaps_vec[i]);
  }
  Donors->ndonorvecsize = Donors->ndonorpops;
  //  Donors->nrecipinds = recipon;
}

void populateDonors(struct donor_t *Donors,struct ids_t *Ids,struct ids_t *Ids2,struct data_t *Data,struct param_t *Par) {
  // Create donor populations in all-vs-all mode
  if (!Par->allvsall){
    fprintf(Par->out,"Logic error: we should not create donor labels when not in all-vs-all mode!\n"); stop_on_error(1,Par->errormode,Par->err);
  }
  int i;
  
  Donors->ndonorpops=(Data->condhaps)/Data->hapsperind;
  if(Par->geno2_find){
    if(Par->vverbose) fprintf(Par->out,"PopulateDonors: creating pseudo donor file...\n");
    Donors->nrecipinds=(Data->reciphaps)/Data->hapsperind;
    Donors->ndonorvecsize=Donors->ndonorpops;
    Donors->nrecippops=1;
    if(Par->vverbose) fprintf(Par->out,"PopulateDonors: Using %i recipient individuals, %i donor populations\n",Donors->nrecipinds,Donors->ndonorpops);
  }else{
    Donors->ndonorvecsize=Donors->ndonorpops;
    Donors->nrecipinds=Donors->ndonorpops;
    if(Donors->nrecipinds != Ids->nind_inc){
      fprintf(Par->out,"Mismatch between the number of individuals in the data (%d) and the ID file (%d); terminating!\n",Donors->nrecipinds,Ids->nind_inc);
      stop_on_error(1,Par->errormode,Par->err);
    }
    Donors->nrecippops=Donors->nrecipinds;
  }

  if(Par->vverbose) fprintf(Par->out,"Assigning memory for %i recipients, %i donorlabels\n",Donors->nrecipinds,Donors->ndonorvecsize);
  Donors->reciplabels=malloc(Donors->nrecipinds * sizeof(char *));
  Donors->recippoplabels=malloc(Donors->nrecippops * sizeof(char *));
  Donors->donorlabels=malloc(Donors->ndonorvecsize * sizeof(char *));
  Donors->ndonorhaps_vec = malloc(Donors->ndonorvecsize * sizeof(int));
  Donors->recippops = malloc(Donors->nrecipinds * sizeof(char *));

  for(i=0;i<Donors->nrecipinds;i++){ // in all-vs-all mode, every ind is a population
    Donors->recippops[i]=i;
  }
  // alocate memory
  // NOTE: the self population is a donor here
  for(i=0;i<Donors->ndonorvecsize;i++) {
    Donors->donorlabels[i]=malloc(2047*sizeof(char));
  }
  for(i=0;i<Donors->nrecipinds;i++) {
    Donors->reciplabels[i]=malloc(2047*sizeof(char));
  }
  for(i=0;i<Donors->nrecippops;i++) {
    Donors->recippoplabels[i]=malloc(2047*sizeof(char));
  }
  int ii=-1;

  Donors->ndonorhaps=0;
  for(i=0;i<Donors->ndonorvecsize;i++){
    Donors->ndonorhaps+=Data->hapsperind;
    Donors->ndonorhaps_vec[i]=Data->hapsperind;
    do{
      ++ii;
    }while(Ids->include_ind_vec[ii]==0);
    strcpy(Donors->donorlabels[i],Ids->id[ii]);
    if(!Par->geno2_find) {
      strcpy(Donors->reciplabels[i],Ids->id[ii]);
      strcpy(Donors->recippoplabels[i],Ids->id[ii]);
    }
    //    sprintf(Donors->donorlabels[i],"IND%i",i+1);
    //    sprintf(Donors->reciplabels[i],"IND%i",i+1);
    if(Par->vverbose) fprintf(Par->out,"populateDonors: Assigned donor pop %i the name %s with %i haplotype(s)\n",i+1,Donors->donorlabels[i],Donors->ndonorhaps_vec[i]);
    fflush(Par->out);
  }
  // set up recipients
  if(Par->vverbose) fprintf(Par->out,"populateDonors: Assigning %i recipient individuals\n",Donors->nrecipinds);
  if(Par->geno2_find){
    for(i=0;i<Donors->nrecipinds;i++){
      do{
	++ii;
      }while(Ids->include_ind_vec[ii]==0);
      strcpy(Donors->reciplabels[i],Ids->id[ii]);
      if(Par->vverbose) fprintf(Par->out,"populateDonors: Assigned donor ind %i the name %s\n",i+1,Donors->reciplabels[i]);
    } 
    reallocateDataHaps(Data,Donors,Ids,Par);
  }

  // assign recipient haplotypes
  if(Par->geno2_find){    
    ii=0;
    int k;
    for(i=0;i<Donors->nrecipinds;++i){
      for(k=0;k<Ids->nind_tot;++k){
	if(Ids->include_ind_vec[k] && strcmp(Ids->id[k],Donors->reciplabels[i])==0){
	  Data->recipnums[ii]=k;
	  ++ii;
	  if(Par->vverbose){
	    fprintf(Par->out,"Donors::reallocateDataHaps Assigned data structures for recipient individual %i, location %i in the ID file (called %s)\n",ii,k,Ids->id[k]);
	  }
	}
      }
    }
    Par->end_val=Par->start_val+Donors->nrecippops;
  }
  
  if(Par->vverbose) fprintf(Par->out,"populateDonors: Assigned  %i donor haplotype(s)\n",Donors->ndonorhaps);
}    




void clearDonors(struct donor_t *Donors,struct param_t *Par) 
{
  // NEED TO TIDY UP EVERYTHING HERE!
  
  if(Donors==NULL) return;
  //* for some reason this is buggy?
  int i;
  for(i=0;i<Donors->ndonorvecsize;i++)
    free(Donors->donorlabels[i]);
  for(i=0;i<Donors->nrecipinds;i++)
    free(Donors->reciplabels[i]);
  for(i=0;i<Donors->nrecippops;i++)
    free(Donors->recippoplabels[i]);

  free(Donors->donorlabels);
  free(Donors->reciplabels);
  free(Donors->recippoplabels);
  free(Donors->recippops);

  free(Donors->ndonorhaps_vec);
  if (Par->prior_donor_probs_ind==1)  free(Donors->ndonorprobs);
  if (Par->mutation_rate_ind==1)   free(Donors->ndonormutrates);
  free(Donors);
}
