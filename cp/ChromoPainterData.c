#include "ChromoPainterData.h"
#include "ChromoPainterMutEM.h"

//#define DEBUG

int codeHaplotypes(char in){
  if(in=='0')
    return(0);
  if(in=='1')
    return(1);
  if(in=='A')
    return(2);
  if(in=='C')
    return(3);
  if(in=='G')
    return(4);
  if(in=='T')
    return(5);
  if(in=='a')
    return(2);
  if(in=='c')
    return(3);
  if(in=='g')
    return(4);
  if(in=='t')
    return(5);
  if(in=='-')
    return(8);
  if(in=='8')
    return(8);
  if(in=='9')
    return(9);
  if(in=='N')
    return(9);
  if(in=='?')
    return(9);
  if(in=='.')
    return(9);
  return(-1);
}

int * getallelic_type_count_vec(struct data_t * Data){
                  // find number of alleles per snp (this is NOT every used, but perhaps should be to get default mutation rate correct):
      int * allelic_type_count_vec = malloc(Data->nsnps*sizeof(int));
      int * found_vec = malloc(6*sizeof(int));
      int i=0,j=0;
      for (j=0; j < Data->nsnps; j++)
	{
	  allelic_type_count_vec[j] = 0;
	  for (i=0; i < 6; i++)
	    found_vec[i]=0;
	  for (i=0; i < Data->current_donor_nind*Data->hapsperind; i++)
	    {
	  if ((Data->cond_chromosomes[i][j] == 0) && (found_vec[0] == 0))
		{
		  allelic_type_count_vec[j] = allelic_type_count_vec[j] + 1;
		  found_vec[0] = 1;
		}
	      if ((Data->cond_chromosomes[i][j] == 1) && (found_vec[1] == 0))
		{
		  allelic_type_count_vec[j] = allelic_type_count_vec[j] + 1;
		  found_vec[1] = 1;
		}
	      if ((Data->cond_chromosomes[i][j] == 2) && (found_vec[2] == 0))
		{
		  allelic_type_count_vec[j] = allelic_type_count_vec[j] + 1;
		  found_vec[2] = 1;
		}
	      if ((Data->cond_chromosomes[i][j] == 3) && (found_vec[3] == 0))
		{
		  allelic_type_count_vec[j] = allelic_type_count_vec[j] + 1;
		  found_vec[3] = 1;
	}
	      if ((Data->cond_chromosomes[i][j] == 4) && (found_vec[4] == 0))
		{
		  allelic_type_count_vec[j] = allelic_type_count_vec[j] + 1;
		  found_vec[4] = 1;
		}
	      if ((Data->cond_chromosomes[i][j] == 5) && (found_vec[5] == 0))
		{
		  allelic_type_count_vec[j] = allelic_type_count_vec[j] + 1;
		  found_vec[5] = 1;
		}
	    }
	}
      free(found_vec);
      return(allelic_type_count_vec);
}

void setIndAsRecipient(int m, struct data_t *Data,struct ids_t *Ids, struct param_t *Par)
{
  int h,hh,i=0,ii=0;
  free(Data->current_donor_inds);
  free(Data->current_donor_haps);
  // assign recipient
  Data->currentind=Data->recipnums[m];
  for(h=0;h<Data->hapsperind;h++){
    int mm;
    mm=(Data->recipnums[m]*Data->hapsperind)+h;
    Data->ind_chromosomes[h]=Data->all_chromosomes[mm];
    
    if(Par->vverbose){
      fprintf(Par->out,"Data::setIndAsRecipient: Will be painting individual %i, haplotype %i\n",m,mm);
      for(i=0;i<10;++i)fprintf(Par->out,"%i",Data->ind_chromosomes[h][i]);
      fprintf(Par->out,"\n");
    }
  }

  // assign donors
  Data->current_donor_nind = 0;
  for(i=0;i<Data->condhaps/Data->hapsperind;++i){
    if(Data->condnums[i]!=Data->currentind) Data->current_donor_nind++;
  }
  Data->current_donor_inds=malloc(Data->current_donor_nind*sizeof(int));
  Data->current_donor_haps=malloc(Data->current_donor_nind*Data->hapsperind*sizeof(int));

  ii=0; hh=0;
  for(i=0;i<Data->condhaps/Data->hapsperind;++i){
    if(Data->condnums[i]!=Data->currentind) {
      Data->current_donor_inds[ii++] = Data->condnums[i];
      for(h=0;h<Data->hapsperind;++h){
	Data->current_donor_haps[hh] = Data->condnums[i]*Data->hapsperind + h;
	Data->cond_chromosomes[hh] = Data->all_chromosomes[Data->current_donor_haps[hh]];
	++hh;
      }
    }else{
      if(Par->vverbose) fprintf(Par->out,"Omitting individual indexed %d (individual called %s at position %i in the ID file) from the donor pool\n",i,Ids->id[Data->currentind],Data->currentind);
    }
  }

}

struct ids_t *CreateIdsFromPhase(char *filename, struct param_t *Par) {
  // fd is the PHASE file
  char line1[100000];
  char line2[100000];
  char line3[100000];
  struct ids_t *Ids = malloc(sizeof(struct ids_t));
  int i=0;

  FILE *fd=fopen(filename,"r"); // phase file, remember
  // First line: EITHER ndonors or ninds
  if(fgets(line1,2047,fd) ==NULL) {// get the whole line
    fprintf(Par->out,"error reading %s\n",filename); stop_on_error(1,Par->errormode,Par->err);	
  }
  if(fgets(line2,2047,fd) ==NULL) {// get the whole line
    fprintf(Par->out,"error reading %s\n",filename); stop_on_error(1,Par->errormode,Par->err);	
  }
  if(fgets(line3,2047,fd) ==NULL) {// get the whole line
    fprintf(Par->out,"error reading %s\n",filename); stop_on_error(1,Par->errormode,Par->err);	
  }
  fclose(fd);
  if(line3[0]=='P'){ // New style datafile without the initial "donor haplotypes" line
    sscanf(line1,"%d",&Ids->nind_tot);// read the number fof HAPLOTYPES rom the first line
    Ids->nind_tot/=(2-Par->haploid_ind);
  }else{
    sscanf(line2,"%d",&Ids->nind_tot);// read the number from the second line
  }
  Ids->nind_inc=Ids->nind_tot;

  // allocate memory
  Ids->include_ind_vec = malloc(Ids->nind_tot*sizeof(int));
  Ids->id = malloc(Ids->nind_tot*sizeof(char *));
  Ids->popid = malloc(Ids->nind_tot*sizeof(char *));
  for(i=0;i<Ids->nind_tot;i++) {
    Ids->id[i]=malloc(1000*sizeof(char));
    Ids->popid[i]=malloc(1000*sizeof(char));
  }

  // create the id data

  for(i=0;i<Ids->nind_tot;i++) {
    Ids->include_ind_vec[i]=1;
    sprintf(Ids->id[i],"IND%d",i+1);
    strcpy(Ids->popid[i],Ids->id[i]);
    Ids->include_ind_vec[i]=1;
  }
  return(Ids);
}

struct ids_t *ReadIdfile(char *filename, struct param_t *Par) {
  // fd is the ID file
  struct ids_t *Ids = malloc(sizeof(struct ids_t));
  char line[10000];
  char *step;
  char waste[1000];
  FILE *fd;
  int i=0;
  // Count the number of entries in the file
  fd= fopen(filename,"r");
  if (fd == NULL) { fprintf(Par->out,"error opening %s\n",filename); stop_on_error(1,Par->errormode,Par->err);}
  Ids->nind_tot=0;
  Ids->nind_inc=0;
  while(!feof(fd))
    {
      if (fgets(line,2047,fd)!=NULL)
	++ Ids->nind_tot;
    }
  fclose(fd);
  if(Par->vverbose)fprintf(Par->out,"ReadIds: Detected %d IDs\n",Ids->nind_tot);
  // allocate memory
  Ids->include_ind_vec = malloc(Ids->nind_tot*sizeof(int));
  Ids->id = malloc(Ids->nind_tot*sizeof(char *));
  Ids->popid = malloc(Ids->nind_tot*sizeof(char *));
  for(i=0;i<Ids->nind_tot;i++) {
    Ids->id[i]=malloc(1000*sizeof(char));
    Ids->popid[i]=malloc(1000*sizeof(char));
  }
  // Read the names and inclusion status of each ind
  fd= fopen(filename,"r");
  if (fd == NULL) { fprintf(Par->out,"error opening %s\n",filename); stop_on_error(1,Par->errormode,Par->err);}
  for(i=0;i<Ids->nind_tot;i++){
      if(fgets(line,2047,fd) ==NULL) {// get the whole line
	fprintf(Par->out,"error reading %s\n",filename); stop_on_error(1,Par->errormode,Par->err);	
      }
      step=line; 
      if(reading(&step,"%s",waste)) {
	if(strlen(waste)>1000){
	  fprintf(Par->out,"error reading %s: id %i too long?\n",filename,i);
	  stop_on_error(1,Par->errormode,Par->err);
	}
	strcpy(Ids->id[i],waste);
      }else {fprintf(Par->out,"Invalid line %d in ID file %s\n",i,filename);stop_on_error(1,Par->errormode,Par->err);}
      if(reading(&step,"%s",waste)) {// there is a pop ID column
	strcpy(Ids->popid[i],waste);
      }else{ // no more fields
	strcpy(Ids->popid[i],Ids->id[i]);
      }
      Ids->include_ind_vec[i]=1;
      if(reading(&step,"%s",waste)) {// there is an inclusion column
	if (strcmp(waste,"0")==0){
	  Ids->include_ind_vec[i]=0;
	}
      }
      if(Ids->include_ind_vec[i]==1){
	++ Ids->nind_inc;	
      }
      if(Par->vverbose)fprintf(Par->out,"ReadIds: ID #%i Label %s Pop %s including %i\n",i,Ids->id[i],Ids->popid[i],Ids->include_ind_vec[i]);

    }
  fclose(fd);
  return(Ids);
}

struct ids_t *CreateMergedIds(struct ids_t *Ids1, struct ids_t *Ids2, struct param_t *Par){
  struct ids_t *Ids = malloc(sizeof(struct ids_t));
  int i;
  Ids->nind_tot = Ids1->nind_tot + Ids2->nind_tot;
  Ids->nind_inc = Ids1->nind_inc + Ids2->nind_inc;

  if(Par->vverbose) fprintf(Par->out,"Assigning space for %i IDs (retaining %i)\n",Ids->nind_tot,Ids->nind_inc);
  //
  Ids->include_ind_vec = malloc(Ids->nind_tot*sizeof(int));
  Ids->id = malloc(Ids->nind_tot*sizeof(char *));
  Ids->popid = malloc(Ids->nind_tot*sizeof(char *));
  for(i=0;i<Ids->nind_tot;i++) {
    Ids->id[i]=malloc(1000*sizeof(char));
    Ids->popid[i]=malloc(1000*sizeof(char));
  }
  
  for(i=0;i<Ids1->nind_tot;i++) {
    Ids->include_ind_vec[i] = Ids1->include_ind_vec[i];
    strcpy(Ids->id[i],Ids1->id[i]);
    strcpy(Ids->popid[i],Ids1->popid[i]);
  }
  
  int ii=Ids1->nind_tot;
  for(i=0;i<Ids2->nind_tot;i++) {
    Ids->include_ind_vec[ii] = Ids2->include_ind_vec[i];
    strcpy(Ids->id[ii],Ids2->id[i]);
    strcpy(Ids->popid[ii],Ids2->popid[i]);
    ++ii;
  }
  return(Ids);
}

struct data_t *ReadDataHeader(FILE *fd, struct param_t *Par){
  struct data_t *Data;
  
  char line1[10000];
  char line2[10000];
  char * line = malloc(100000000 * sizeof(char));
  char *step;
  char waste[10];
  int i,j;
  
  char *ret;

  Data=malloc(sizeof(struct data_t));
  Data->hapsperind = 2-Par->haploid_ind; // number of haploids per individual

  if(Par->vverbose) fprintf(Par->out,"ReadData: Starting read of header\n");

  //////////////////////////
  // Read the header

  // Is this the hap line or the ind line?
  if(fgets(line1,2047,fd)==NULL) { fprintf(Par->out,"ReadData: error with PHASE-style input file; first (haplotype donor) line empty\n"); stop_on_error(1,Par->errormode,Par->err);}

  // Ind line or SNP line?
  if(fgets(line2,2047,fd)==NULL) { fprintf(Par->out,"ReadData: error with PHASE-style input file; second (individual) line empty\n"); stop_on_error(1,Par->errormode,Par->err);}

  // SNP line or position line?
  if(fgets(line,100000000,fd)==NULL) { fprintf(Par->out,"ReadData: error with PHASE-style input file; third (snps() line empty\n"); stop_on_error(1,Par->errormode,Par->err);}

  // Figure out whether we have the header line
  if(line[0]!='P'){ // Old style format
    sscanf(line1,"%d",&Data->nhaps_startpop);// read the number
    sscanf(line2,"%lf",&Data->nind);
    sscanf(line,"%d",&Data->nsnps_raw);
    Data->nhapstotal=Data->nind * Data->hapsperind;
    if(fgets(line,100000000,fd) ==NULL) { fprintf(Par->out,"ReadData: error with PHASE-style input file; fourth (SNP location) line empty\n"); stop_on_error(1,Par->errormode,Par->err);}
  }else{ // New style format
    Data->nhaps_startpop=0;
    sscanf(line1,"%d",&Data->nhapstotal);
    Data->nind = floor(Data->nhapstotal/Data->hapsperind);
    sscanf(line2,"%d",&Data->nsnps_raw);
  }

  //////////////////////////////
  // Check the header information
  // nhaps
  if (Data->nhaps_startpop < 0) { fprintf(Par->out,"Number of donor haplotypes must be >= 0. Found %d. Exiting...\n",Data->nhaps_startpop); stop_on_error(1,Par->errormode,Par->err);
  }else if ((Par->allvsall==0) && (Data->nhaps_startpop==0) && (Par->donorlist_find==0)) {
    fprintf(Par->out," Number of donor haplotypes can only be 0 if you choose  all-versus-all '-a' options)\n"); stop_on_error(1,Par->errormode,Par->err);
  }

  // nind
  if (Data->nind < 0) { fprintf(Par->out,"Number of total individuals must be > 0. Exiting...\n"); stop_on_error(1,Par->errormode,Par->err);}
  // No check here, due to the possibility of Ids being made up from two files
  // Check is in the donors section
  /* if(Data->nind != Ids->nind_tot){ */
  /*   fprintf(Par->out,"Mismatch between the number of individuals in the data (%f) and the ID file (%d); terminating!\n",Data->nind,Ids->nind_inc); */
  /*   stop_on_error(1,Par->errormode,Par->err); */
  /* } */

  // nsnps
  Data->nsnps=Data->nsnps_raw;
  if (Data->nsnps_raw <= 0) { fprintf(Par->out,"Number of sites must be > 0. Exiting...\n"); stop_on_error(1,Par->errormode,Par->err);}
  if(Par->emloci>0) {  /// choose a random region
    if(Data->nsnps_raw>Par->emloci){
      Par->startlocus=randInt(0,Data->nsnps_raw-Par->emloci);
      Par->endlocus =Par->startlocus + Par->emloci;
    }
    Par->emloci=-1;
    if(Par->verbose) fprintf(Par->out,"Choosing the random region %d-%d (of %d total loci)\n",Par->startlocus,Par->endlocus,Data->nsnps_raw);
  }
  if(Par->endlocus<=0) Par->endlocus=Data->nsnps_raw;
  if(Par->startlocus>0 || Par->endlocus<Data->nsnps_raw){ /// Restrict the data to a region
    if(Par->endlocus>Data->nsnps_raw){
      fprintf(Par->out,"Trying to use the final locus %d which is greater than the number of snps %d. Exiting\n",Par->endlocus,Data->nsnps_raw);
      stop_on_error(1,Par->errormode,Par->err);
    }
    Data->nsnps = Par->endlocus - Par->startlocus;
  }

  if(Par->vverbose) fprintf(Par->out,"Found %d SNPs in phase file, using %d of them\n",Data->nsnps_raw,Data->nsnps);

  // This is only true for all-vs-all mode; in donor mode we will reallocate after seeing the donors

  if(Par->vverbose) fprintf(Par->out,"ReadData: Reading %d positions\n",Data->nsnps);

  /* Positions */
  Data->positions=malloc(Data->nsnps*sizeof(double));
  Data->lambda=malloc((Data->nsnps-1)*sizeof(double));

  step=line;
  reading(&step,"%s",waste);
  for (i=0; i<Data->nsnps_raw; i++)
    {
      int ii=i-Par->startlocus;
      //      fprintf(Par->out,"i=%d, ii=%d sl=%d el=%d\n",i,ii,Par->startlocus,Par->endlocus);
      double tmp;
      reading(&step,"%lf",&tmp);

      if(i<Par->startlocus || i>=Par->endlocus) continue;
      if(Par->vverbose && ii<10) fprintf(Par->out,"Using data location at position %i in the file, %i SNP, at location %f\n",i,ii,tmp);
      Data->positions[ii]=tmp;
      if (Data->positions[ii] < 0) { fprintf(Par->out,"Basepair positions must be >= 0. Exiting...\n"); stop_on_error(1,Par->errormode,Par->err);}
      if(ii>0) if (Data->positions[ii] == Data->positions[ii-1])
	{
	  if(Par->jitter_locations==0) { 
	    fprintf(Par->out,"Basepair positions must be increasing (at SNPs %i-%i with positions %f-%f). Rerun with option \"-J\" to ignore. Exiting...\n",i-1,i,Data->positions[i-1],Data->positions[i]); stop_on_error(1,Par->errormode,Par->err);
	  }else{
	    double newloc=Data->positions[ii-1]+1;
	    Data->positions[ii]=newloc;
	  }
	}
      if (ii < (Data->nsnps-1)) Data->lambda[ii] = 1.0;
    }

  free(line);
  if(Par->vverbose) fprintf(Par->out, "ReadData: Successfully read header.\n");
  return(Data);
}

struct data_t *ReadData(struct infiles_t *Infiles, struct ids_t *Ids,struct param_t *Par){
  struct data_t *Datamain,*Datarecip,*Data;

  char line1[10000];
  char line2[10000];
  char * line = malloc(100000000 * sizeof(char));
  char *step;
  //  char waste[10];
  int i,j;
  int ninds_main=0,ninds_recip=0,ninds_switch=0;
  FILE *fphase=NULL, *fphase2=NULL;
  char *ret;

  if(Par->vverbose)fprintf(Par->out,"ReadData: Reading main phase file %s\n",Infiles->phase);
  //  if(!openPhase(Infiles,Par->out)){
  fphase=fopen(Infiles->phase,"r"); // phase file. This is LOCAL TO THIS FUNCTION so NOT PART OF Infiles!
  if(fphase==NULL) { 
    fprintf(Par->out,"error opening file %s\n",Infiles->phase); stop_on_error(1,Par->errormode,Par->err);
  }
  Datamain=ReadDataHeader(fphase,Par); // Information from the main data file

  if(Par->geno2_find) {
    // read the recipient file
    if(Par->vverbose)fprintf(Par->out,"ReadData: Reading recipient phase file %s\n",Infiles->phase2);
    fphase2=fopen(Infiles->phase2,"r"); // phase file, remember
    Datarecip=ReadDataHeader(fphase2,Par); // Information from the recipient data file
    // Count the number of retained mainfile individuals and recipient file individuals
    for(i=0;i<Datamain->nhapstotal/Datamain->hapsperind;++i) {
      ninds_main += Ids->include_ind_vec[i];
      ninds_switch++;
    }
    for(i=Datamain->nhapstotal/Datamain->hapsperind;i<Ids->nind_tot;++i) ninds_recip += Ids->include_ind_vec[i];
    // Create the merged data
    Data=malloc(sizeof(struct data_t));
    Data->positions=Datamain->positions;
    Data->lambda=Datamain->lambda;
    Data->hapsperind=Datamain->hapsperind;
    Data->nhaps_startpop = Datamain->nhaps_startpop + Datarecip->nhaps_startpop; // *** NOT SURE ABOUT THIS
    Data->nind = Datamain->nind + Datarecip->nind;
    Data->nhapstotal = Datamain->nhapstotal + Datarecip->nhapstotal;
    if(Datamain->nsnps_raw !=Datarecip->nsnps_raw){
      fprintf(Par->out,"ReadData: error matching main and recipient phase input files SNP difference (main: %i vs recipient %i)\n",Datamain->nsnps_raw ,Datarecip->nsnps_raw); stop_on_error(1,Par->errormode,Par->err);
    }
    Data->nsnps_raw = Datamain->nsnps_raw;
    Data->nsnps = Datamain->nsnps;
    // Update quantities that depend on the merging
    Data->nhaps = Ids->nind_inc * Data->hapsperind;  // h*N
    Data->condhaps = 0;
    for(i=0;i<Datamain->nhapstotal;++i){
      if(Ids->include_ind_vec[(int)floor(i/Data->hapsperind)]) ++Data->condhaps;
    }
    //    Data->condhaps = Data->nhaps;
    Data->reciphaps = 0;
    for(i=Datamain->nhapstotal;i<Data->nhapstotal;++i){
      if(Ids->include_ind_vec[(int)floor((i+Datamain->nhapstotal)/Data->hapsperind)]) ++Data->reciphaps;
    }
    // remove the data objects we don't need
    free(Datarecip->positions);
    free(Datarecip->lambda);
    free(Datarecip);
    free(Datamain);
    if(Data->nhapstotal !=Ids->nind_tot*Data->hapsperind){ fprintf(Par->out,"Number of total haps (in phase files %s + %s) must be equal to the number of individuals in the ID files (%s + %s) multiplied by the ploidy. You have specified %i haplotypes in the data, with %i IDs with a ploidy of %i\n",Infiles->phase,Infiles->phase2,Infiles->id,Infiles->id2,Data->nhapstotal,Ids->nind_tot,Data->hapsperind);  stop_on_error(1,Par->errormode,Par->err);}
  }else{
    if(Datamain->nhapstotal !=Ids->nind_tot*Datamain->hapsperind){ fprintf(Par->out,"Number of total haps (in phase file %s) must be equal to the number of individuals in the ID file (%s) multiplied by the ploidy. You have specified %i haplotypes in the data, with %i IDs with a ploidy of %i\n",Infiles->phase,Infiles->id,Datamain->nhapstotal,Ids->nind_tot,Datamain->hapsperind); stop_on_error(1,Par->errormode,Par->err);}
    Data=Datamain;
  // Compute the number of haplotypes we will use
    Data->nind = Ids->nind_inc; /// Set this to the number BEING USED
    //    Data->nhapstotal = Data->nind*Data->hapsperind;
    Data->nhaps = Data->hapsperind*Data->nind;  // h*N
    Data->condhaps = Data->nhaps;
    Data->reciphaps = Data->nhaps;
  }
    
  if ((Data->nhaps-Data->hapsperind) < Data->nhaps_startpop) { fprintf(Par->out,"Number of total haps must be greater than or equal to number of donor haplotypes plus one extra individual.\n"); stop_on_error(1,Par->errormode,Par->err);}

  /////////////////////////////
  // Memory allocation

  if(Par->vverbose) {
    fprintf(Par->out,"ReadData: Allocating memory: %i\n",			    Data->hapsperind);
    fprintf(Par->out,"ReadData: Allocating memory: %i\n",			    Data->condhaps);
    fprintf(Par->out,"ReadData: Allocating memory: %i\n",			    Data->nhapstotal);
    fprintf(Par->out,"ReadData: Allocating memory: %i\n",			    Data->nhaps);

    fprintf(Par->out,"ReadData: Allocating memory: %i recipient, %i condition haplotypes, %i in file, using %i\n",
			    Data->hapsperind,
			    Data->condhaps,
			    Data->nhapstotal,
			    Data->nhaps);
  }
  Data->all_chromosomes=malloc(Data->nhapstotal*sizeof(int *));
  Data->cond_chromosomes=malloc(Data->condhaps*sizeof(int *));
  Data->ind_chromosomes=malloc(Data->hapsperind*sizeof(int *));

  // Condition and recipient individuals
  Data->condnums = malloc((Data->condhaps)/Data->hapsperind*sizeof(int));
  Data->recipnums = malloc((Data->reciphaps)/Data->hapsperind*sizeof(int));

  // We allocate memory for these, but they get reallocated in the setIndAsRecipient call
  Data->current_donor_nind=0;
  Data->current_donor_inds=malloc(Data->current_donor_nind*sizeof(int));
  Data->current_donor_haps=malloc(Data->current_donor_nind*Data->hapsperind*sizeof(int));
  
  j=0;
  if(Par->geno2_find){
    for (i=0; i<Data->nhapstotal/Data->hapsperind; i++){
      if(Ids->include_ind_vec[i]) {
	if(j<Data->condhaps/Data->hapsperind) {
	  Data->condnums[j++] = i;
	}else{
	  Data->recipnums[j-Data->condhaps/Data->hapsperind] = i;
	  ++j;
	}
      }
    }
  }else{
    for (i=0; i<Data->nhapstotal/Data->hapsperind; i++){
      if(Ids->include_ind_vec[i]) {
	Data->condnums[j] = i;
	Data->recipnums[j++] = i;
      }
    }
  }
  // Allocate the all_chromosomes memory
  // Note that cond and ind just point to this
  
  for (i=0; i<Data->nhapstotal; i++){
    if(Ids->include_ind_vec[(int)floor(i/Data->hapsperind)]) {
      Data->all_chromosomes[i]=malloc(Data->nsnps*sizeof(int));
    }else Data->all_chromosomes[i]=malloc(0*sizeof(int));// malloc for easy free later
  }


  ///////////////////////
  // Read the Data
  if(Par->vverbose) fprintf(Par->out,"ReadData: Reading %i haplotype lines\n",Data->nhapstotal);

  int ii=-1;
  FILE *fd=fphase;
  for (i= 0; i < Data->nhapstotal; i++)
    {
      if(Par->geno2_find && i==ninds_switch*Data->hapsperind) {
	fd=fphase2;
	if(Par->vverbose) fprintf(Par->out,"(switching to recipient file %s after %d haplotypes)\n",Infiles->phase2,i);
      }
      if(fgets(line,100000000,fd) ==NULL) { fprintf(Par->out,"ReadData: error with PHASE-style input file at haplotype %i of %i\n",i, Data->nhapstotal); stop_on_error(1,Par->errormode,Par->err);}
      
      if(Ids->include_ind_vec[(int)floor(i/Data->hapsperind)]==0) continue; // skip individuals which are ignored
      ++ii;
      if(Par->vverbose) fprintf(Par->out,"Retaining individual %d as haplotype %i  (%d after removing unused haps) for ID %s ",(int)floor(i/Data->hapsperind),i,ii,Ids->id[(int)floor(i/Data->hapsperind)]);
      if(strlen(line)!=(Data->nsnps_raw+1)) { fprintf(Par->out,"ReadData: error with PHASE-style input file at haplotype %i of %i : Incorrect number of SNPS (expected %d, received %d)\n",i, Data->nhapstotal,Data->nsnps_raw,(int) strlen(line)-1); stop_on_error(1,Par->errormode,Par->err);}

      for (j=0; j<Data->nsnps_raw;j++) {
	if(j<(Par->startlocus) || j>=Par->endlocus) continue;
	int jj=j-Par->startlocus;
	if(jj>Data->nsnps){
	  fprintf(Par->out,"Trying to read the %dth snp, yet we can only have %d. Exiting...\n",jj,Data->nsnps);
	  stop_on_error(1,Par->errormode,Par->err);
	}
	Data->all_chromosomes[i][jj]=codeHaplotypes(line[j]);
	if(Data->all_chromosomes[i][jj]<0) {
	  fprintf(Par->out,"Allele-type \"%c\" invalid for hap%d, snp%d. Exiting...\n",line[j],i+1,j+1);
	  stop_on_error(1,Par->errormode,Par->err);
	}
      }
      if(Par->vverbose){ fprintf(Par->out," : ");
	  int t;
	  	for(t=0;t<10;++t)fprintf(Par->out,"%i",Data->all_chromosomes[i][t]);
	fprintf(Par->out," ...\n");
	}
    }
  
  ///////////////////////
  // Assign the cond and ind haplotypes
  if(Par->vverbose) fprintf(Par->out,"ReadData: Assigning haplotypes\n");
  //setIndAsRecipient(0,Data,Par);

  free(line);

  // True for allvsall, otherwise we fix it later
  if(!Par->geno2_find){
    Data->reciphaps = Data->nind * Data->hapsperind - Data->nhaps_startpop;
  }
  Data->recipinds = Data->reciphaps/Data->hapsperind;

  // Information about the haplotyes
  if(Par->vverbose && Par->allvsall) {
    fprintf(Par->out,"Painting haplotypes ");
    for(i=0;i<Data->hapsperind;i++) fprintf(Par->out,"%i, ",i);
    fprintf(Par->out,"against all other haplotypes\n");
  }
  fclose(fphase);
  if(fphase2!=NULL) fclose(fphase2);
  /* closePhase(Infiles); */
  /* closePhase2(Infiles); */
  return Data;
}


void DestroyData(struct data_t *Data)
{
  if(Data==NULL) return;
  int i;
  for (i=0; i<Data->nhapstotal; i++)
    free(Data->all_chromosomes[i]);
  // remember: cond and ind just point to entries in all. Also, we never assigned space for the unused data, but we did malloc its pointer.

  free(Data->cond_chromosomes);
  free(Data->ind_chromosomes);
  free(Data->all_chromosomes);
  free(Data->positions);
  free(Data->lambda);
  free(Data->condnums);
  free(Data->recipnums);
  free(Data->current_donor_inds);
  free(Data->current_donor_haps);
  free(Data);
}


void DestroyIds(struct ids_t *Ids)
{
  int i;
  for (i=0; i<Ids->nind_tot; i++) {
    free(Ids->id[i]);
    free(Ids->popid[i]);
  }
  free(Ids->id);
  free(Ids->popid);
  free(Ids->include_ind_vec);
  free(Ids);
}
