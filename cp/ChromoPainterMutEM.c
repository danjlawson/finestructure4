
#include <setjmp.h>
#include "ChromoPainterMutEM.h"
#include "ChromoPainterConstants.h"



/*****************************************************************
// ChromoPainter algorithm for painting halpotypes

// to compile:  ./configure && make
// (old: gcc -Wall -o ChromoPainter ChromoPainter.c -lm -lz )
// To get detailed help, see ./ChromoPainter -h
// or refer to the helpfilestring in ChromoPainterConstants.h

*******************************************************************/

void usage(FILE *mainout)
{
  fprintf(mainout,"%s\n",helpfilestring);
}

void initRNG(struct param_t *Par)
{
  if(Par->rseed<0) {
    srand((unsigned)time(NULL));
  }else{
    srand((unsigned)Par->rseed);
  }
}

int randInt(int min,int max)
{ ///returns a random integer in the range [min,max], i.e. inclusive of min and inclusive of max
  
  return(min + (int) (((double) rand()/(RAND_MAX+1.0))*(1+max-min)));
}

void validateData(struct infiles_t *Infiles,struct param_t *Par,struct data_t *Data) {
  if (Par->start_val >= Data->recipinds)
    {
      fprintf(Par->out,"Your '-a' switch specifies to start with individual %d but there are only %d individuals in %s. Exiting....\n",Par->start_val+1,Data->recipinds,Infiles->phase);
      stop_on_error(1,Par->errormode,Par->err);
    }
  if ((Data->nhaps_startpop != 0) && (Par->donorlist_find==0))
    {
      fprintf(Par->out,"You have specified >0 donor haplotypes but have no donor-list (-f) file. Exiting....\n");
      stop_on_error(1,Par->errormode,Par->err);
    }
  if ((Par->allvsall) && (Par->donorlist_find))
    {
      fprintf(Par->out,"You have specified all-vs-all (-a) and a donor file (-f) which is not allowed.\n");
      stop_on_error(1,Par->errormode,Par->err);
    }
  if ((Data->nhaps_startpop != 0) && (Par->allvsall==1))
    {
      fprintf(Par->out,"You have specified >0 donor haplotypes but have also specified you want everyone to copy everyone ('-a')? If you want to use the '-a' switch, please make the first row of %s a 0. Exiting....\n",Infiles->phase);
      stop_on_error(1,Par->errormode,Par->err);
    }
  if ((Data->nhaps_startpop == 0) && (Par->allvsall==0) &&(Par->donorlist_find==0))
    {
      fprintf(Par->out,"You have specified 0 donor haplotypes, which is only allowed when specifying you want everyone to copy from everyone ('-a' switch) or you specify donors ('-f' switch). Exiting....\n");
      stop_on_error(1,Par->errormode,Par->err);
    }
  if ((Data->nhaps_startpop == 0) && (Par->donorlist_find==1))
    {
      //fprintf(Par->out,"You have specified 0 donor haplotypes but have also specified '-f' file. Are you sure you want 0 donor haplotypes?? Exiting....\n");
      //stop_on_error(1,Par->errormode,Par->err);
    }
  if (Data->nhaps_startpop >= Data->nhapstotal)
    {
      fprintf(Par->out,"You've specified wanting %d donor haps, but you have only %d total haps\n",Data->nhaps_startpop,Data->nhapstotal);
      stop_on_error(1,Par->errormode,Par->err);
    }
  if (((int) Data->nind) != Data->nind)
    {
      fprintf(Par->out,"You cannot have fractions of individuals, but you have specified %lf total individuals.\n",Data->nind);
      stop_on_error(1,Par->errormode,Par->err);
    }
  if ((floor(Data->reciphaps/2.0) != Data->reciphaps/2.0) && (Par->haploid_ind==0))
    {
      fprintf(Par->out,"You have specified %d recipient haplotypes, but this number must be even for diploid individuals. (Maybe you want to use '-j' to specify haploid individuals?)\n",Data->reciphaps);
      stop_on_error(1,Par->errormode,Par->err);
    }
}

void assignParameters(struct param_t *Par,struct infiles_t *Infiles,struct files_t *Files, int argc, char *argv[]) {
  /************************************************************/

  int i=0;
  for (i=1; i < argc; i++)
    {
    if (strcmp(argv[i],"--internalerrors")==0)
	{
	    Par->errormode=1;
	}
    if (strcmp(argv[i],"--noexitonerrors")==0)
	{
	    Par->errormode=2;
	}
    if ((strcmp(argv[i],"-help")==0) || (strcmp(argv[i],"--help")==0) || (strcmp(argv[i],"-h")==0))
	{
	  usage(Par->out);
	  stop_on_error(0,Par->errormode,Par->err);
	}
      if (strcmp(argv[i],"-g")==0)
	{
	  Par->geno_find=1;
	}
      if (strcmp(argv[i],"-r")==0)
	{
	  Par->recom_find=1;
	}
      if (strcmp(argv[i],"-f")==0)
	{
	  Par->donorlist_find=1;
	}
      if (strcmp(argv[i],"-i")==0)
	{
	  Par->EMiter_find=1;
	}
      if (strcmp(argv[i],"-s")==0)
	{
	  Par->numsamples_find=1;
	}
       if (strcmp(argv[i],"-n")==0)
	{
	  Par->ne_find=1;
	}
       if (strcmp(argv[i],"-M")==0)
	{
	  Par->mut_find=1;
	}
       if (strcmp(argv[i],"-k")==0)
	{
	  Par->region_size_find=1;
	}
       if (strcmp(argv[i],"-o")==0)
	{
	  Par->outfile_find=1;
	}
       if (strcmp(argv[i],"-t")==0)
	{
	  Par->idfile_find=1;
	}
      if (strcmp(argv[i],"-Rg")==0)
	{
	  Par->geno2_find=1;
	}
      if (strcmp(argv[i],"-Rt")==0)
	{
	  Par->idfile2_find=1;
	}
       if (strcmp(argv[i],"--emfilesonly")==0)
	 {
	   emoutfiles(Files);
	}
       if (strcmp(argv[i],"-ip")==0)
	 Par->copy_prop_em_find=1;
       if (strcmp(argv[i],"-in")==0)
	 Par->recom_em_find=1;
       if (strcmp(argv[i],"-im")==0)
	 Par->mutation_em_find=1;
       if (strcmp(argv[i],"-iM")==0)
	 Par->mutationALL_em_find=1;
       if (strcmp(argv[i],"-a")==0)
	 Par->allvsall=1;
       if (strcmp(argv[i],"-j")==0)
	 Par->haploid_ind=1;
       if (strcmp(argv[i],"-u")==0)
	 Par->unlinked_ind=1;
       if (strcmp(argv[i],"-p")==0)
	 Par->prior_donor_probs_ind=1;
       if (strcmp(argv[i],"-b")==0) {
	 Par->print_file9_ind=1;
	 Files->usingFile[8]=1;
       }
       if (strcmp(argv[i],"-m")==0)
	 {
	   Par->mutation_rate_ind=1;
	 }
       if(strcmp(argv[i],"-J")==0)
	 {
	   Par->jitter_locations=1;
	 }
      if (strcmp(argv[i],"-v")==0)
	{
	  Par->verbose=1;
	}
      if (strcmp(argv[i],"-vv")==0)
	{
	  Par->verbose=1;
	  Par->vverbose=1;
	}
      if (strcmp(argv[i],"-d")==0)
	{
	  Par->printnorecprobs=1;
	  Files->usingFile[9]=1;
	}
      if (strcmp(argv[i],"-e")==0)
	{ //// Use a limited data range for EM estimation
	  if (argv[(i+1)][0] == '-')
	    {
	      fprintf(Par->out,"Something wrong with input command line at \"-e\" (missing arguments?). Exiting....\n");
	      usage(Par->out);
	      stop_on_error(1,Par->errormode,Par->err);
	    }
	  Par->emloci = atoi(argv[(i+1)]);
	  if(Par->vverbose) fprintf(Par->out,"Read -e (using a block of %i loci)\n",Par->emloci);
	}
      if (strcmp(argv[i],"-S")==0)
	{ //// Use a limited data range for EM estimation
	  if (argv[(i+1)][0] == '-')
	    {
	      fprintf(Par->out,"Something wrong with input command line at \"-S\" (missing arguments?). Exiting....\n");
	      usage(Par->out);
	      stop_on_error(1,Par->errormode,Par->err);
	    }
	  Par->rseed = atoi(argv[(i+1)]);
	  if(Par->vverbose) fprintf(Par->out,"Read -S (using a seed of %i)\n",Par->rseed);
	}
      if (strcmp(argv[i],"-l")==0)
	{ //// Use a limited data range for EM estimation
	  if ((argv[(i+1)][0] == '-') || (argv[(i+2)][0] == '-'))
	    {
	      fprintf(Par->out,"Something wrong with input command line at \"-l\" (missing arguments?). Exiting....\n");
	      usage(Par->out);
	      stop_on_error(1,Par->errormode,Par->err);
	    }
	  Par->startlocus = atoi(argv[(i+1)])-1;
	  Par->endlocus = atoi(argv[(i+2)])-1;
	  if(Par->vverbose) fprintf(Par->out,"Read -l (process data starting at (0 indexed) locus %i and ending at locus %i containing %i SNPs)\n",Par->startlocus,Par->endlocus,Par->endlocus-Par->startlocus);
	}
    }
  if(Par->vverbose) fprintf(Par->out,"Reached end of first parameter pass\n");
  /// End first loop
  for (i=1; i < argc; i++)
    {
      if (strcmp(argv[i],"-g")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      fprintf(Par->out,"Something wrong with input command line (missing arguments?). Exiting....\n");
	      fprintf(Par->out,"In argument %i (%s)\n",i,argv[(i+1)]);
	      usage(Par->out);
	      stop_on_error(1,Par->errormode,Par->err);
	    }
	  setPhase(argv[(i+1)],Infiles);
	  if(Par->verbose) fprintf(Par->out,"Detected -g %s\n",Infiles->phase);

	}
      if (strcmp(argv[i],"-t")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      fprintf(Par->out,"Something wrong with input command line (missing arguments?). Exiting....\n");
	      fprintf(Par->out,"In argument %i (%s)\n",i,argv[(i+1)]);
	      usage(Par->out);
	      stop_on_error(1,Par->errormode,Par->err);
	    }
	  setIdfile(argv[(i+1)],Infiles);
	  if(Par->verbose) fprintf(Par->out,"Detected -t %s\n",Infiles->id);

	}
      if (strcmp(argv[i],"-Rg")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      fprintf(Par->out,"Something wrong with input command line (missing arguments?). Exiting....\n");
	      fprintf(Par->out,"In argument %i (%s)\n",i,argv[(i+1)]);
	      usage(Par->out);
	      stop_on_error(1,Par->errormode,Par->err);
	    }
	  setPhase2(argv[(i+1)],Infiles);
	  if(Par->verbose) fprintf(Par->out,"Detected -Rg %s\n",Infiles->phase2);

	}
      if (strcmp(argv[i],"-Rt")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      fprintf(Par->out,"Something wrong with input command line (missing arguments?). Exiting....\n");
	      fprintf(Par->out,"In argument %i (%s)\n",i,argv[(i+1)]);
	      usage(Par->out);
	      stop_on_error(1,Par->errormode,Par->err);
	    }
	  setIdfile2(argv[(i+1)],Infiles);
	  if(Par->verbose) fprintf(Par->out,"Detected -t %s\n",Infiles->id2);

	}
      if (strcmp(argv[i],"-r")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      fprintf(Par->out,"Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage(Par->out);
	      stop_on_error(1,Par->errormode,Par->err);
	    }
	  setRecmap(argv[(i+1)],Infiles);
	  if(Par->verbose) fprintf(Par->out,"Detected -r file %s\n",Infiles->recmap);
	}
      if (strcmp(argv[i],"-f")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      fprintf(Par->out,"Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage(Par->out);
	      stop_on_error(1,Par->errormode,Par->err);
	    }
	  setDonorlist(argv[(i+1)],Infiles);
	  Par->start_val = atoi(argv[(i+2)]);
	  Par->end_val = atoi(argv[(i+3)]);
	  if ((Par->end_val < Par->start_val) || (Par->start_val < 0) || (Par->end_val < 0))
	    {
	      fprintf(Par->out,"Invalid start_ind/stop_ind vals ('-f' switch). If you want to condition each recipient individual on each donor individual, use '-f 0 0'. Exiting...\n");
	      stop_on_error(1,Par->errormode,Par->err);
	    }
	  if (Par->start_val > 0) Par->start_val=Par->start_val-1;

	  if(Par->verbose) fprintf(Par->out,"Detected -f file %s\n",Infiles->donorlist);

	}
     if (strcmp(argv[i],"-i")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      fprintf(Par->out,"Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage(Par->out);
	      stop_on_error(1,Par->errormode,Par->err);
	    }
	  Par->EMruns = atoi(argv[(i+1)]);
	  if(Par->verbose) fprintf(Par->out,"Detected -i %i\n",Par->EMruns);

	  if (Par->EMruns < 0)
	    {
	      fprintf(Par->out,"Number of EM runs must be at least 0. Exiting...\n");
	      stop_on_error(1,Par->errormode,Par->err);
	    }
	  if ((Par->EMruns>0) && (Par->copy_prop_em_find==0) && (Par->recom_em_find==0) && (Par->mutation_em_find==0) && (Par->mutationALL_em_find==0))
	    {
	      fprintf(Par->out,"You have specified to perform E-M iterations, but have not specified which parameter(s) to maximize. If using '-i' switch, please specify at least one of '-in', '-ip', '-im', and/or '-iM'. Exiting...\n");
	      stop_on_error(1,Par->errormode,Par->err);
	    }
	}
      if (strcmp(argv[i],"-s")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      fprintf(Par->out,"Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage(Par->out);
	      stop_on_error(1,Par->errormode,Par->err);
	    }
	  Par->samplesTOT = atoi(argv[(i+1)]);
	  if(Par->verbose) fprintf(Par->out,"Detected -s %i\n",Par->samplesTOT);
	  if (Par->samplesTOT < 0)
	    {
	      fprintf(Par->out,"Number of samples must be >= 0. Exiting...\n");
	      stop_on_error(1,Par->errormode,Par->err);
	    }
	}
       if (strcmp(argv[i],"-n")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      fprintf(Par->out,"Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage(Par->out);
	      stop_on_error(1,Par->errormode,Par->err);
	    }
	  Par->readN_e = atof(argv[(i+1)]);
	  Par->N_e=Par->readN_e;
	  if (Par->readN_e <= 0)
	    {
	      fprintf(Par->out,"Recombination scaling parameter N_e must be > 0. Exiting...\n");
	      stop_on_error(1,Par->errormode,Par->err);
	    }
	}
       if (strcmp(argv[i],"-M")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      fprintf(Par->out,"Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage(Par->out);
	      stop_on_error(1,Par->errormode,Par->err);
	    }
	  Par->GlobalMutRate = atof(argv[(i+1)]);
	  if (Par->GlobalMutRate==0) Par->GlobalMutRate=-9;
	  if (Par->GlobalMutRate < 0)
	    fprintf(Par->out,"Mutation (emission) parameter must be > 0. Using Li & Stephens (2003) version of Watterson's estimate instead of user-supplied value...\n");
	}
       if (strcmp(argv[i],"-k")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      fprintf(Par->out,"Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage(Par->out);
	      stop_on_error(1,Par->errormode,Par->err);
	    }
	  Par->region_size = atof(argv[(i+1)]);
	  if (Par->region_size < 1)
	    {
	      fprintf(Par->out,"Region_size must be >= 1. Exiting...\n");
	      stop_on_error(1,Par->errormode,Par->err);
	    }
	}
       if (strcmp(argv[i],"-m")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      fprintf(Par->out,"Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage(Par->out);
	      stop_on_error(1,Par->errormode,Par->err);
	    }
	  Par->mut_rate_self = atof(argv[(i+1)]);
	}
        if (strcmp(argv[i],"-a")==0)
	{
	  if ((argv[(i+1)][0] == '-') || (argv[(i+2)][0] == '-'))
	    {
	      fprintf(Par->out,"Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage(Par->out);
	      stop_on_error(1,Par->errormode,Par->err);
	    }
	  Par->start_val = atoi(argv[(i+1)]);
	  Par->end_val = atoi(argv[(i+2)]);
	  if(Par->verbose) fprintf(Par->out,"Detected -a %i %i\n",Par->start_val,Par->end_val);
	  if ((Par->end_val < Par->start_val) || (Par->start_val < 0) || (Par->end_val < 0))
	    {
	      fprintf(Par->out,"Invalid start_ind/stop_ind vals ('-a' switch). If you want to condition each individual on every other individual, use '-a 0 0'. Exiting...\n");
	      stop_on_error(1,Par->errormode,Par->err);
	    }
	  if (Par->start_val > 0) Par->start_val=Par->start_val-1;
	}
       if (strcmp(argv[i],"-o")==0)
	 {
	   if (argv[(i+1)][0] == '-')
	     {
	       fprintf(Par->out,"Something wrong with input command line (missing arguments?). Exiting....\n");
	       usage(Par->out);
	       stop_on_error(1,Par->errormode,Par->err);
	     }
	   if(Par->verbose) fprintf(Par->out,"Detected -o %s\n",argv[(i+1)]);
	   specifyOutfiles(argv[i+1],Files);
	 }
    }
}

void printInformation(struct files_t *Outfiles,struct infiles_t *Infiles,struct param_t *Par,struct data_t *Data){
  if (Par->allvsall && !Par->geno2_find) fprintf(Par->out,"Will condition each individual on every other individual...\n");
  if (Par->allvsall && Par->geno2_find) fprintf(Par->out,"Will condition each individual in file %s on every individual in file %s...\n",Infiles->phase2,Infiles->phase);
  if (Par->donorlist_find) fprintf(Par->out,"Will use donors and recipients specified in %s...\n",Infiles->donorlist);
  if (Par->donorlist_find && Par->geno2_find) fprintf(Par->out,"Will use additional recipients specified in %s...\n",Infiles->phase2);
  if (Par->geno2_find && !Par->idfile2_find) fprintf(Par->out,"Will name additional recipients IND...\n");
  if (Par->geno2_find && Par->idfile2_find) fprintf(Par->out,"Will name additional recipients from file %s...\n",Infiles->id2);
  
  //  if (Data->nhaps_startpop==0 && Par->donorlist_find==1) fprintf(Par->out,"Will use %s for population labels in output (even though there are no donor haplotypes)\n",Infiles->donorlist);

  if (Par->haploid_ind==1) fprintf(Par->out,"Assuming all inds are haploid....\n");
  if (Par->unlinked_ind==1)
    {
      fprintf(Par->out,"Assuming sites are unlinked....\n");
      //Par->EMruns=0;
    }
  if (Par->prior_donor_probs_ind==1)
    fprintf(Par->out,"Using specified prior donor probs from input file....\n");
  if ((Par->mutation_rate_ind==1) && (Par->mutationALL_em_find==0))
    fprintf(Par->out,"Using specified mutation rates from input file....\n");
  if (Par->copy_prop_em_find==1)
    fprintf(Par->out,"Running E-M to estimate copying proportions....\n");
  if (Par->recom_em_find==1)
    fprintf(Par->out,"Running E-M to estimate N_e....\n");
  if (Par->mutation_em_find==1)
    fprintf(Par->out,"Running E-M to estimate mutation (emission) probabilities....\n");
  if (Par->mutationALL_em_find==1)
    fprintf(Par->out,"Running E-M to estimate global mutation (emission) probability....\n");
  if (Par->jitter_locations==1)
    fprintf(Par->out,"Moving SNP locations where collisions occur....\n");
  fprintf(Par->out," Number of EM-runs = %d\n Number of samples = %d\n N_e value = %lf\n Region size = %lf\n",Par->EMruns,Par->samplesTOT,Par->N_e,Par->region_size);
  fprintf(Par->out," Global mutation value = %lf\n",Par->GlobalMutRate);
  //  fprintf(Par->out," Number of potential donor haplotypes = %d\n Number of recipient haplotypes = %d\n",Data->nhaps,Data->nhapstotal-Data->nhaps_startpop);

  if(Outfiles->usingFile[0]) gzprintf(*Outfiles->fout1, "EM_iter = %d (N_e = %d / copy_prop = %d / mutation = %d / mutationGLOBAL = %d), nsamples = %d, N_e_start = %lf, region_size = %lf, genotype dataset = %s, genmap dataset = %s, donor-list dataset = %s\n", Par->EMruns, Par->recom_em_find, Par->copy_prop_em_find, Par->mutation_em_find, Par->mutationALL_em_find, Par->samplesTOT, Par->N_e, Par->region_size, Infiles->phase,Infiles->recmap,Infiles->donorlist);
}

void cleanup(int retval,struct param_t *Par,struct donor_t *Donors,
	     struct copyvec_t *Copyvec,struct files_t *Outfiles,
	     struct infiles_t *Infiles,struct data_t *Data,struct ids_t *Ids){
  if(Par->verbose) fprintf(Par->out,"Cleanup...\n");
  fflush(Par->out);

  if(Par==NULL) return;
  if(Par->verbose) fprintf(Par->out,"Freeing memory\n");
  fflush(Par->out);
  clearDonors(Donors,Par);
  if(Par->verbose) fprintf(Par->out,"Clearing Copyvec\n");
  fflush(Par->out);
  clearCopyvec(Copyvec);

  if(Par->verbose) fprintf(Par->out,"Closing files\n");
  fflush(Par->out);
  closeOutfiles(Outfiles);
  if(Par->vverbose) fprintf(Par->out,"Closing infiles\n");
  fflush(Par->out);
  closeInfiles(Infiles);

  if(Par->vverbose) fprintf(Par->out,"Freeing infiles\n");
  fflush(Par->out);
  freeInfiles(Infiles);
  if(Par->vverbose) fprintf(Par->out,"Freeing outfiles\n");
  fflush(Par->out);
  freeOutfiles(Outfiles);
  if(Par->vverbose) fprintf(Par->out,"Freeing data\n");
  fflush(Par->out);
  DestroyData(Data);

  if(Par->vverbose) fprintf(Par->out,"Freeing ids\n");
  fflush(Par->out);
  DestroyIds(Ids);

  if(retval==0){
    fprintf(Par->out,"%s\n",cpsuccesstext);
  }else{
    fprintf(Par->out,"Chromopainter ended on error\n");
  }
  fflush(Par->out);
  DestroyParam(Par);
}


int chromopainter(int argc, char *argv[])
{
  FILE *mainout=stdout;
  struct param_t *Par=DefaultParam();
  Par->out=mainout;
  struct donor_t *Donors=NULL;
  struct copyvec_t *Copyvec=NULL;
  struct infiles_t *Infiles = defaultInfiles();
  struct ids_t *Ids=NULL; 
  struct data_t *Data=NULL; 
  struct files_t *Outfiles=defaultOutfiles();

  if(setjmp(Par->err)){ // this is c error handling; its mad
    // google setjmp error handling for details
    // we only play with it if --noexitonerrors is used
    fprintf(Par->out,"Cleaning up after an error was thrown.\n");
    cleanup(1,Par,Donors,Copyvec,Outfiles,Infiles,Data,Ids);
    return(1); // return an error
  }
  unsigned int c1=0;
  while (c1<argc-1) {
    if (strcmp(argv[c1],"--noexitonerrors")==0)
      if(argv[(c1+1)][0] != '-') {
	mainout = fopen(argv[(c1+1)],"w");
	fflush(mainout);
      }
    c1++;
  }
  if(mainout!=stdout)  Par->out=mainout;
  if(argc<2) {
    usage(mainout);
    return(0);
  }
  
  c1=0;
  while (c1<argc){ fprintf(mainout,"%s\n",argv[c1]); c1++;};

  fprintf(Par->out,"Assigning parameters from command line\n");
  assignParameters(Par,Infiles,Outfiles,argc,argv);
  
  if(Par->vverbose) fprintf(Par->out,"Opening output files\n");
  openOutfiles(Outfiles);
  if(Par->vverbose) fprintf(Par->out,"Validating output files\n");
  if(!validateOutfiles(Outfiles,Par->out)){
    fprintf(Par->out,"Failed to open output files; terminating!\n");
    stop_on_error(1,Par->errormode,Par->err);
  }
  if(Par->verbose) fprintf(Par->out,"Created output files\n");


  int i,j;
  int log_lik_check;
  initRNG(Par);

  /***********************************************************/
  // DEFAULT VALUES:

  if(Par->verbose) fprintf(Par->out,"Checking parameter validity\n");
  parameterCheck(Par);
  // END CHECKs
  if(Par->verbose) fprintf(Par->out,"Parameters pass validity checks\n");

  ///////////////////////////////
  // Create the IDs
  if(Par->idfile_find==1){
    if(Par->verbose) fprintf(Par->out,"Opening ID file %s\n",Infiles->id);
    Ids=ReadIdfile(Infiles->id,Par);
  }else{ 
    if(Par->verbose) fprintf(Par->out,"Creating IDs from phase file %s\n",Infiles->phase);
    Ids=CreateIdsFromPhase(Infiles->phase,Par);
  }
  // Read in the additional ids, if present
  struct ids_t *Ids2;
  if(Par->geno2_find){
    struct ids_t *Ids1=Ids;
    if(Par->idfile2_find){
      if(Par->verbose) fprintf(Par->out,"Reading recipient IDs from file %s\n",Infiles->id2);
      Ids2=ReadIdfile(Infiles->id2,Par);
    }else{
      if(Par->verbose) fprintf(Par->out,"Creating recipient IDs from phase file %s\n",Infiles->phase2);
      Ids2=CreateIdsFromPhase(Infiles->phase2,Par);
    }
    if(Par->vverbose) fprintf(Par->out,"Merging IDs\n");
    Ids=CreateMergedIds(Ids1,Ids2,Par); 
    if(Par->vverbose) fprintf(Par->out,"Removing unneeded ID objects\n");
    DestroyIds(Ids1);
  }

  //
  if(Par->verbose) {
    fprintf(Par->out,"Using %d individuals from %d total in data\n",Ids->nind_inc,Ids->nind_tot);
  }
  ///////////////////////////////
  
  // Read the phase file
      // open first file (to get information on haplotype numbers)
  if(Par->verbose) fprintf(Par->out,"Opening phase file %s\n",Infiles->phase);

  if(!openPhase(Infiles,Par->out)) { fprintf(Par->out,"error opening phase file %s\n",Infiles->phase); stop_on_error(1,Par->errormode,Par->err);}

  if(Par->verbose) fprintf(Par->out,"Reading data %s\n",Infiles->phase);
  Data = ReadData(Infiles,Ids,Par);
  // also closes phase files

  // Some processing of the data
  if(Par->verbose) fprintf(Par->out,"Performing additional checks\n");
  // More error checking
  validateData(Infiles,Par,Data);


  ////////////////////////////////
  /// DONORS!!!!
  //////////////////////////////////////////////
  Donors=createDonors(Infiles,Ids,Ids2,Data,Par);
  if(Par->idfile2_find){ // no longer needed
    DestroyIds(Ids2);
  }
  
  setNe(Data->condhaps,Par);
  // Note: this Ne is possibly the wrong one to use, as unused haplotypes in the file change the result? Not a problem if we do EM...
  // More error checking
  printInformation(Outfiles, Infiles,Par,Data);

  /* if(Par->vverbose) { */
  /*   fprintf(Par->out,"Main: Created %i donor populations containing %i haplotype(s):\n",Donors->ndonorpops,Data->condhaps); */
  /*   for(i=0;i<Donors->ndonorpops;i++){ */
  /*     fprintf(Par->out,"  %i: %s, %i\n",i+1,Donors->donorlabels[i],Donors->ndonorhaps_vec[i]); */
  /*   } */
  /* } */

  fprintf(Par->out," num donor pops = %d\n",Donors->ndonorpops);
  fprintf(Par->out," Using individuals [%d - %d] as recipients\n",Par->start_val+1,Par->end_val);

  
  // define copy_prob, MutProb_vec, and pop_vec:
  if(Par->vverbose) fprintf(Par->out,"Initialising copy vectors\n");
  Copyvec = initializeCopyVecs(Donors,Data,Par);

  if(Par->vverbose) fprintf(Par->out,"Creating recombination map...\n");
  assignRecMap(Copyvec,Infiles,Data,Par);

  Par->EMruns = Par->EMruns + 1; // the final em run is treated as a likelihood calculation only

  ////////////////////////////////
  // Constructing output files
                 /* print-out headers for copy-props, chunk counts, lengths, and differences: */
  if(Par->verbose) fprintf(Par->out,"Constructing skeleton output files\n");
  makeHeaders(Outfiles,Donors,Par);
  makeSNPbasedHeaders(Outfiles,Data->positions,Data->nsnps);

  //  if (Par->end_val==0 || Par->end_val > Data->recipinds) Par->end_val=Data->recipinds;
  if(Par->verbose) fprintf(Par->out,"Beginning processing\n");

  // Destroying data
  //  if(Par->vverbose) fprintf(Par->out,"Destroying data (will reread later)\n");

  // Do EVERYTHING!!!!!
  log_lik_check = loglik(Copyvec, Donors, Data, Ids, Infiles, Outfiles, Par);
  //  log_lik_check = loglik(Copyvec,Donors,
  //Data->nhaps_startpop, &Data->nsnps, Data->nhapstotal, Par->N_e, Copyvec->recom_map, Copyvec->MutProb_vec, Copyvec->copy_prob, Copyvec->copy_probSTART, Copyvec->pop_vec, Infiles->phase, Outfiles, Par);
  
  if(Par->verbose) fprintf(Par->out,"Processing complete\n");
  // Check it worked
  if (log_lik_check != 1)
    {
      fprintf(Par->out,"Algorithm failed. Check input files and parameters. Exiting....\n");
      stop_on_error(1,Par->errormode,Par->err);
    }

  cleanup(0,Par,Donors,Copyvec,Outfiles,Infiles,Data,Ids);
  return 0;
}
