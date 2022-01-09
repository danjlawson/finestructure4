
#include "ChromoPainterOutfiles.h"

struct files_t *defaultOutfiles()  {
  struct files_t *Outfiles;
  char * tfilename = malloc(1000 * sizeof(char));
  int i;
  // Create and define memory for files
  Outfiles=malloc(sizeof(struct files_t));
  Outfiles->filenameOUT = malloc(1000 * sizeof(char));
  Outfiles->filenames= malloc(10 * sizeof(char*));
  for(i=0;i<10;i++){
    Outfiles->filenames[i]=malloc(1000 * sizeof(char));
  }
  
  for(i=0;i<8;i++){
    Outfiles->usingFile[i]=1;
  }
  Outfiles->usingFile[8]=0; // by default, don't use zipped copyprobsperlocus
  Outfiles->usingFile[9]=0; // by default, don't use zipped notransitionprobs

  Outfiles->Donors=NULL;
  Outfiles->fout1=NULL;
  Outfiles->fout2=NULL;
  Outfiles->fout3=NULL;
  Outfiles->fout4=NULL;
  Outfiles->fout5=NULL;
  Outfiles->fout6=NULL;
  Outfiles->fout7=NULL;
  Outfiles->fout8=NULL;
  Outfiles->fout9=NULL;
  Outfiles->fout10=NULL;
  return(Outfiles);
}

void emoutfiles(struct files_t *Outfiles) {
  int i;
  for(i=0;i<10;i++){
    Outfiles->usingFile[i]=0;
  }
  Outfiles->usingFile[2]=1;
}

void specifyOutfiles(char *filenameOUT,struct files_t *Outfiles)  {

  char * tfilename = malloc(1000 * sizeof(char));
  int i;
  // Create and define memory for files

  strcpy(Outfiles->filenameOUT ,filenameOUT);
  strcpy(tfilename,filenameOUT);

  strcat(tfilename,".samples.out.gz");
  strcpy(Outfiles->filenames[0],tfilename);
  strcpy(tfilename,filenameOUT);
  strcpy(Outfiles->filenames[1],strcat(tfilename,".prop.out"));
  strcpy(tfilename,filenameOUT);
  strcpy(Outfiles->filenames[2],strcat(tfilename,".EMprobs.out"));
  strcpy(tfilename,filenameOUT);
  strcpy(Outfiles->filenames[3],strcat(tfilename,".chunkcounts.out"));
  strcpy(tfilename,filenameOUT);
  strcpy(Outfiles->filenames[4],strcat(tfilename,".chunklengths.out"));
  strcpy(tfilename,filenameOUT);
  strcpy(Outfiles->filenames[5],strcat(tfilename,".mutationprobs.out"));
  strcpy(tfilename,filenameOUT);
  strcpy(Outfiles->filenames[6],strcat(tfilename,".regionchunkcounts.out"));
  strcpy(tfilename,filenameOUT);
  strcpy(Outfiles->filenames[7],strcat(tfilename,".regionsquaredchunkcounts.out"));
  strcpy(tfilename,filenameOUT);
  strcpy(Outfiles->filenames[8],strcat(tfilename,".copyprobsperlocus.out.gz"));
  strcpy(tfilename,filenameOUT);
  strcpy(Outfiles->filenames[9],strcat(tfilename,".transitionprobs.out.gz"));
  free(tfilename);
}

void openOutfiles(struct files_t *Outfiles){
  // Open all the files that are requested
  if(Outfiles->usingFile[0]){
    Outfiles->tfout1 = gzopen(Outfiles->filenames[0], "w");
    Outfiles->fout1 = &Outfiles->tfout1;
  }
  if(Outfiles->usingFile[1]){Outfiles->fout2 = fopen(Outfiles->filenames[1], "w");  }
  if(Outfiles->usingFile[2]){Outfiles->fout3 = fopen(Outfiles->filenames[2], "w");  }
  if(Outfiles->usingFile[3]){Outfiles->fout4 = fopen(Outfiles->filenames[3], "w");  }
  if(Outfiles->usingFile[4]){Outfiles->fout5 = fopen(Outfiles->filenames[4], "w");  }
  if(Outfiles->usingFile[5]){Outfiles->fout6 = fopen(Outfiles->filenames[5], "w");  }
  if(Outfiles->usingFile[6]){Outfiles->fout7 = fopen(Outfiles->filenames[6], "w");  }
  if(Outfiles->usingFile[7]){Outfiles->fout8 = fopen(Outfiles->filenames[7], "w");  }
  if(Outfiles->usingFile[8]){
    Outfiles->tfout9 = gzopen(Outfiles->filenames[8], "w");
    Outfiles->fout9 = &Outfiles->tfout9;	      
  }
  if(Outfiles->usingFile[9]){
    Outfiles->tfout10 = gzopen(Outfiles->filenames[9], "w");
    Outfiles->fout10 = &Outfiles->tfout10;
  }
}

void closeOutfiles(struct files_t *Outfiles){
  if(Outfiles==NULL) return;
  // Close the output files
  if(Outfiles->usingFile[0]){gzclose(*Outfiles->fout1);}
  if(Outfiles->usingFile[1]){fclose(Outfiles->fout2);}
  if(Outfiles->usingFile[2]){fclose(Outfiles->fout3);}
  if(Outfiles->usingFile[3]){fclose(Outfiles->fout4);}
  if(Outfiles->usingFile[4]){fclose(Outfiles->fout5);}
  if(Outfiles->usingFile[5]){fclose(Outfiles->fout6);}
  if(Outfiles->usingFile[6]){fclose(Outfiles->fout7);}
  if(Outfiles->usingFile[7]){fclose(Outfiles->fout8);}
  if(Outfiles->usingFile[8]){gzclose(*Outfiles->fout9);}
  if(Outfiles->usingFile[9]){gzclose(*Outfiles->fout10);}
}

int validateOutfiles(struct files_t *Outfiles,FILE *mainout){
  if (Outfiles->fout1 == NULL && Outfiles->usingFile[0]) {fprintf(mainout,"error opening file1\n"); return(0);}
  if (Outfiles->fout2 == NULL && Outfiles->usingFile[1]) {fprintf(mainout,"error opening file2\n"); return(0);}
  if (Outfiles->fout3 == NULL && Outfiles->usingFile[2]) {fprintf(mainout,"error opening file3\n"); return(0);}
  if (Outfiles->fout4 == NULL && Outfiles->usingFile[3]) {fprintf(mainout,"error opening file4\n"); return(0);}
  if (Outfiles->fout5 == NULL && Outfiles->usingFile[4]) {fprintf(mainout,"error opening file5\n"); return(0);}
  if (Outfiles->fout6 == NULL && Outfiles->usingFile[5]) {fprintf(mainout,"error opening file6\n"); return(0);}
  if (Outfiles->fout7 == NULL && Outfiles->usingFile[6]) {fprintf(mainout,"error opening file7\n"); return(0);}
  if (Outfiles->fout8 == NULL && Outfiles->usingFile[7]) {fprintf(mainout,"error opening file8\n"); return(0);}
  if (Outfiles->fout9 == NULL && Outfiles->usingFile[8]) {fprintf(mainout,"error opening file9\n"); return(0); }
  if (Outfiles->fout10 == NULL && Outfiles->usingFile[9]) {fprintf(mainout,"error opening file10\n"); return(0); }
  return(1);
}

void makeHeaders(struct files_t *Outfiles, struct donor_t *Donors,struct param_t *Par) {
  int i;
  int nlabels=Donors->ndonorpops;
  Outfiles->Donors=Donors;
  if(Outfiles->usingFile[1]) fprintf(Outfiles->fout2,"Recipient");
  if(Outfiles->usingFile[3]) fprintf(Outfiles->fout4,"Recipient");
  if(Outfiles->usingFile[4]) fprintf(Outfiles->fout5,"Recipient");
  if(Outfiles->usingFile[5]) fprintf(Outfiles->fout6,"Recipient");
  if(Outfiles->usingFile[6]) fprintf(Outfiles->fout7,"Recipient num.regions");
  if(Outfiles->usingFile[7]) fprintf(Outfiles->fout8,"Recipient num.regions");
  for (i=0; i < nlabels; i++)
    {
      if(Outfiles->usingFile[1]) fprintf(Outfiles->fout2," %s",Donors->donorlabels[i]);
      if(Outfiles->usingFile[3]) fprintf(Outfiles->fout4," %s",Donors->donorlabels[i]);
      if(Outfiles->usingFile[4]) fprintf(Outfiles->fout5," %s",Donors->donorlabels[i]);
      if(Outfiles->usingFile[5]) fprintf(Outfiles->fout6," %s",Donors->donorlabels[i]);
      if(Outfiles->usingFile[6]) fprintf(Outfiles->fout7," %s",Donors->donorlabels[i]);
      if(Outfiles->usingFile[7]) fprintf(Outfiles->fout8," %s",Donors->donorlabels[i]);
    }

  if(Outfiles->usingFile[1]) fprintf(Outfiles->fout2,"\n");
  if(Outfiles->usingFile[3]) fprintf(Outfiles->fout4,"\n");
  if(Outfiles->usingFile[4]) fprintf(Outfiles->fout5,"\n");
  if(Outfiles->usingFile[5]) fprintf(Outfiles->fout6,"\n");
  if(Outfiles->usingFile[6]) fprintf(Outfiles->fout7,"\n");
  if(Outfiles->usingFile[7]) fprintf(Outfiles->fout8,"\n");

  // Do the copyprobsperlocus file
  if (Outfiles->usingFile[8]) {
    gzprintf(* Outfiles->fout9,"pos");
    for (i=0; i < nlabels; i++) {
      gzprintf(*Outfiles->fout9," %s",Donors->donorlabels[i]);
    }
    gzprintf(*Outfiles->fout9,"\n");
  }
}

void printCopyProbs(double * exp_copy_pop,int ind_val, double tpos, struct files_t *Outfiles,struct param_t *Par){
  if(!Outfiles->usingFile[8]) return;
  struct donor_t *Donors=Outfiles->Donors;

  int i;
  gzprintf(*Outfiles->fout9,"%.0lf",tpos);//pos[(*p_Nloci-1)]);
 {
   for (i=0; i < Donors->ndonorpops; i++)
     {
       //       if (i == ind_val && Par->allvsall&&!Par->geno2_find) gzprintf(*Outfiles->fout9," 0.0");
       gzprintf(*Outfiles->fout9," %lf",exp_copy_pop[i]);
     }
  }
  gzprintf(*Outfiles->fout9,"\n");
}

void makeSNPbasedHeaders(struct files_t *Outfiles,double * snp_locations,int nloci){
  printTransitionProbHeader(snp_locations,nloci,Outfiles);
  printSamplesHeader(snp_locations,nloci,Outfiles);
}

void printSamplesHeader(double * snp_locations,int nloci,struct files_t *Outfiles){
  if(!Outfiles->usingFile[0]) return;
  // WE NO LONGER WRITE THE SNP HEADER, FOR COMPATABILITY WITH cpv2
  /*
  int locus;
  gzprintf(*Outfiles->fout1,"%.0lf",snp_locations[0]);
  for(locus=1;locus<nloci;locus++){
    gzprintf(*Outfiles->fout1," %.0lf",snp_locations[locus]); 
  }
  gzprintf(*Outfiles->fout1,"\n");
  */

}

void printTransitionProbHeader(double * snp_locations,int nloci,struct files_t *Outfiles){
  if(!Outfiles->usingFile[9]) return;

  int locus;
  gzprintf(*Outfiles->fout10,"HAP POS_%.0lf",snp_locations[0]);
  for(locus=1;locus<nloci;locus++){
    gzprintf(*Outfiles->fout10," POS_%.0lf",snp_locations[locus]); 
  }
  gzprintf(*Outfiles->fout10,"\n");
}

void printTransitionProb(double * exp_trans_prob,int ind_val,int nloci,struct files_t *Outfiles){
  if(!Outfiles->usingFile[9]) return;

  int locus;
  gzprintf(*Outfiles->fout10,"%i 0.0",ind_val);
  for(locus=1;locus<nloci;locus++){
    gzprintf(*Outfiles->fout10," %f",exp_trans_prob[locus-1]);
  }
  gzprintf(*Outfiles->fout10,"\n"); 
}

void printSummary(int m,int num_regions_tot,double **copy_prob_pop,double *total_counts,double *total_lengths,double *total_differences,double *total_region_counts,double *total_squared_region_counts,struct files_t *Outfiles, struct param_t *Par) {
  struct donor_t *Donors=Outfiles->Donors;
  int j=0;

  if(Outfiles->usingFile[1])  fprintf(Outfiles->fout2,"%s",Donors->reciplabels[m]);
  if(Outfiles->usingFile[3])  fprintf(Outfiles->fout4,"%s",Donors->reciplabels[m]);
  if(Outfiles->usingFile[4])  fprintf(Outfiles->fout5,"%s",Donors->reciplabels[m]);
  if(Outfiles->usingFile[5])  fprintf(Outfiles->fout6,"%s",Donors->reciplabels[m]);
  if(Outfiles->usingFile[6])  fprintf(Outfiles->fout7,"%s",Donors->reciplabels[m]);
  if(Outfiles->usingFile[7])  fprintf(Outfiles->fout8,"%s",Donors->reciplabels[m]);
 
  if(Outfiles->usingFile[6])  fprintf(Outfiles->fout7," %d",num_regions_tot);
  if(Outfiles->usingFile[7])  fprintf(Outfiles->fout8," %d",num_regions_tot);

  int ndonoreff=Donors->ndonorpops;
  for (j=0; j < ndonoreff; j++) {
    int jj=j;
      if(Outfiles->usingFile[1])  fprintf(Outfiles->fout2," %lf",copy_prob_pop[0][jj]);
      if(Outfiles->usingFile[3])  fprintf(Outfiles->fout4," %lf",total_counts[jj]);
      if(Outfiles->usingFile[4])  fprintf(Outfiles->fout5," %lf",total_lengths[jj]);
      if(Outfiles->usingFile[5])  fprintf(Outfiles->fout6," %lf",total_differences[jj]);
      if(Outfiles->usingFile[6])  fprintf(Outfiles->fout7," %lf",total_region_counts[jj]);
      if(Outfiles->usingFile[7])  fprintf(Outfiles->fout8," %lf",total_squared_region_counts[jj]);
  }

  if(Outfiles->usingFile[1])  fprintf(Outfiles->fout2,"\n");
  if(Outfiles->usingFile[3])  fprintf(Outfiles->fout4,"\n");
  if(Outfiles->usingFile[4])  fprintf(Outfiles->fout5,"\n");
  if(Outfiles->usingFile[5])  fprintf(Outfiles->fout6,"\n");
  if(Outfiles->usingFile[6])  fprintf(Outfiles->fout7,"\n");
  if(Outfiles->usingFile[7])  fprintf(Outfiles->fout8,"\n");

}


void freeOutfiles(struct files_t *Outfiles) {
  if(Outfiles==NULL) return;
  int i;
  for(i=0;i<10;i++){
    free(Outfiles->filenames[i]);
  }
  free(Outfiles->filenameOUT);
  free(Outfiles->filenames);
  free(Outfiles);
  
}

