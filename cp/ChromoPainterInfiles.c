
#include "ChromoPainterInfiles.h"

struct infiles_t *defaultInfiles()  {
  struct infiles_t *Files = malloc(sizeof(struct infiles_t));
  int i;
  // Create and define memory for files
  Files->recmap = malloc(1000 * sizeof(char));
  Files->phase = malloc(1000 * sizeof(char));
  Files->donorlist = malloc(1000 * sizeof(char));
  Files->id = malloc(1000 * sizeof(char));
  Files->phase2 = malloc(1000 * sizeof(char));
  Files->id2 = malloc(1000 * sizeof(char));
  for(i=0;i<6;i++){
    Files->usingFile[i]=0;
    Files->fopen[i]=0;
  }
  strcpy(Files->recmap,"NULL");
  strcpy(Files->donorlist,"NULL");
  strcpy(Files->phase,"NULL");
  strcpy(Files->id,"NULL"); 
  strcpy(Files->phase2,"NULL");
  strcpy(Files->id2,"NULL");
 
  Files->frecmap=NULL;
  Files->fphase=NULL;
  Files->fdonorlist=NULL;
  Files->fid=NULL;
  Files->fphase2=NULL;
  Files->fid2=NULL;
  
  return(Files);
}

void setPhase(char *filenameGEN, struct infiles_t *Files){
  strcpy(Files->phase ,filenameGEN);
  Files->usingFile[0]=1;
}

void setRecmap(char *filename, struct infiles_t *Files){
  strcpy(Files->recmap ,filename);
  Files->usingFile[1]=1;
}

void setDonorlist(char *filenameDONORLIST, struct infiles_t *Files){
  strcpy(Files->donorlist ,filenameDONORLIST);
  Files->usingFile[2]=1;
}

void setIdfile(char *filenameID, struct infiles_t *Files){
  strcpy(Files->id ,filenameID);
  Files->usingFile[3]=1;
}

void setPhase2(char *filenameGEN, struct infiles_t *Files){
  strcpy(Files->phase2 ,filenameGEN);
  Files->usingFile[4]=1;
}

void setIdfile2(char *filenameID, struct infiles_t *Files){
  strcpy(Files->id2 ,filenameID);
  Files->usingFile[5]=1;
}

int openPhase(struct infiles_t *Files,FILE *mainout){
  if(Files->usingFile[0] && !Files->fopen[0]) {
    Files->fphase = fopen(Files->phase, "r");  
    Files->fopen[0] =1;
  }else return(0);
  return(validateInfiles(Files,mainout));
}

int openRecmap(struct infiles_t *Files, FILE *mainout){
  if(Files->usingFile[1] && !Files->fopen[1]) {
    Files->frecmap = fopen(Files->recmap, "r");  
    Files->fopen[1] =1;
  }else return(0);
  return(validateInfiles(Files,mainout));
}

int openDonorlist(struct infiles_t *Files, FILE *mainout){
  if(Files->usingFile[2] && !Files->fopen[2]) {
    Files->fdonorlist = fopen(Files->donorlist, "r");  
    Files->fopen[2] =1;
  }else return(0);
  return(validateInfiles(Files,mainout));
}

int openId(struct infiles_t *Files, FILE *mainout){
  if(Files->usingFile[3] && !Files->fopen[3]) {
    Files->fid = fopen(Files->id, "r");  
    Files->fopen[3] =1;
  }else return(0);
  return(validateInfiles(Files,mainout));
}

int openPhase2(struct infiles_t *Files,FILE *mainout){
  if(Files->usingFile[4] && !Files->fopen[4]) {
    Files->fphase2 = fopen(Files->phase2, "r");  
    Files->fopen[4] =1;
  }else return(0);
  return(validateInfiles(Files,mainout));
}
int openId2(struct infiles_t *Files, FILE *mainout){
  if(Files->usingFile[5] && !Files->fopen[5]) {
    Files->fid2 = fopen(Files->id2, "r");  
    Files->fopen[5] =1;
  }else return(0);
  return(validateInfiles(Files,mainout));
}

void closePhase(struct infiles_t *Files){
  // Close the output files
  if(Files->fopen[0]){fclose(Files->fphase);}
  Files->fopen[0]=0;
}

void closeRecmap(struct infiles_t *Files){
  // Close the output files
  if(Files->fopen[1]){fclose(Files->frecmap);}
  Files->fopen[1]=0;
}

void closeDonorlist(struct infiles_t *Files){
  // Close the output files
  if(Files->fopen[2]){fclose(Files->fdonorlist);}
  Files->fopen[2]=0;
}

void closeId(struct infiles_t *Files){
  // Close the output files
  if(Files->fopen[3]){fclose(Files->fid);}
  Files->fopen[3]=0;
}

void closePhase2(struct infiles_t *Files){
  // Close the output files
  if(Files->fopen[4]){fclose(Files->fphase2);}
  Files->fopen[4]=0;
}

void closeId2(struct infiles_t *Files){
  // Close the output files
  if(Files->fopen[5]){fclose(Files->fid2);}
  Files->fopen[5]=0;
}

void closeInfiles(struct infiles_t *Files){
  if(Files==NULL) return;
  // Close the input files
  closePhase(Files);
}

int validateInfiles(struct infiles_t *Files, FILE *mainout){
  // check that all files that think they are open, are in fact open
  if (Files->fphase == NULL && Files->fopen[0]) {fprintf(mainout,"error opening phase file\n"); return(0);}
  if (Files->frecmap == NULL && Files->fopen[1]) {fprintf(mainout,"error opening recmap\n"); return(0);}
  if (Files->fdonorlist == NULL && Files->fopen[2]) {fprintf(mainout,"error opening donor file\n"); return(0);}
  if (Files->fid == NULL && Files->fopen[3]) {fprintf(mainout,"error opening ID file\n"); return(0);}
  if (Files->fphase2 == NULL && Files->fopen[4]) {fprintf(mainout,"error opening recipient phase file\n"); return(0);}
  if (Files->fid2 == NULL && Files->fopen[5]) {fprintf(mainout,"error opening recipient ID file\n"); return(0);}
  return(1);
}

void freeInfiles(struct infiles_t *Files) {
  if(Files==NULL) return;
  closeInfiles(Files);
  free(Files->recmap);
  free(Files->phase);
  free(Files->donorlist);
  free(Files->id);
  free(Files->phase2);
  free(Files->id2);
  free(Files);
}
