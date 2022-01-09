#ifdef __cplusplus
extern "C" {
#endif

  // Start the header
#ifndef CHROMOPAINTERINFILES_H
#define CHROMOPAINTERINFILES_H

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <time.h>
#include <string.h>

#include "zlib.h"
#include "ChromoPainterError.h"

  struct infiles_t {
    char *recmap; // was filename
    char *phase; // was filenameGEN
    char *donorlist; // was filenameDONORLIST
    char *id; // id file
    char *phase2; // name of the recipient phase file
    char *id2; // name of the recipient id file
    int usingFile[6];// 0 if we are not using each inputfile (phase,rec,donor,id,phase2,id2)
    int fopen[6];// 0 if inputfile is open
    FILE *frecmap; // rec file 
    FILE *fphase; // Phase file , was fd
    FILE *fdonorlist; // Donor file , was fd3
    FILE *fid; // Id file, was fd4
    FILE *fphase2; // Phase file for recipients
    FILE *fid2; // Id file for recipients
  };

  struct infiles_t *defaultInfiles();

  void setRecmap(char *filename, struct infiles_t *Files);
  void setPhase(char *filenameGEN, struct infiles_t *Files);
  void setDonorlist(char *filenameDONORLIST, struct infiles_t *Files);
  void setIdfile(char *filenameID, struct infiles_t *Files);
  void setPhase2(char *filenameGEN, struct infiles_t *Files);
  void setIdfile2(char *filenameID, struct infiles_t *Files);

  int openPhase(struct infiles_t *Files, FILE *mainout); // open the phase file (if it is supposed to be used)
  int openRecmap(struct infiles_t *Files, FILE *mainout); // open the recmap file (if it is supposed to be used)
  int openDonorlist(struct infiles_t *Files, FILE *mainout); // open the donorlist file (if it is supposed to be used)
  int openId(struct infiles_t *Files, FILE *mainout); // open the ID file (if it is supposed to be used)
  int openPhase2(struct infiles_t *Files, FILE *mainout); // open the phase file (if it is supposed to be used)
  int openId2(struct infiles_t *Files, FILE *mainout); // open the ID file (if it is supposed to be used)

  void closePhase(struct infiles_t *Files); // close the phase file
  void closeRecmap(struct infiles_t *Files); // close the phase file
  void closeDonorlist(struct infiles_t *Files); // close the phase file
  void closeId(struct infiles_t *Files); // close the id file
  void closePhase2(struct infiles_t *Files); // close the phase file
  void closeId2(struct infiles_t *Files); // close the id file

  void closeInfiles(struct infiles_t *Files); // close all input files
  int validateInfiles(struct infiles_t *Files, FILE *mainout); // check files are open
  void freeInfiles(struct infiles_t *Files); // free (and close if necessary) infiles

#endif
  // End the header

#ifdef __cplusplus
}
#endif
