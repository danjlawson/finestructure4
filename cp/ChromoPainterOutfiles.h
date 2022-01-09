#ifdef __cplusplus
extern "C" {
#endif

  // Start the header
#ifndef CHROMOPAINTEROUTFILES_H
#define CHROMOPAINTEROUTFILES_H

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <time.h>
#include <string.h>

#include "zlib.h"
#include "ChromoPainterError.h"
#include "ChromoPainterPar.h"
#include "ChromoPainterDonors.h"

  struct files_t {///< Hold all of the information required to create output files
    char *filenameOUT; // root of the output file names
    char **filenames; // actual output file names
    FILE *fout2, *fout3, *fout4, *fout5, *fout6, *fout7, *fout8; // each output file
    gzFile tfout1,tfout9,tfout10;// the gzipped samples, copyprobsfile and notransitionprob (Stored like this for convenience?)
    gzFile *fout1,*fout9,*fout10;// the gzipped samples, copyprobsfile and notransitionprob
    int usingFile[10];// 0 if we are not using each outputfile
    struct donor_t *Donors;///< The donor information
  };

  struct files_t *defaultOutfiles();
  void emoutfiles(struct files_t *Outfiles); // just use EM outfile
  void specifyOutfiles(char *filenameOUT,struct files_t *Outfiles);
  void openOutfiles(struct files_t *Outfiles); // open all output files
  void closeOutfiles(struct files_t *Outfiles); // close all output files
  int validateOutfiles(struct files_t *Outfiles, FILE *mainout); // check files are open

  void makeHeaders(struct files_t *Outfiles, struct donor_t *Donors,struct param_t *Par);// make the headers for the normal files
  void makeSNPbasedHeaders(struct files_t *Outfiles,double * snp_locations,int nloci);

  void printTransitionProbHeader(double * snp_locations,int nloci,struct files_t *Outfiles);
  void printSamplesHeader(double * snp_locations,int nloci,struct files_t *Outfiles);
  void printTransitionProb(double * exp_trans_prob,int ind_val,int nloci,struct files_t *Outfiles);
  void printSummary(int m,int num_regions_tot,double **copy_prob_pop,double *total_counts,double *total_lengths,double *total_differences,double *total_region_counts,double *total_squared_region_counts,struct files_t *Outfiles, struct param_t *Par);

  void printCopyProbs(double * exp_copy_pop,int ind_val, double tpos,struct files_t *Outfiles,struct param_t *Par);

  void freeOutfiles(struct files_t *Outfiles);

#endif
  // End the header

#ifdef __cplusplus
}
#endif
