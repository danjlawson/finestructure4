#ifdef __cplusplus
extern "C" {
#endif

  // Start the header
#ifndef CHROMOPAINTERDATA_H
#define CHROMOPAINTERDATA_H

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <time.h>
#include <string.h>

#include "ChromoPainterReading.h"
#include "ChromoPainterInfiles.h"
#include "ChromoPainterPar.h"
  
  struct ids_t {
    int nind_tot; // total number of individuals present
    int nind_inc; // total number of individuals USED
    int *include_ind_vec;// which individuals are included in the analysis
    char ** id;// ID of each individual
    char ** popid; // population identifier for each individual
    int nind_file_main;// total number of individuals present in the main data file
  };// All information from the id file

  struct data_t {
    
    int nhaps_startpop; // number of haps specified in the phase file
    double nind; // number of inds specified in the phase file
    int nsnps_raw; // number of snps specified in the phase file

    int nsnps; // number of snps to be used
    int nhaps; // number of haplotypes USED = condhaps+hapsperind (in a-v-a mode)
    int condhaps; // number of haplotypes being conditioned on
    int hapsperind; // the ploidy of the data = number of chromosomes being painted

    int nhapstotal; // number of haplotypes IN TOTAL, ie number of lines in the phase file
    int reciphaps; // number of haplotypes being treated as recipient (was cond_nhaps)
    int recipinds; // number of individuals being treated as recipient (was cond_nhaps). NOT the size of cond_chromosomes, which are read individually!
    
    double *positions; // positions of the SNPs
    double *lambda; // recombination map between each SNP

    int *condnums; // the line number of the conditioned haplotypes in the sorting of the phase file (starting at 0 for the first haplotype)
    int *recipnums; // the line number of the painted haplotypes in the sorting of the phase file (starting at 0 for the first haplotype)

    int **all_chromosomes; // all the chromosomes

    // Things for the current painting
    int currentind; /// The current individual being painted, for which cond_chromosomes and ind_chromosomes are configured

    int current_donor_nind; // number of donor individuals
    int* current_donor_inds; // indices of donor individuals
    int* current_donor_haps; // indices of donor haplotypes
    int **cond_chromosomes; // just those chromosomes being conditioned on NOW, i.e. length of the ploidy of the data
    int **ind_chromosomes; // just those chromosomes being painted
  };


  int codeHaplotypes(char in);///< code SNP values into integers

  struct ids_t *ReadIdfile(char *filename, struct param_t *Par); // Read ID file
  struct ids_t *CreateIdsFromPhase(char *filename, struct param_t *Par); // Create ID from phase file
  struct ids_t *CreateMergedIds(struct ids_t *Ids1, struct ids_t *Ids2, struct param_t *Par); // Merge the ID files (1 then 2), copied into a further id file
  
  struct data_t *ReadDataHeader(FILE *fd, struct param_t *Par); // Read the header of a phase file
  struct data_t *ReadData(struct infiles_t *Infiles, struct ids_t *Ids,struct param_t *Par); // Read chromopainter style phase format
  void DestroyData(struct data_t *Data); // destroy the data
  void DestroyIds(struct ids_t *Ids); // destroy the ids
  void setIndAsRecipient(int m, struct data_t *Data,struct ids_t *Ids,struct param_t *Par); // set individual m to be the recipient

  int * getallelic_type_count_vec(struct data_t * Data); // count the number of each allele

#endif
  // End the header

#ifdef __cplusplus
}
#endif
