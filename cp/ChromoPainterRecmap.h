#ifdef __cplusplus
extern "C" {
#endif

  // Start the header
#ifndef CHROMOPAINTERRECMAP_H
#define CHROMOPAINTERRECMAP_H

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <time.h>
#include <string.h>

#include "ChromoPainterDonors.h"
#include "ChromoPainterReading.h"
#include "ChromoPainterData.h"
#include "ChromoPainterInfiles.h"
#include "ChromoPainterPar.h"

  struct copyvec_t {
    int recom_map_size; // number of SNPs used in the recom map
    int ndonors; // number of donors

    // quantities given per  donor population:
    int * pop_vec;
    double * MutProb_vec;
    double * copy_prob;
    double * copy_probSTART;

    double * recom_map; // recom map
  };

  struct copyvec_t * initializeCopyVecs(struct donor_t *Donors,struct data_t *Data, struct param_t * Par) ;
  void assignRecMap(struct copyvec_t *Copyvec, struct infiles_t *Files,struct data_t *Data, struct param_t *Par);
  void clearCopyvec(struct copyvec_t *Copyvec);

#endif
  // End the header

#ifdef __cplusplus
}
#endif
