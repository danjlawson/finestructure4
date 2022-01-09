#ifdef __cplusplus
extern "C" {
#endif

#ifndef CHROMOPAINTEREM_H
#define CHROMOPAINTEREM_H


  
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <time.h>
#include <string.h>

#include "zlib.h"
#include "ChromoPainterPar.h"
#include "ChromoPainterInfiles.h"
#include "ChromoPainterData.h"
#include "ChromoPainterRecmap.h"
#include "ChromoPainterSampler.h"
#include "ChromoPainterOutfiles.h"
#include "ChromoPainterError.h"

static const char * const cpsuccesstext="ChromoPainter completed successfully.";

  void usage(FILE *mainout);
  int chromopainter(int argc, char *argv[]);
  void validateData(struct infiles_t *Infiles,struct param_t *Par,struct data_t *Data);
  void printInformation(struct files_t *Outfiles,struct infiles_t *Infiles,struct param_t *Par,struct data_t *Data);

  void assignParameters(struct param_t *Par,struct infiles_t *Infiles,struct files_t *Files, int argc, char *argv[]);// Assign parameters based on the command line arguments

  void initRNG(struct param_t *Par);
  int randInt(int min,int max); ///returns a random integer in the range [min,max), i.e. inclusive of min and exclusive of max


#endif

#ifdef __cplusplus
}
#endif 
