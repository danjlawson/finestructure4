#ifndef INFEXTRACT2_H
#define INFEXTRACT2_H
#include <vector>
#include "state.h"
#include "prior.h"
#include "data.h"
#include "rng.h"
#include "node.h"
#include "fsxml.h"

namespace fines
{

/**
    @brief Inference algorithm: extraction of highest posterior state
*/
class InfExtract2
{
public:
  InfExtract2(my_rng * rng,Data *d,FsXml *fs,Prior *prior,int v=0,Data *dlength=NULL,Data *dref=NULL,int datainference=INFDATA_COUNTS,int modeltype=MODELTYPE_FINESTRUCTURE);///<Constructor 
    ~InfExtract2();
    inline State * getState(){return(state);}///< Returns the state

protected:
    Data *data;
    State *state;
    int verbose;
    my_rng * rng;
};

} // end namespace fines
#endif
