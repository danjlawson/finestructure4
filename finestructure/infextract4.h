#ifndef INFEXTRACT4_H
#define INFEXTRACT4_H
#include <vector>
#include "state.h"
#include "data.h"
#include "prior.h"
#include "rng.h"
#include "node.h"
#include "fsxml.h"

namespace fines
{

/**
    @brief Inference algorithm: extraction of posterior distribution
*/
class InfExtract4
{
public:
  InfExtract4(my_rng * rng,Data *d,FsXml *fs,Prior *prior,int v=0);///<Constructor 
    ~InfExtract4();
    inline std::vector<double> getPosteriors(){return(posteriors);};///< gets the posterior probabilities
protected:
    Data *data;
    State *state;
    int verbose;
    std::vector<double> posteriors;
    long counts;
    my_rng * rng;
};

} // end namespace fines
#endif
