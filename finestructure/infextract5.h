#ifndef INFEXTRACT5_H
#define INFEXTRACT5_H
#include <vector>
#include "state.h"
#include "data.h"
#include "prior.h"
#include "rng.h"
#include "node.h"
#include "fsxml.h"
#include "infadmixture.h"

namespace fines
{

/**
    @brief Inference algorithm: extraction of likelihood distribution
*/
class InfExtract5
{
public:
  InfExtract5(my_rng * rng,Data *d,FsXml *fs,Prior *prior,long reps=1,int v=0);///<Constructor 
    ~InfExtract5();
    inline std::vector<double> getLikelihoods(){return(likelihoods);};///< gets the likelihoods
    double calcLikelihood(std::vector<std::vector<double> > * P);///< Calculate the likelihood
protected:
    Data *data;
    State *state;
    int verbose;
    std::vector<double> likelihoods;
    long counts;
    my_rng * rng;
};

} // end namespace fines
#endif
