#ifndef INFEXTRACTDONORS_H
#define INFEXTRACTDONORS_H
#include <vector>
#include "state.h"
#include "data.h"
#include "rng.h"
#include "node.h"
#include "fsxml.h"

namespace fines
{

/**
    @brief Inference algorithm: extraction of donor file
*/
class InfExtractDonor
{
public:
    InfExtractDonor(State *s,int v=0);///<Constructor 
    void printDonorFile(std::ostream * out);///
    ~InfExtractDonor();
protected:
    State *state;
    int verbose;
};

} // end namespace fines
#endif
