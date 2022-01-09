#ifndef INFEXTRACT3_H
#define INFEXTRACT3_H
#include <vector>
#include "state.h"
#include "data.h"
#include "rng.h"
#include "node.h"
#include "fsxml.h"

namespace fines
{

/**
    @brief Inference algorithm: extraction of tree split support
*/
class InfExtract3
{
public:
    InfExtract3(my_rng * rng,Data *d,FsXml *fs,std::vector<Node*> nodes, int v=0);///<Constructor 
    ~InfExtract3();
    double getSingleSplit(State *state, std::vector<int> tips);///< returns the proportion of agreement between the tips can be constructed out of populations in the state, 0 if they can't
    double getSingleSplitStrict(State *state, std::vector<int> tips);///< returns 1 if the tips can be constructed out of populations in the state, 0 if they can't
    std::vector<double> getSplits(State *state, std::vector<Node*> nodes,bool strict=true);///< gets the vector of splits for all nodes above startnode
    inline State * getState(){return(state);}///< Returns the state
protected:
    Data *data;
    State *state;
    int verbose;
    std::vector<double> weights;
    std::vector<double> weights2;
    int counts;
    my_rng * rng;
};

} // end namespace fines
#endif
