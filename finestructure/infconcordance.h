#ifndef INFCONCORDANCE_H
#define INFCONCORDANCE_H
#include <vector>
#include "math.h"
#include "state.h"
#include "prior.h"
#include "rng.h"
#include "node.h"
#include "infextract.h"

namespace fines
{

/**
    @brief Inference algorithm for maximising the mcmc samples co-occurance concordance
*/
class InfConcordance
{
public:
  InfConcordance(my_rng * rng,
		 State *refstate,
		 int maxiters,
		 Data *d,
		 FsXml *fs,
		 Prior *prior,
		 int v=0,
		 Data *dlength=NULL,
		 Data *dref=NULL,
		 int datainference=INFDATA_COUNTS,
		 int modeltype=MODELTYPE_FINESTRUCTURE);///<Constructor 
    ~InfConcordance();
    inline State * getState(){return(refstate);}///< Returns the state
    vector< vector <double> > obtainScoreMatrix(State *refstate);///<calculate the score matrix for a given reference state
    int updateState();///< Update the state given a computed score matrix
    vector <int> getAssignments(State * state);///< get the assignments of individuals into populations
    void runUntilConvergence();///< run the iterative process
    void updateReferenceState(); ///< update the reference state given the score matrix
    vector< vector<double> > calculateScore(State *state);///< Calculate the score matrix as used in each iteration of the concordance step
    inline vector< vector <double> > getScoreMatrix(){
      return(refscore);
    }
    void printScoreMatrix(ostream * out);///< Prints the score matrix
protected:
    my_rng * rng;
    Data *data;
    FsXml *fs;
    Prior *prior;
    Data *dlength;
    Data *dref;
    int datainference;
    int modeltype;
    int verbose;

    State *refstate; ///< The current reference state, which starts out given and ends up as our best state
    int maxiters; ///< maximum number of times we move individuals around
    int iteron; ///< Iterations actually used
    int converged; ///< Whether we've converged yet
    vector < vector <double> > refscore;///< Scores for each assignment at the last iteration to be processed
    vector <int> refassignments; ///< Assignments of individuals to populations
};

} // end namespace fines
#endif
