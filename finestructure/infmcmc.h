#ifndef INFMCMC_H
#define INFMCMC_H
#include <vector>
#include <iomanip>
#include "state.h"
#include "fines.h"
#include "rng.h"

#define MCMCPROB 31.0
#define MCMCNUMVARS 100.0

namespace fines
{

/**
    @brief Inference algorithm using mcmc and SAMS
*/
class InfMCMC
{
public:
  InfMCMC(my_rng * rng,State *s,int datainference=INFDATA_COUNTS,double pcaprob=0, int v=0);///<Constructor if everything already known
    void fixK();///< Fixes the number of populations to the inital number by disabling SAMS moves
    bool moveMergeAndSplit(bool greedy=false);///<merge and split in one move on SAMS like rules, return success or failure
    bool moveSAMS(bool greedy=false);///<merge and split using SAMS.  returns success or failure
    bool moveIndiv(bool greedy=false);///< moves a single random individual between populations
    bool moveDelta(int a,bool greedy=false);///< Moves the delta hyperparameter
    bool moveF(int a,bool greedy=false);///< Moves the F hyperparameter
    bool moveFRef(int a,bool greedy=false);///< Move the F hyperparameter for the reference panel
  bool moveHyperParamLength(int par,bool greedy=false);///<Moves the parameter for lengths from 0-3 as specified
    bool moveHyperParamTotLength(int par,bool greedy=false);/// Moves the parameter for total lengths from 0-1 as specified
    double moveHyper(bool greedy=false);///< Moves all the Fs and deltas
    void runone(long iter,long thin,std::ostream * fileout);/// a single MCMC iteration
    void metropolis(long prevints,long numiters,long thin=1,std::ostream * fileout=NULL,long totallength=-1);///<Metropolis Hastings algorithm
    void hillClimb(long prevints,long numiters, long thin=1,std::ostream * fileout=NULL);///< Greedy hill climbing

    void exportXmlIter(std::ostream * fileout,int iter);///<Saves the state
    inline void exportXmlComment(std::ostream * fileout,std::string comment){
	    *fileout<< "<comment>"<<comment<<"</comment>"<<std::endl;
    };///<Print a comment tag
    void exportXmlHead(std::ostream * fileout,long x=-1,long y=-1,long z=-1);///<Xml start and header
    void exportXmlTail(std::ostream * fileout);///<Xml end
    inline State* getState(){return(state);}///<Gets the state
    inline void resetCounters(){
	numSAMS=0;numIndiv=0;numMergeSplit=0;numHyper=0;
	accSAMS=0;accIndiv=0;accMergeSplit=0;accHyper=0;
    }///<Resets the acceptance rate counters 
    ~InfMCMC();
protected:
    State *state;
    Prior *prior;
    int verbose;
    int silent;
    double logstateprob;
    long numSAMS,numIndiv,numMergeSplit,numHyper;
    long accSAMS,accIndiv,accMergeSplit,accHyper;
    double P_SAMS,P_MERGESPLIT,P_IND;
    int datainference;
    double nullprob;
    my_rng * rng;
    
};

} // end namespace fines
#endif
