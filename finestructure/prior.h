#ifndef PRIOR_H
#define PRIOR_H

#include <vector>
#include <string>

#include "fines.h"


using namespace std;

namespace fines
{
  
  class Prior
  {
    // This is the HYPER prior, that is, all of the parts of the prior that are constant.
    // The state contains all of the non-constant parts
  public:
    Prior();

    // Things to do with alpha
    inline void setAlpha(double alpha){
      this->alpha=alpha;
    };///< Specify the prior on the number of clusters via DP(alpha)

    // Things to do with the prior on data
    void setPriorCounts(int betamodel);///< Specify the model type for the counts 
    void setPriorCounts(string tmp);   ///< Specify the parameters for the counts
    void setPriorCounts(vector<double> bvec);   ///< Specify the parameters for the counts
    void setPriorRef(int betamodel);   ///< Specify the model type for the reference
    void setPriorRef(string tmp);      ///< Specify the parameters for the reference
    void setPriorRef(vector<double> bvec);      ///< Specify the parameters for the reference

    void setPriorLength(int betamodel);///< Specify the model type for the lengths
    void setPriorTotLength(int betamodel);   ///< Specify the model type for the total lengths

    inline int getBetaModel(){ return(betamodel);}
    inline int getBetaModelRef(){ return(betamodelref);}
    inline int getBetaModelLength(){ return(betamodellength);}
    inline int getBetaModelTotLength(){ return(betamodeltotlength);}

    inline vector<double> getBetaPar(){ return(bvec);}
    inline vector<double> getBetaParRef(){ return(bvecref);}
    inline vector<double> getBetaParLength(){ return(bveclength);}
    inline vector<double> getBetaParTotLength(){ return(bvectotlength);}

    inline double getAlpha(){return(alpha);}
    ~Prior();
  private:
    double alpha;
    int betamodel;
    int betamodelref;
    int betamodellength; ///<priors for the chunk length model (alpha_0,beta_0,delta_alpha,delta_beta)
    int betamodeltotlength;
    vector<double> bvec;
    vector<double> bvecref;    
    vector<double> bveclength;
    vector<double> bvectotlength;
  }; ///* The fixed, top level priors and the information needed to construct priors of all levels
}
#endif
