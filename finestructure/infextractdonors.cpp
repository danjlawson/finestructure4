#include "infextractdonors.h"
#include "data.h"
#include <iomanip>
#include <iostream>

using namespace std;
namespace fines
{

  InfExtractDonor::InfExtractDonor(State *s,int v)
  {
    this->state=s;
    this->verbose=v;
  }


  void InfExtractDonor::printDonorFile(std::ostream *out)
{
  /* It would be useful to be able to do something sensible to make a) minimum population sizes, and b) remove admixed individuals as donors*/
  if(verbose) cout<<"Extracting donors from state with "<<state->getP()<<" Populations"<<endl;
  Data * data=state->getData();
  int n =state->getN();
  for(int ind=0;ind<n;++ind){
    *out<<data->getnames(ind)<<"\tPop"<<state->getPop(ind)<<"\t1"<<endl;
  }
}

InfExtractDonor::~InfExtractDonor()
{
}


} // end namespace fines
