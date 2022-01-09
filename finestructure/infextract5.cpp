#include "infextract5.h"
#include <iomanip>
#include <iostream>

using namespace std;
namespace fines
{

  InfExtract5::InfExtract5(my_rng * rng,Data *d,FsXml *fs,Prior *prior,long reps,int v)
{
	data=d;
	verbose=(v==1);
	streampos fpos;
	counts=0;
	if(verbose) cout<<"Creating likelihood sample ..."<<endl;
	fpos=fs->gotoLineContaining("<Iteration>");
	while(!fs->eof() && fpos>=0) {
	  if(fpos>0) state=new  State(rng,data,fs,prior,true);
//cout<<"state beta[0,1]="<<state->getBeta(0,1)<<endl;
		InfAdmixture infad(rng,d,state,1.0,verbose,0.0);
		for(long c1=0;c1<reps;c1++){
		  infad.samplePs(false);
//		  infad.printPs(&cout);
		  likelihoods.push_back(calcLikelihood(infad.getP()));
//		  likelihoods.push_back(state->admixtureLogLikelihood(infad.getQ(),infad.getP(),infad.getQsums()));
		}
		fpos=fs->gotoNextLineContaining("<Iteration>");
		counts+=reps;
	}
	if(verbose) cout<<"done."<<endl;
}

double InfExtract5::calcLikelihood(std::vector<std::vector<double> > * P)
{
	double ret=0.0;
/*	for(unsigned int i=0;i<data->getDim();i++){
	  for(unsigned int j=0;j<data->getDim();j++){
		double nhat=state->getPsize(state->getPop(j));
		if(state->getPop(i)==state->getPop(j))nhat-=1.0;
		if(nhat>0) ret+=log(P->at(state->getPop(i))[state->getPop(j)]/nhat)*data->get(i,j,true);
	  }
	}*/
	for(int a=0;a<state->getP();a++){
	  for(int b=0;b<state->getP();b++){
		double nhat=state->getPsize(b);
		if(a==b)nhat-=1.0;
		if(nhat>0) ret+=log(P->at(a)[b]/nhat)*state->sumXab(a,b);
	  }
	}
	return(ret);
}

InfExtract5::~InfExtract5()
{
}


} // end namespace fines
