#include "infextract2.h"
#include <iomanip>
#include <iostream>

using namespace std;
namespace fines
{

  InfExtract2::InfExtract2(my_rng * rng,Data *d,FsXml *fs,Prior *prior,int v,Data *dlength,Data *dref,int datainference,int modeltype)
{
  this->rng=rng;
	data=d;
	verbose=v;
	streampos fpos;
	State *tstate=NULL;
	if(verbose) cout<<"Creating highest posterior state..."<<endl;
	fpos=fs->gotoLineContaining("<Iteration>");
	if(fpos<0){ // its a population file, not a complete xml file
	  state=new State(rng,d,fs,prior,true,dlength,dref,datainference,modeltype);
		if(verbose) cout<<"InfExtract2: Extracting single population"<<endl;
	}else{// is an xml output file
		int counts=0,maxit=0;
		double bestpost=0,curpost;
		while(!fs->eof() && fpos>=0) {
		  if(fpos>0) tstate=new State(rng,data,fs,prior,true,dlength,dref,datainference,modeltype);
			
			curpost=tstate->posteriorProb();
			if(counts==0 || curpost> bestpost){
				if(counts>0) delete(state);
				state=tstate;
				bestpost=curpost;
				maxit=counts;
			}else {
				delete(tstate);
			}
			counts++;
			fpos=fs->gotoNextLineContaining("<Iteration>");
		}
		if(verbose) cout<<setprecision(32)<<"Best state is "<<maxit<<" with posterior "<<bestpost<<endl;
	}
	if(verbose) cout<<"done"<<endl;
}

InfExtract2::~InfExtract2()
{
}


} // end namespace fines
