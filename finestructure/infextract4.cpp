#include "infextract4.h"
#include <iomanip>
#include <iostream>

using namespace std;
namespace fines
{

  InfExtract4::InfExtract4(my_rng * rng,Data *d,FsXml *fs,Prior *prior,int v)
{
	data=d;
	verbose=v;
	streampos fpos;
	State *tstate=NULL;
	counts=0;
	if(verbose) cout<<"Creating posterior probability list ..."<<endl;
	fpos=fs->gotoLineContaining("<Iteration>");
	while(!fs->eof() && fpos>=0) {
	  if(fpos>0) tstate=new  State(rng,data,fs,prior,true);
		posteriors.push_back(tstate->posteriorProb());
		fpos=fs->gotoNextLineContaining("<Iteration>");
		counts++;
	}
	if(verbose) cout<<"done."<<endl;
}

InfExtract4::~InfExtract4()
{
}


} // end namespace fines
