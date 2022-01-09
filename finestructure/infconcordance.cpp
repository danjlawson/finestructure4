#include "infconcordance.h"
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <cassert>

#define MIN_LINK 0.00001
#define BASE_LINK 1
using namespace std;
namespace fines
{

  InfConcordance::InfConcordance(my_rng * rng, State *refstate, int maxiters,Data *d,FsXml *fs,Prior *prior,int v,Data *dlength,Data *dref,int datainference,int modeltype)
{
  // 
  this->rng=rng;
  this->data=d;
  this->verbose=v;
  this->fs = fs;
  this->prior=prior;
  this->verbose = v;
  this->dlength = dlength;
  this->dref = dref;
  this->datainference = datainference;
  this->modeltype = modeltype;

  this->refstate = refstate;
  this->maxiters = maxiters;

  this->iteron=0;
  this->converged=0;
  runUntilConvergence();
}
  
vector< vector <double> > InfConcordance::obtainScoreMatrix(State *refstate) {
  // Takes a reference state and computes the score matrix from it
  streampos fpos;
  State *tstate=NULL;
  vector< vector<double> > totalscore;
  fs->rewind();
  fpos=fs->gotoLineContaining("<Iteration>");
  if(fpos<0){ // its a population file, not a complete xml file
    // ERROR: Not a valid mcmc file
    cerr<<"ERROR: Cannot run concordance calculation on this xml file, which does not appear to contain any iterations."<<endl;
    exit(0);
  }
  
  int counts=0,maxit=0;
  double bestpost=0,curpost;
  while(!fs->eof() && fpos>=0) {
    if(fpos>0) tstate=new State(rng,data,fs,prior,true,dlength,dref,datainference,modeltype);
    
    vector< vector<double> > currentscore = calculateScore(tstate); 
    // if(verbose){
    //   cout<<"DEBUG: Iteration "<<counts<<" Counted "<<currentscore.size()<<" individuals with ";
    //   if(currentscore.size()>0) {cout<<currentscore[0].size()<<" populations."<<endl;
    //   }else cout<<"no populations."<<endl;
    // }
    for(unsigned int c1=0;c1<currentscore.size();++c1){
      while(totalscore.size()<=c1) totalscore.push_back(vector<double>());
      for(unsigned int c2=0;c2<currentscore[c1].size();++c2){
	while(totalscore[c1].size()<=c2) totalscore[c1].push_back(0);
	totalscore[c1][c2]+=currentscore[c1][c2];
      }
    }

    delete(tstate);
    counts++;
    fpos=fs->gotoNextLineContaining("<Iteration>");
  }
  if(verbose){
    cout<<" Counted "<<totalscore.size()<<" individuals with ";
    if(totalscore.size()>0) cout<<totalscore[0].size()<<" populations."<<endl;
    else cout<<"no populations."<<endl;
  }
  // normalise each entry by the number of mcmc iterations
  for(unsigned int c1=0;c1<totalscore.size();++c1){
    for(unsigned int c2=0;c2<totalscore[c1].size();++c2){
      totalscore[c1][c2]/=counts;
    }
  }
  
  return(totalscore);
}

int InfConcordance::updateState(){
  // obtain the score matrix, move individuals appropriately, and then report whether the state has changed.
  vector <int> before=getAssignments(refstate); // population assignments before we change anything
  refscore=obtainScoreMatrix(refstate); // updated score matrix
  updateReferenceState(); // update the reference state by putting everyone in their favorite state
  refassignments= getAssignments(refstate); // get the population assignments now
  
  assert(before.size()==refassignments.size());
  
  for(unsigned int i=0;i<before.size();++i){
    if(before[i]!=refassignments[i]) return(0);
  }
  return(1);
}

vector <int> InfConcordance::getAssignments(State * state){
  vector <int> ret;
  for(int i=0;i<state->getDim();++i){
    ret.push_back(state->getPop(i));
  }
  return(ret);
}

void InfConcordance::runUntilConvergence(){
  if(verbose) cout<<"Running MCMC Concordance state calculation with up to "<<maxiters<<" iterations."<<endl;
  while((!converged) && (iteron<maxiters)){
    if(verbose) cout<<"Iteration "<<iteron<<endl;
    converged  = updateState();
    ++iteron;
  }
  if(converged){
    if(verbose) cout<<"Converged after "<<iteron<<" Iterations."<<endl;
  }else{
    cout<<"WARNING: Concordance calculation not converged after "<<iteron<<" Iterations, but terminating anyway!"<<endl;
  }
  
}

vector< vector<double> > InfConcordance::calculateScore(State *state){
  /// calculates the score used in the pobi final step
  // The score is y_{ik}/x_i for ind i in pop k
  // where x_i is the number of other individuals in the same population as i in this mcmc sample state
  // we avoid explicitly storing y_{ik} = the number of such individuals that are also in the population k in the reference state
  // then ret = y_{ik}/x_i is returned.
  vector <double> x; // 
  vector <vector<double> > ret;
  for(int i=0;i<state->getDim();++i){
    x.push_back(state->getPsize(state->getPop(i)));
    ret.push_back(vector<double>(refstate->getP(),0.0) );

    vector<int> otherinds=state->getIndInPop(state->getPop(i));
    for(unsigned int ti=0;ti<otherinds.size();++ti){
      int otherspop =refstate->getPop(otherinds[ti]);
      ret[i][otherspop]+=1.0/x[i];
    }
  }
  return(ret);
}

void InfConcordance::updateReferenceState(){
  // Updates the reference state according to the refscore, by putting everyone in their favorite populations.
  //  cout<<"DEBUG: refscore="<<refscore.size()<<" Dim ="<<refstate->getDim()<<endl;
  assert((int)refscore.size()==refstate->getDim());
  
  for(unsigned int i=0;i<refscore.size();++i){
    //    cout<<"DEBUG: refscore["<<i<<".size()="<<refscore[i].size()<<endl;
    // for(unsigned int c1=0;c1<refscore[i].size();++c1){
    //   cout<<refscore[i][c1]<<",";
    // }
    // cout<<endl;
    int newp=distance(refscore[i].begin(),
		      max_element(refscore[i].begin(),
				  refscore[i].begin()+refscore[i].size()));
    if(refstate->getPop(i) != newp){
      //      cout<<"DEBUG: Moving individual "<<i<<" to population "<<newp<<endl;
      refstate->moveInd(i,newp);
    }
  }
  //  cout<<"DEBUG: Removing empty populations"<<endl;
  refstate->removeEmptyPops();
  //  cout<<"DEBUG: Done removing"<<endl;
}

  void InfConcordance::printScoreMatrix(ostream * out){
    Data *d=data;
    if(d==NULL) d=dlength;
    if(d==NULL) d=dref;
    if(d==NULL) {
      cerr<<"ERROR: No data defined!"<<endl;
      exit(1);
    }
    //    if(refscore.size()==0)
    refscore=obtainScoreMatrix(refstate);
    *out<<"COINCIDENCE";
    if(refscore.size()>0) for(unsigned int i=0;i<refscore[0].size();i++) {
      *out<<", POP"<<i;
    }
    *out<<endl;
    
    for(unsigned int i=0;i<refscore.size();i++) {
      *out<<d->getnames(i);
      for(unsigned int j=0;j<refscore[i].size();j++) {
	*out<<", "<<refscore[i][j];
      }
      *out<<endl;
    }
  }
  
InfConcordance::~InfConcordance()
{
}


} // end namespace fines
