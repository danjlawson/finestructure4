#include "infmcmc.h"
#include "prior.h"
#include "math.h"

using namespace std;
namespace fines
{

  InfMCMC::InfMCMC(my_rng * rng,State *s,int datainference,double pcaprob, int v)
{
  this->rng=rng;
	this->datainference=datainference;
	state=s;
	verbose= (v==1);
	silent=0;if(v>1)silent=1;
	logstateprob=state->posteriorProb();
	resetCounters();
	P_SAMS=MCMCPROB/MCMCNUMVARS;
	P_MERGESPLIT=P_SAMS+MCMCPROB/MCMCNUMVARS;
	P_IND=P_MERGESPLIT+MCMCPROB/MCMCNUMVARS;
	//P_HYPER =1-sum(these)
	this->prior=state->getPrior();
}

  void InfMCMC::fixK()
  {
	P_SAMS=0;
	P_MERGESPLIT=P_SAMS+MCMCPROB/(MCMCNUMVARS-1.0);
	P_IND=P_MERGESPLIT+MCMCPROB/(MCMCNUMVARS-1.0);
    
  }


bool InfMCMC::moveMergeAndSplit(bool greedy)
{
	if(state->getP()==1) return(0);// Cannot merge only one population
	double logpofsplit1=0.0,logpofsplit2=0.0;
	double logposteriorold=logstateprob,logposteriornew;
	int i=RandomInteger(rng,0,state->getDim()-1),j=i;
	// force a merge of two different populations
	while(i==j || state->getPop(i)==state->getPop(j)) j=RandomInteger(rng,0,state->getDim()-1);
	if(verbose){cout<<"MergeAndSplit: INITIAL STATE:"<<endl;
		state->setprint(&cout);}
	if(state->getPlength(state->getPop(i))==1 && state->getPlength(state->getPop(j))==1){
		if(verbose) {cout<<"Rejected: Single size populations"<<endl;
			state->setprint(&cout);}
		return(0);// will always give identical populations
	}
 	logpofsplit1= state->probOfSplitSAMS(i,j);
	State *oldstate=new State(state);
	state->merge(state->getPop(i),state->getPop(j));
	logpofsplit2 = state->splitSAMS(i,j,false);
	logposteriornew=state->posteriorProb();

	  if(verbose){cout<<"PROPOSED STATE:"<<endl;
		state->setprint(&cout);
	  	cout<<"PROPOSED MERGE&SPLIT LOG PROB="<<logpofsplit1<<"-"<<logpofsplit2<<endl;
	  }

	if((!greedy && log(rnd(rng))<logposteriornew-logposteriorold + logpofsplit2-logpofsplit1 && fabs(logposteriorold-logposteriornew)>1e-07) || (greedy && logposteriornew>logposteriorold)) {
// insist we only accept if the states are different
		if(verbose) {cout<<"Accepted:"<<"LL="<<logposteriornew<<" versus oldLL="<<logposteriorold<<endl;
			state->setprint(&cout);}
		logstateprob=logposteriornew;
		delete(oldstate);
		return(true);
	}else{
		if(verbose) {cout<<"Rejected:"<<"LL="<<logposteriornew<<" versus oldLL="<<logposteriorold<<endl;
			state->setprint(&cout);}
		delete(state);
		state=oldstate;
		return(false);
	}
}

bool InfMCMC::moveSAMS(bool greedy)
{
	if(state->getDim()==1) return(0);// cannot merge-split on one individual
	double logpofsplit=0.0;
	double logposteriorold=logstateprob,logposteriornew;
// Choose two random individuals
	int i=RandomInteger(rng,0,state->getDim()-1),j=i;
	while(i==j) j=RandomInteger(rng,0,state->getDim()-1);
// Split or merge as appropiate
	if(verbose){cout<<"Move SAMS BEFORE MOVE"<<endl;state->setprint(&cout);}
 	if(state->getPop(i)==state->getPop(j)) {// try a split
	  logpofsplit = state->splitSAMS(i,j,false);
	  logposteriornew=state->posteriorProb();
	  if(verbose){cout<<"PROPOSED STATE:"<<endl;
		state->setprint(&cout);
	  	cout<<"PROPOSED SPLIT WITH LOG PROB="<<logpofsplit<<endl;
	  }
	  if((!greedy && log(rnd(rng))<logposteriornew-logposteriorold - logpofsplit) || (greedy && logposteriornew>logposteriorold)) {// accept split
		if(verbose) {cout<<"Accepted:"<<"LL="<<logposteriornew<<" versus oldLL="<<logposteriorold <<"PofSplit(-)"<<logpofsplit <<endl;
			state->setprint(&cout);}
		logstateprob=logposteriornew;
		return(true);
	  }else {
		state->merge(state->getPop(i),state->getPop(j));
		if(verbose) {cout<<"Rejected:"<<"LL="<<logposteriornew<<" versus oldLL="<<logposteriorold <<"PofSplit(-)"<<logpofsplit <<endl;
			state->setprint(&cout);}
		return(false);
	  }
	}else { // try a merge
	  logpofsplit= state->probOfSplitSAMS(i,j);
	  State *oldstate=new State(state);
	  state->merge(state->getPop(i),state->getPop(j));
	  logposteriornew=state->posteriorProb();
	  if(verbose){cout<<"PROPOSED STATE:"<<endl;
		state->setprint(&cout);
	  	cout<<"PROPOSED MERGE WITH split LOG PROB="<<logpofsplit<<endl;
	  }
	  if(log(rnd(rng))<logposteriornew-logposteriorold + logpofsplit || (greedy && logposteriornew>logposteriorold)) {// accept merge
		logstateprob=logposteriornew;
		if(verbose) {cout<<"Accepted:"<<"LL="<<logposteriornew<<" versus oldLL="<<logposteriorold <<"PofSplit(+)"<<logpofsplit <<endl;
			state->setprint(&cout);}
		delete(oldstate);
		return(true);
	  }else {
		delete(state);
		state=oldstate;
		if(verbose) {cout<<"Rejected:"<<"LL="<<logposteriornew<<" versus oldLL="<<logposteriorold <<"PofSplit(+)"<<logpofsplit <<endl;
			state->setprint(&cout);}
		return(false);
	  }
	}
	// should never get here:
	return(false);
}

bool InfMCMC::moveIndiv(bool greedy)
{
	if(state->getP()==1) return(false);// can't move if only one population
// Choose a random individual
	int i=RandomInteger(rng,0,state->getDim()-1);
	if(state->getPlength(state->getPop(i))==1) return(false);//cant undo moving single individual; thats a split/merge move
	if(verbose) {cout<<"MOVEINDIV: BEFORE ";state->setprint(&cout);}
	  
// choose a random other population
	int oldpop = state->getPop(i),newpop=oldpop;
	while(newpop==oldpop) newpop=RandomInteger(rng,0,state->getP()-1);

	vector <int> tmppop=state->getIndInPop(state->getPop(i));
	int otherindina=-1;
	for(int j=0;j<(int)tmppop.size();j++) if(tmppop[j]!=i) {otherindina=tmppop[j];break;}
	if(otherindina<0) throw(string("Error in moveIndiv: Can't find any other individuals!"));
	State *mergestateold=new State(state);
	mergestateold->merge(mergestateold->getPop(i),newpop);
	
	double logposteriorold=state->posteriorProb(mergestateold);
	//double logposteriorold=state->posteriorProb(NULL);
	state->moveInd(i,newpop);
	if(verbose) {cout<<"MOVEINDIV: PROPOSED ";state->setprint(&cout);}
	State *mergestatenew=new State(state);
	mergestatenew->merge(mergestatenew->getPop(otherindina),newpop);
	double logposteriornew=state->posteriorProb(mergestatenew);
	//double logposteriornew=state->posteriorProb(NULL);
	delete(mergestateold);
	delete(mergestatenew);
// Acceptance/rejection step
	if(verbose) {cout<<"MOVEINDIV: log(p) = "<<logposteriornew<<" - "<<logposteriorold<<endl;}
	if( (!greedy && log(rnd(rng))< logposteriornew-logposteriorold) || (greedy && logposteriornew>logposteriorold)) {
		//accept
	if(verbose) cout<<"MOVEINDIV ACCEPTED"<<endl;
		logstateprob=logposteriornew;
		return(true);
	}else {
	if(verbose) cout<<"MOVEINDIV REJECTED"<<endl;
		state->moveInd(i,oldpop);
		return(false);
	}
}

bool InfMCMC::moveF(int a,bool greedy)
{
  if(verbose) cout<<"MOVE F"<<endl;
    double d0=state->getBetaF(a);
    double d1=exp(log(d0) + (rnd(rng)*2.0-1.0)*1);
    while(d1<=0 || d1>=1) {
	if(d1<0) d1=-d1;
	if(d1>=1) d1=2-d1;
	if(d1==0) d1=rnd(rng)*0.001;
	if(d1==1) d1=1.0*rnd(rng) *0.001;
    }
    double p0=state->posteriorProb();
    state->setBetaF(a,d1);
    double p1=state->posteriorProb();
    double postdiff=p1-p0 + state->LGammaDistProb(prior->getBetaPar()[0],prior->getBetaPar()[2],d1) - state->LGammaDistProb(prior->getBetaPar()[0],prior->getBetaPar()[2],d0);
    if((!greedy && log(rnd(rng))< postdiff) || (greedy &&postdiff>0)) {
	if(verbose) cout<<"MOVEF ACCEPTED: from "<<d0<<" to "<<d1<<endl;
	return(true);
    }else {
	state->setBetaF(a,d0);
	if(verbose) cout<<"MOVEF REJECTED"<<endl;
	return(false);
    }
}

  bool InfMCMC::moveFRef(int a,bool greedy)
{
  if(verbose) cout<<"MOVE F Reference"<<endl;
    double d0=state->getRefBetaF(a);
    double d1=exp(log(d0) + (rnd(rng)*2.0-1.0)*1);
    while(d1<=0 || d1>=1) {
	if(d1<0) d1=-d1;
	if(d1>=1) d1=2-d1;
	if(d1==0) d1=rnd(rng)*0.001;
	if(d1==1) d1=1.0*rnd(rng) *0.001;
    }
    double p0=state->posteriorProb();
    state->setRefBetaF(d1);
    double p1=state->posteriorProb();
    double postdiff=p1-p0 + state->LGammaDistProb(prior->getBetaParRef()[0],prior->getBetaParRef()[1],d1) - state->LGammaDistProb(prior->getBetaParRef()[0],prior->getBetaParRef()[1],d0);
    if((!greedy && log(rnd(rng))< postdiff) || (greedy &&postdiff>0)) {
	if(verbose) cout<<"MOVEFREF ACCEPTED: from "<<d0<<" to "<<d1<<endl;
	return(true);
    }else {
	state->setRefBetaF(d0);
	if(verbose) cout<<"MOVEFREF REJECTED"<<endl;
	return(false);
    }
}

bool InfMCMC::moveDelta(int a,bool greedy)
{
//   state->setDelta(a,state->sampleBetaP(true));
//   return(true);
  if(verbose) cout<<"MOVE DELTA"<<endl;
    double d0=state->getDelta(a);
    double d1=exp(log(d0) + (rnd(rng)*2.0-1.0)*1);
    if(d1<0) d1=-d1;
    double p0=state->posteriorProb();
    state->setDelta(a,d1);
    double p1=state->posteriorProb();
    double postdiff=p1-p0 + state->LGammaDistProb(prior->getBetaPar()[1],prior->getBetaPar()[3],d1) - state->LGammaDistProb(prior->getBetaPar()[1],prior->getBetaPar()[3],d0);
    if((!greedy && log(rnd(rng))< postdiff)||(greedy && postdiff>0)) {
	if(verbose) cout<<"MOVEDELTA ACCEPTED: from "<<d0<<" to "<<d1<<endl;
	return(true);
    }else {
	state->setDelta(a,d0);
	if(verbose) cout<<"MOVEDELTA REJECTED"<<endl;
	return(false);
    }
}

bool InfMCMC::moveHyperParamLength(int par,bool greedy)
{
    if(verbose) cout<<"MOVE HYPER PARAM LENGTH"<<endl;
    
    double d0=state->getpriorLengthsParam(par);
    double d1=exp(log(d0) + (rnd(rng)*2.0-1.0)*1);
    if(d1<0) d1=-d1;
    double p0=state->posteriorProb();
    state->setpriorLengthsParam(par,d1);
    double p1=state->posteriorProb();
    vector<double>hypervec=prior->getBetaParLength();
    double postdiff=p1-p0 + state->LGammaDistProb(hypervec[par],hypervec[par],d1) - state->LGammaDistProb(hypervec[par],hypervec[par],d0);
    if((!greedy && log(rnd(rng))< postdiff)||(greedy && postdiff>0)) {
	if(verbose) cout<<"MOVEHYPERPARAMLENGTH "<<par<<" ACCEPTED: from "<<d0<<" to "<<d1<<endl;
	return(true);
    }else {
	state->setpriorLengthsParam(par,d0);
	if(verbose) cout<<"MOVEHYPERPARAMLENGTH "<<par<<" REJECTED"<<endl;
	return(false);
    }
}

bool InfMCMC::moveHyperParamTotLength(int par,bool greedy)
{
  if(verbose) cout<<"MOVE HYPER PARAM TOT LENGTH"<<endl;

  double d0=state->getpriorNumChunksParam(par);
    double d1=exp(log(d0) + (rnd(rng)*2.0-1.0)*1);
    if(d1<0) d1=-d1;
    double p0=state->posteriorProb();
    state->setpriorNumChunksParam(par,d1);
    double p1=state->posteriorProb();
    vector<double>hypervec=prior->getBetaParTotLength();
    double postdiff=p1-p0 + state->LGammaDistProb(hypervec[par],hypervec[par],d1) - state->LGammaDistProb(hypervec[par],hypervec[par],d0);
    if((!greedy && log(rnd(rng))< postdiff)||(greedy && postdiff>0)) {
	if(verbose) cout<<"MOVEHYPERPARAMTOTLENGTH "<<par<<" ACCEPTED: from "<<d0<<" to "<<d1<<endl;
	return(true);
    }else {
	state->setpriorNumChunksParam(par,d0);
	if(verbose) cout<<"MOVEHYPERPARAMTOTLENGTH "<<par<<" REJECTED"<<endl;
	return(false);
    }
}

double InfMCMC::moveHyper(bool greedy)
{
  if(verbose) cout<<"MOVE HYPER"<<endl;

  double pacc=0;
	double denom=0;
	vector<double> hypercounts,hyperref,hyperlen,hypersums;
	if(state->useCounts()) {
	  hypercounts=prior->getBetaPar();
	  for(unsigned int c1=0;c1< hypercounts.size()/2;c1++) if(hypercounts[c1]>0) denom++;
	}
	if(state->useRef()) {
	  hyperref=prior->getBetaParRef();
	  for(unsigned int c1=0;c1< hyperref.size()/2;c1++) if(hyperref[c1]>0) denom++;
	}
	//////////// *****************************
	// I have a suspicion that these have thr wrong length...
	if(state->useLengths()) {
	  hyperlen=prior->getBetaParLength();
	  for(unsigned int c1=0;c1< hyperlen.size();c1++) if(hyperlen[c1]>0) denom++;
	}
	if(state->useSums()) {
	  hypersums=prior->getBetaParTotLength();
	  for(unsigned int c1=0;c1< hypersums.size();c1++) if(hypersums[c1]>0) denom++;
	}

	
	if(state->useCounts()){
	  if(prior->getBetaModel()==BETAMOD_F2) {
	    if(hypercounts[0]>0) pacc+=moveF(0)/denom;
	    if(hypercounts[1]>0) pacc+=moveDelta(0)/denom; 
	  }
	  else if(prior->getBetaModel()==BETAMOD_F) {
	    for(unsigned int i=0;i<hypercounts.size()/2;i++) {
	      if(hypercounts[0]>0) pacc+=moveF(i)/denom;
	      if(hypercounts[1]>0) pacc+=moveDelta(i)/denom; 
	    }
	  }
	}// end if state->useCounts

	if(state->useRef()){
	  if(prior->getBetaModelRef()==BETAMOD_F2) {
	    if(hyperref[0]>0) pacc+=moveFRef(0)/denom; 
	  }
	}
	if(state->useLengths()){
	  for(unsigned int i=0;i<hyperlen.size();i++)  if(hyperlen[i]>0) pacc+=moveHyperParamLength(i)/denom; 
	}
	if(state->useSums()){
	  for(unsigned int i=0;i<hypersums.size();i++) if(hypersums[i]>0) pacc+=moveHyperParamTotLength(i)/denom;
	}
	return(pacc);
}

void InfMCMC::runone(long iter,long thin,std::ostream * fileout){
  	double x=rnd(rng);
	if(x<P_SAMS) {
		numSAMS++;
		accSAMS+=(int)moveSAMS();
	}else if(x<P_MERGESPLIT) {
		numMergeSplit++;
		accMergeSplit+=(int)moveMergeAndSplit();
	}else if(x<P_IND) {
		numIndiv++;
		accIndiv+=moveIndiv();
	}else {
		numHyper++;
		accHyper+=(int)moveHyper();
	}
	cout<<"DEBUG C1"<<endl<<flush;
	cerr<<"CERR DEBUG C1"<<endl<<flush;
	if(iter % thin==0 && fileout!=NULL) exportXmlIter(fileout,iter);
	cerr<<"CERR DEBUG C2"<<endl<<flush;
	cout<<"DEBUG C2"<<endl<<flush;
}

void InfMCMC::metropolis(long prevints,long numints, long thin,std::ostream * fileout,long totallength)
{
   if(totallength<0)totallength=numints;
	for(long c1=prevints;c1<prevints+numints;c1++) {
		if (totallength>50 && (c1)%((totallength)/50)==0)
		{
			if (100l*c1/(totallength)<=1){
			  if(!silent) cout<<"#  "<<100l*(double)c1/(numints)<<"%"<<flush;
			}else if (100l*c1/(totallength)<10){
				if(!silent) cout<<"\b\b\b\b#  "<<(int)(100l*(double)c1/(totallength))<<"%"<<flush;
			}else
				if(!silent) cout<<"\b\b\b\b# "<<(int)(100l*(double)c1/(totallength))<<"%"<<flush;
		}
		if (c1+1==totallength)
			if(!silent) cout<<"\b\b\b\b# 100%"<<endl<<flush;
		cerr<<"CERR DEBUG B1"<<endl<<flush;
		cout<<"DEBUG B1"<<endl<<flush;
		runone(c1,thin,fileout);
		cerr<<"CERR DEBUG B2"<<endl<<flush;
		cout<<"DEBUG B2"<<endl<<flush;
	}
	cerr<<"CERR DEBUG B3"<<endl<<flush;
	cout<<"DEBUG B3"<<endl<<flush;
}


void InfMCMC::hillClimb(long prevints,long numints, long thin,std::ostream * fileout)
{
	for(long c1=0;c1<numints;c1++) {
		if (numints>50 && (c1)%((numints)/50)==0)
		{
			if (100l*c1/(numints)<=1){
				cout<<"#  "<<100l*(double)c1/(numints)<<"%"<<flush;
			}else if (100l*c1/(numints)<10)
				cout<<"\b\b\b\b#  "<<(int)(100l*(double)c1/(numints))<<"%"<<flush;
			else
				cout<<"\b\b\b\b# "<<(int)(100l*(double)c1/(numints))<<"%"<<flush;
		}
		if (c1+1==numints)
			cout<<"\b\b\b\b# 100%"<<endl<<flush;
		double x=rnd(rng);
		if(x<P_SAMS) {
			numSAMS++;
			accSAMS+=(int)moveSAMS(true);
		}else if(x<P_MERGESPLIT) {
			numMergeSplit++;
			accMergeSplit+=(int)moveMergeAndSplit(true);
		}else if(x<P_IND) {
			numIndiv++;
			accIndiv+=moveIndiv(true);
		}else {
			numHyper++;
			accHyper+=(int)moveHyper(true);
		}
		if(c1 % thin==0 && fileout!=NULL) exportXmlIter(fileout,c1+prevints);
	}
}

  void InfMCMC::exportXmlHead(std::ostream * fileout,long x,long y,long z)
{
	*fileout<<"<?xml version = '1.0' encoding = 'UTF-8'?>"<<endl<<"<outputFile>"<<endl;
  *fileout<< "<header>"<<endl;
  *fileout<< "<runtype>MCMC</runtype>"<<endl;
  *fileout<<setprecision(12);
  if(getState()->getData()!=NULL) *fileout<< "<inflation>"<<getState()->getData()->getCfactor()<< "</inflation>"<<endl;
  if(getState()->getRefData()!=NULL) *fileout<< "<inflationRef>"<<getState()->getRefData()->getCfactor()<< "</inflationRef>"<<endl;
  if(x>=0) *fileout<< "<burnin>"<<x<< "</burnin>"<<endl;
  if(y>=0) *fileout<< "<mcmclength>"<<y<< "</mcmclength>"<<endl;
  if(z>=0) *fileout<< "<skip>"<<z<< "</skip>"<<endl;
  *fileout<<"<datafilename>"<<state->getFileNames()<<"</datafilename>"<<endl;
  *fileout<< "<copymodel>"<<state->getModelType()<< "</copymodel>"<<endl;    
  *fileout<<"</header>"<<endl;	
}

void InfMCMC::exportXmlTail(std::ostream * fileout)
{
	*fileout<<"</outputFile>"<<endl;
}

void InfMCMC::exportXmlIter(std::ostream * fileout,int iter)
{
//	state->printBeta(&cout);
	*fileout<<"<Iteration>"<<endl;
	state->setprint(fileout);
	*fileout<<setprecision (32)<<"<Posterior>"<<logstateprob<<"</Posterior>"<<setprecision (6)<<endl;
	*fileout<<"<K>"<<state->getP()<<"</K>"<<endl;
	*fileout<<"<alpha>"<<state->getAlpha()<<"</alpha>"<<endl;
	if(state->useCounts()){
	  *fileout<<"<beta>"<<state->getSumBeta(0)<<"</beta>"<<endl;
	  if(prior->getBetaModel()==BETAMOD_F2){
	    *fileout<<"<delta>"<<state->getDelta(0)<<"</delta>"<<endl;
	    *fileout<<"<F>"<<state->getBetaF(0)<<"</F>"<<endl;
	  }
	}
	if(state->useRef()){
	  *fileout<<"<RefBeta>"<<state->getSumRefBeta(0)<<"</RefBeta>"<<endl;
	  if(prior->getBetaModelRef()==BETAMOD_F2){
	    *fileout<<"<RefF>"<<state->getRefBetaF(0)<<"</RefF>"<<endl;
	  }
	}
	if(state->useLengths()){
	*fileout<<"<lengthalpha0>"<<state->getpriorLengthsParam(0)<<"</lengthalpha0>"<<endl;
	*fileout<<"<lengthbeta0>"<<state->getpriorLengthsParam(1)<<"</lengthbeta0>"<<endl;	
	*fileout<<"<lengthdeltaalpha>"<<state->getpriorLengthsParam(2)<<"</lengthdeltaalpha>"<<endl;	
	*fileout<<"<lengthdeltabeta>"<<state->getpriorLengthsParam(3)<<"</lengthdeltabeta>"<<endl;		
	}
	if(state->useSums()){
	*fileout<<"<meanmualpha>"<<state->getpriorNumChunksParam(0)<<"</meanmualpha>"<<endl;
	*fileout<<"<meanmubeta>"<<state->getpriorNumChunksParam(1)<<"</meanmubeta>"<<endl;	
	*fileout<<"<meanmugamma>"<<state->getpriorNumChunksParam(2)<<"</meanmugamma>"<<endl;	
	}
	double ns=0,ni=0,nms=0,nh=0;
	if (numSAMS>0) ns=((double)accSAMS)/numSAMS;
	if (numIndiv>0) ni=((double)accIndiv)/numIndiv;
	if (numMergeSplit>0) nms=((double)accMergeSplit)/numMergeSplit;
	if (numHyper>0) nh=((double)accHyper)/numHyper;
	*fileout<<"<AccSAMS>"<<ns<<"</AccSAMS>"<<endl;
	*fileout<<"<AccMS>"<<nms<<"</AccMS>"<<endl;
	*fileout<<"<AccIndiv>"<<ni<<"</AccIndiv>"<<endl;
	*fileout<<"<AccHyper>"<<nh<<"</AccHyper>"<<endl;	*fileout<<"<Number>"<<iter<<"</Number>"<<endl;
	*fileout<<"</Iteration>"<<endl;
}

InfMCMC::~InfMCMC()
{
	if(state!=NULL) delete(state);
}


} // end namespace fines
