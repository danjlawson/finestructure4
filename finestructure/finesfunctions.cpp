#include "finesfunctions.h"

namespace fines
{
  /*  std::vector<double> addToBvec(vector<double> bvec,int betamodel){
    int bvecmin=1;
    if(bvec.size()==0 && betamodel==BETAMOD_CONST) bvec.push_back(1.0);
    if(betamodel==BETAMOD_F) {
      while(bvec.size()<4) bvec.push_back(-0.001); 
      bvecmin=4;
    }
    else if(betamodel==BETAMOD_F2 || betamodel==BETAMOD_F2_COPYMAT) {
      while(bvec.size()<2) bvec.push_back(2);
	while(bvec.size()<4) bvec.push_back(0.01);
	bvecmin=4;
    }
    if(betamodel==BETAMOD_COPYMAT || betamodel==BETAMOD_F2_COPYMAT) {
	Data *d2;
        try{d2=new Data(betapriorstring,ignorelines,xhead,yhead);
	//bvec.clear();//=vector<double>(d2->getN()*d2->getN(),0);
	for(int c1=0;c1<d2->getDim();c1++) {
	  for(int c2=0;c2<d2->getDim();c2++) bvec.push_back(d2->get(c1,c2,true));
	}
	}catch(string x){cerr<<"Error in betamodel creation:"<<x<<endl;exit(0);}
	bvecmin+=d2->getDim()*d2->getDim();
	delete(d2);	
    }
    return(bvec);
    
  }
  
  std::vector<double> getBvec(int betamodel,int refbetamodel,int datainference,vector<double> bvec,
			      string betapriorstring,long ignorelines,bool xhead,bool yhead) {
    // Add the counts square matrix parameters
    addToBvec(bvec,betamodel, betapriorstring, ignorelines, xhead, yhead);
    // Add the counts reference matrix parameters
    while(bvec.size()<2) bvec.push_back(2);
	while(bvec.size()<4) bvec.push_back(0.01);

    //    addToBvec(bvec,refbetamodel, betapriorstring, ignorelines, xhead, yhead);    
    // Add any other parameters
    int bvecmin=bvec.size();
    int bvecmin2=0;
    if(datainference==INFDATA_LENGTHS || datainference==INFDATA_ALL){
	bvecmin2=NUMHYPERPARAMLENGTH*2;
	while((long)bvec.size()<bvecmin+NUMHYPERPARAMLENGTH) bvec.push_back(2);
	while((long)bvec.size()<bvecmin+NUMHYPERPARAMLENGTH) bvec.push_back(0.01);
    }
    if(datainference==INFDATA_TOTALLENGTHS || datainference==INFDATA_ALLNOTLENGTHS || datainference==INFDATA_ALL){
	while((long)bvec.size()<bvecmin+bvecmin2+NUMHYPERPARAMTOTLENGTH) bvec.push_back(2);
	while((long)bvec.size()<bvecmin+bvecmin2+2*NUMHYPERPARAMTOTLENGTH) bvec.push_back(0.01);
    }
    return(bvec);
}*/

  Inf1 mergeTree(my_rng * rng,int treetype, Data *d, string fs,long testmax,long hcsteps, Prior *prior,
	       int datainference,int modeltype, State *startstate, Data *dlength, Data *dref, bool havefullxmlinput,bool fixK,int treescale,string oname,int maxconcordance,int verbose) {
	if(verbose==1) cout<<"MERGE PHASE"<<endl;

	FsXml *infile=new FsXml(fs);
//	InfExtract iext(d,infile,verbose);
//	delete(infile);
//	infile=new FsXml(fs);
	Inf1 * inf1_i=NULL;
	InfExtract2 * iext2=NULL;
	InfMCMC * infHillClimb=NULL;
	State * state2;
	if(havefullxmlinput){
	  //	  infile->gotoLineContaining("<Number>9990</Number>");///*******************************
	  iext2=new InfExtract2(rng,d,infile,prior,verbose,dlength,dref,datainference,modeltype);
	  if(treetype==TREETYPE_USECONCORDANCESTATE){// ***
	    state2=new State(iext2->getState()); // this is the original state, stored for later
	    InfConcordance infc(rng, state2,maxconcordance,d,infile,prior,verbose,dlength,dref,datainference,modeltype);

	    /////////////////////////
	    filebuf fb;
	    oname.append("_concordance.csv");
	    try{ fb.open (oname.c_str(),ios::out);
	    }catch(std::string x){
	      cerr<<"Error opening file!"<<endl<<x<<endl; exit(1);}
	    ostream os (&fb);
	    infc.printScoreMatrix(&os);
	    fb.close();
	    /////////////////////////
	    
	    inf1_i=new Inf1(rng,state2,datainference,verbose,testmax,treescale);
	  }else if(treetype==TREETYPE_USEMERGESTATE) {
	    inf1_i=new Inf1(rng,iext2->getState(),datainference,verbose,testmax,treescale);
	    try{inf1_i->mergeHillClimb(NULL,true,false);}catch(std::string x){cout<<x<<endl;exit(0);}
	    state2=new State(inf1_i->getState());
	    
		}else if(treetype==TREETYPE_USEHILLCLIMBSTATE) {
	    try{infHillClimb=new InfMCMC(rng,iext2->getState(),datainference,0,verbose);
			//infHillClimb->hillClimb(0,opt().burnin);
			if(fixK)  infHillClimb->fixK();
			infHillClimb->metropolis(0,hcsteps);
			state2=new State(infHillClimb->getState());
			}catch(std::string x){cout<<x<<endl;exit(0);}
		}else {cerr<<"Invalid tree type (-T option): "<<treetype<<endl;exit(0);}
		delete(infile);
	}else {
	  if(treetype==TREETYPE_USECONCORDANCESTATE){// ***
	    // THIS IS NOT OK AND SHOULD GIVE AN ERROR
	    //inf1_i=new Inf1(rng,d,startstate,dlength,datainference,verbose,testmax,treescale);
	    // state2=new State(startstate);
	    cerr<<"Invalid tree type (concordance) when not providing a valid mcmc xml file."<<endl;exit(0);
		}else if(treetype==TREETYPE_USEMERGESTATE) {
	    inf1_i=new Inf1(rng,startstate,datainference,verbose,testmax,treescale);
			try{inf1_i->mergeHillClimb(NULL,true,false);}catch(std::string x){cout<<x<<endl;exit(0);}
			state2=new State(inf1_i->getState());

		}else if(treetype==TREETYPE_USEHILLCLIMBSTATE){
	    try{infHillClimb=new InfMCMC(rng,startstate,datainference,0,verbose);
			infHillClimb->hillClimb(0,hcsteps);
			state2=new State(infHillClimb->getState());
			}catch(std::string x){cout<<x<<endl;exit(0);}
		}else {cerr<<"Invalid tree type (-t option)."<<endl;exit(0);}
	};
	//	state2->iterPrint(&cout);
	Inf1 inf1(rng,state2,datainference,verbose,testmax,treescale);
	//	inf1.getState()->iterPrint(&cout);
	return(inf1);
/*	inf1.exportXmlHead(&os);
	inf1.exportXmlComment(&os,comment);
	try{inf1.mergeHillClimb(NULL,false,treesuper);}catch(std::string x){cout<<x<<endl;exit(0);}
	if(opt().verbose) cout<<"Assigning certainty"<<endl;
	infile=new FsXml(fs);
	InfExtract3 iext3(d,infile,inf1.getNodes(),opt().verbose);
	delete(infile);
	if(opt().verbose) cout<<"Diagonalise tree"<<endl;
	inf1.diagonaliseOnM(d->getMatrix(),false);
*/
}

  Inf1 GlobalReadTree(my_rng * rng,Data *d,string filename, Prior *prior,Data *dlength,Data *dref,
		      int datainference, int modeltype,int verbose)
{
  FsXml infile(filename);
  
  //  Data *d2=NULL;
  string pop=infile.getParam("Pop");
  string newick=infile.getParam("Tree");
//  cout<<"READING TREE:"<<newick<<endl;
  State *state = new State(rng,d,pop,prior,false,dlength,dref,datainference,modeltype);//*** dref
  Inf1 ret=Inf1(rng,state,newick,verbose);
  return(ret);
}

InfMCMC GlobalRunMCMC(my_rng * rng,State *initstate,ostream *os,long burnin,long additional,long thinin,string comment,
		      int datainference,double pcaprob,bool fixK,int verbose)
{
  //  cout<<"TEST1"<<endl;
  InfMCMC infMCMC(rng,initstate,datainference,pcaprob,verbose);
  try{
    //  cout<<"TEST2"<<endl;
    if(fixK) infMCMC.fixK();
    //   cout<<"TEST2a"<<endl;
    infMCMC.exportXmlHead(os,burnin,additional,thinin);
    //  cout<<"TEST2b"<<endl;
	infMCMC.exportXmlComment(os,comment);
	if(verbose) cout<<"BURN IN PHASE"<<endl;
	infMCMC.metropolis(0,burnin);
	if(verbose) cout<<"MCMC PHASE"<<endl;
	infMCMC.resetCounters();	
	for(long c1=0;c1<additional;c1++){infMCMC.metropolis(c1,1,thinin,os,additional);}
	if(additional % thinin==0) infMCMC.exportXmlIter(os,additional);

  }catch(string x){
    cerr<<"Error in GlobalRunMCMC:"<<x<<endl;
    throw(x);
  }	
	//infMCMC.metropolis(0,additional,thinin,os);
  return(infMCMC);
}

  bool getXmlTagHeader(FsXml *infile, string tag,string &val){
    infile->rewind();
    infile->gotoLineContaining("<header>");
    try{
      val=infile->getParam(tag);
      return(true);
    }catch(string x){
      val="";
      return(false);
    }
    
  }

  void getXmlHeader(string filename, double &cval,double &cvalref,long &burnin, long &mcmclength,long &mcmcskip,string &datafilestr){
    FsXml *infile=new FsXml(filename);
    string val;
    if(getXmlTagHeader(infile,"inflation",val)){
      cval=atof(val.c_str());
    }
    if(getXmlTagHeader(infile,"inflationRef",val)){
      cvalref=atof(val.c_str());
    }
    if(getXmlTagHeader(infile,"burnin",val)){
      burnin=atoi(val.c_str());
    }
    if(getXmlTagHeader(infile,"mcmclength",val)){
      mcmclength=atoi(val.c_str());
    }
    if(getXmlTagHeader(infile,"mcmcskip",val)){
      mcmcskip=atoi(val.c_str());
    }
    if(getXmlTagHeader(infile,"datafilestr",val)){
      datafilestr=val;
    }
    delete(infile);
}

int compareDataFiles(string f1, string f2) {
  	if(f1.size()==0 || f2.size()==0){
//	   cerr<<"WARNING! Cannot confirm data file is the same as the MCMC was run on!"<<endl;
	  return(1);
	}else if(false){
// see if they have the same file name but different directories	  
	  return(2);
	}else if(f1!=f2){
//	  cerr<<"WARNING!  You are trying to build a tree from a differently named datafile than the one used for the MCMC!  This might be due to running it on a different system or might imply that the file is incorrect. This may result in strange behaviour!"<<endl;
	  return(-1);
	}
    return(0);
}

}
