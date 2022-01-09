
#include "prior.h"
#include "data.h"

namespace fines
{ 
  Prior::Prior(){
    alpha=1.0;
    setPriorCounts(BETAMOD_F2);
    setPriorRef(BETAMOD_F2);
    betamodellength=BETAMOD_EQUI;
    betamodeltotlength=BETAMOD_EQUI;    
  }

  void Prior::setPriorCounts(std::vector<double> bvec){
    if(betamodel==BETAMOD_EQUI) {
      if(bvec.size()!=0){
	cerr<<"Beta model with constant parameter requires no parameter but you provided "<<bvec.size()<<endl;
	throw(string("Invalid beta parameter"));
      }
    }else if(betamodel==BETAMOD_CONST) {
      if(bvec.size()!=1){
	cerr<<"Beta model with constant parameter requires 1 parameter but you provided "<<bvec.size()<<endl;
	throw(string("Invalid beta parameter"));
      }
    }else if(betamodel==BETAMOD_F) {
      if(bvec.size()!=4){
	cerr<<"Beta model F requires 4 parameters but you provided "<<bvec.size()<<endl;
	throw(string("Invalid beta parameter"));
      }
    }else if(betamodel==BETAMOD_F2 || betamodel==BETAMOD_F2_COPYMAT) {
      if(bvec.size()!=4){
	cerr<<"Beta model F2 requires 4 parameters but you provided "<<bvec.size()<<endl;
	throw(string("Invalid beta parameter"));
      }
    }else{
      cerr<<"Unimplemented Beta model: "<<betamodel<<endl;
      throw(string("Invalid beta model"));
    }
    this->bvec=bvec;
  }

  void Prior::setPriorCounts(std::string tmp){
    bvec.clear();
    if(betamodel==BETAMOD_COPYMAT || betamodel==BETAMOD_F2_COPYMAT) { // tmp should contain a file name of the beta matrix
      cerr<<"TEMPORARILY DISABLED THE COPY MATRIX PRIOR"<<endl;
      throw(string("Prior invalid"));
      /*      Data *d2;
        try{d2=new Data(betapriorstring,ignorelines,xhead,yhead);
	//bvec.clear();//=vector<double>(d2->getN()*d2->getN(),0);
	for(int c1=0;c1<d2->getDim();c1++) {
	  for(int c2=0;c2<d2->getDim();c2++) bvec.push_back(d2->get(c1,c2,true));
	}
	}catch(string x){cerr<<"Error in betamodel creation:"<<x<<endl;exit(0);}
	delete(d2);	*/
    }else{ // tmp should contain a vector of the hyper parameters
      vector<double> tbvec;
      splitasdouble(tmp,',',tbvec);
      setPriorCounts(tbvec);
    }
  }
  
  void Prior::setPriorCounts(int betamodel){
    this->betamodel=betamodel;
    bvec.clear();

    // Default values
    if(betamodel==BETAMOD_EQUI) {
      // We're good
    }else if(betamodel==BETAMOD_CONST) {
      bvec.push_back(1.0);
    }else if(betamodel==BETAMOD_F) {
      while(bvec.size()<4) bvec.push_back(-0.001); 
    }else if(betamodel==BETAMOD_F2 || betamodel==BETAMOD_F2_COPYMAT) {
      while(bvec.size()<2) bvec.push_back(2);
	while(bvec.size()<4) bvec.push_back(0.01);
    }else{
      cerr<<"Unimplemented Beta model: "<<betamodel<<endl;
      throw(string("Invalid beta model"));
    }   
  }

  ////////////////// REFERENCE MODEL

  void Prior::setPriorRef(std::vector<double> bvec){
    this->bvecref=bvec;
    if(betamodelref==BETAMOD_EQUI) {
      if(bvecref.size()!=0){
	cerr<<"Beta Ref model with constant parameter requires no parameter but you provided "<<bvecref.size()<<endl;
	throw(string("Invalid beta parameter"));
      }
      // We're good
    }else if(betamodelref==BETAMOD_CONST) {
      if(bvecref.size()!=1){
	cerr<<"Beta Ref model with constant parameter requires 1 parameter but you provided "<<bvecref.size()<<endl;
	throw(string("Invalid beta parameter"));
      }
    }else if(betamodelref==BETAMOD_F2) {
      if(bvecref.size()!=2){
	cerr<<"Beta Ref model F2 requires 2 parameters but you provided "<<bvecref.size()<<endl;
	throw(string("Invalid beta parameter"));
      }
    }else{
      cerr<<"Unimplemented Beta model: "<<betamodelref<<endl;
      throw(string("Invalid beta model"));
    }

  }

  void Prior::setPriorRef(std::string tmp){
    bvec.clear();
    vector<double> tbvec;
    splitasdouble(tmp,',',tbvec);
    setPriorCounts(tbvec);
  }

  
  void Prior::setPriorRef(int betamodel){
    this->betamodelref=betamodel;
    bvecref.clear();
    // Default values
    if(betamodelref==BETAMOD_EQUI) {
      // We're good
    }else if(betamodelref==BETAMOD_CONST) {
      bvecref.push_back(1.0);
    }else if(betamodelref==BETAMOD_F2 ) {
      bvecref.push_back(2);
      bvecref.push_back(0.01);
    }else{
      cerr<<"Unimplemented Beta model: "<<betamodelref<<" for Reference data"<<endl;
      throw(string("Invalid beta model"));
    }   
  }

  ////////////////// LENGTHS MODEL
  
  void Prior::setPriorLength(int betamodel){
    this->betamodellength=betamodel;
    bveclength.clear();
    // Default values
    if(betamodellength==BETAMOD_CONST) {
      bveclength.push_back(1.0);
    }else if(betamodellength==BETAMOD_F2 ) {
      bveclength.push_back(2);
      bveclength.push_back(0.01);
    }else{
      cerr<<"Unimplemented Beta model: "<<betamodellength<<" for Lengths data"<<endl;
      throw(string("Invalid beta model"));
    }   
  }
    
    ////////////////// TOTLENGTHS MODEL
  
  void Prior::setPriorTotLength(int betamodel){
    this->betamodeltotlength=betamodel;
    bvectotlength.clear();
    // Default values
    if(betamodeltotlength==BETAMOD_CONST) {
      bvectotlength.push_back(1.0);
    }else if(betamodeltotlength==BETAMOD_F2 ) {
      bvectotlength.push_back(2);
      bvectotlength.push_back(0.01);
    }else{
      cerr<<"Unimplemented Beta model: "<<betamodeltotlength<<" for Total lengths data"<<endl;
      throw(string("Invalid beta model"));
    }   
  }

  Prior::~Prior(){
  }
}
