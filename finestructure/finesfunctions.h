#include <cstring>
#include <fstream>
#include "fines.h"
#include "data.h"
#include "prior.h"
#include "inf1.h"
#include "infextract.h"
#include "infextract2.h"
#include "infextract3.h"
#include "infextract4.h"
#include "infextract5.h"
#include "infmcmc.h"
#include "infadmixture.h"
#include "infconcordance.h"
#include "rng.h"
#include "fsxml.h"
#include "fines.h"

#define TREETYPE_USECONCORDANCESTATE 1
#define TREETYPE_USEMERGESTATE 2
#define TREETYPE_USEHILLCLIMBSTATE 3

using namespace std;
namespace fines
{

  std::vector<double> getBvec(int betamodel,int betamodelref,int datainference,vector<double> bvec=vector<double>(0),
			      string betapriorstring=string(""),long ignorelines=0,bool xhead=true,bool yhead=true);

 Inf1 mergeTree(my_rng * rng,int treetype, Data *d, string fs,long testmax,long hcsteps, Prior *prior,
	       int datainference=INFDATA_COUNTS, int modeltype=MODELTYPE_FINESTRUCTURE,State *startstate=NULL, Data *dlength=NULL, Data *dref=NULL, bool havefullxmlinput=true,bool fixK=false,int  treescale=0,string oname="",int maxconcordance=500,int verbose=0);

Inf1 GlobalReadTree(my_rng * rng,Data *d,string filename, Prior *prior,Data *dlength=NULL,Data *dref=NULL,
		    int datainference=INFDATA_COUNTS, int modeltype=MODELTYPE_FINESTRUCTURE,int verbose=0);
 
InfMCMC GlobalRunMCMC(my_rng * rng,State *initstate,ostream *os,long burnin,long additional,long thinin,string comment,
		      int datainference=INFDATA_COUNTS,double pcaprob=0,bool fixK=false,int verbose=0);

 void getXmlHeader(string filename, double &cval,double &cvalref,long &burnin, long &mcmclength, long &mcmcskip,string &datafilestr);

int compareDataFiles(string f1, string f2);
} // end namespace fines
