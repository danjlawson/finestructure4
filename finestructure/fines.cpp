#include "fines.h"
#include "data.h"
#include "prior.h"
#include "inf1.h"
#include "infextract.h"
#include "infextract2.h"
#include "infextract3.h"
#include "infextract4.h"
#include "infextract5.h"
#include "infextractdonors.h"
#include "infmcmc.h"
#include "infadmixture.h"
#include "rng.h"
#include "fsxml.h"
#include "finesfunctions.h"

using namespace std;
namespace fines
{


// The f model proper has been removed:
//			3/f/F:  F model of Falush et al 2003 (RJ MOVE NOT IMPLEMENTED).
//			(default: -b -0.001,-0.001,-1,-1 for model 3,.
  
static const char * help=
    "\
    Usage: finestructure [OPTIONS] datafile <initialpopfile> outputfile\n\
    	Datafile is a matrix of copy counts.\n\
	initialpopfile (optional) is a population state e.g. an outputfile.\n\
	outputfile is the destination.\n\
IMPORTANT OPTIONS FOR BASIC USAGE:\n\
    -m <method>		Method to use.  Default: oMCMC.\n\
			<method> is either oMCMC (MCMC without tree), \n\
			Tree, or a contraction of either.\n\
    -x <num>		Number of burn in iterations for MCMC method (default: 1000).\n\
    -y <num>		Number of sample iterations for MCMC method (default: 1000).\n\
    -z <num>		Thin interval in the output file, for MCMC method (default: 1).\n\
    -t <num>		Maximum number of tree comparisons for splitting/merging (default: 1500).\n\
    -v          	Verbose mode\n\
    -V          	Print Version info\n\
    -h          	This help message\n\
    -H          	More detailed help with all options listed\n\
IMPORTANT OPTIONS FOR TREE BUILDING:\n\
    -T <type>		When using a merge tree, initialisation can be set to the following:\n\
			1:	Use the \"Maximum concordance State\" as used by Leslie et al 2015 (PoBI paper, Nature).\n\
			3:	Perform full range of moves to to get to best posterior state.\n\
				This is the default.  Set number of attempts with -x <num>, disable\n\
                                any change to the state with -x 0.\n\
    -k <num>		Change the tree building algorithm.\n\
			0:	Discard all ordering and likelihood information (default).\n\
			2:	Maintain ordering and likelihood.\n\
IMPORTANT EXTRACTION FROM THE MCMC OUTPUT:\n\
    -e <name>		Extract details from a state; can be (a unique contraction of):\n\
			X2: the normalised copying matrix (specify <chunkcounts.out> <optional mcmc.xml> <output>)\n\
			X: the copying data matrix for populations (specify <chunkcounts.out> <mcmc.xml> <output>)\n\
			meancoincidence: the mean coincidence matrix (specify <chunkcounts.out> <mcmc.xml> <output>)\n\
			maxstate: maximum observed posterior probability state (specify <chunkcounts.out> <mcmc.xml> <output>)\n\
			popidfile: extract a ChromoPainterv2 idfile containing population information from the best population state. This consists of names\\tpopulation tags\\t1 (the final 1 is for inclusion) (specify <chunkcounts.out> <tree.xml> <output>)\n\
  			range:<from>:<to> gets the iterations in the specified range.   (specify <chunkcount.out> <optional mcmc.xml> <output.xml>)\n\
			thin:<step>: thins the output by step.  (specify <chunkcount.out> <optional mcmc.xml> <output.xml>)\n\
FIXING SOME POPULATIONS:\n\
    -F <filename>	Fix the populations specified in the file.  They should be specified as\n\
			population format, i.e. PopA(ind1,ind2) would fix the data rows ind1 and ind2\n\
			to always be in the same population (they form a 'super individual')\n\
			called PopA. Continents are specified with a * before the name, and are treated\n\
			specially in the tree building phase,  i.e. *ContA(ind1,ind2).  Continents\n\
			are not merged with the rest of the tree.\n\
	\n\
    Examples:\n\
    finestructure -x 100000 -y 100000 -z 1000 datafile.csv out.mcmc.xml\n\
			Infers population structure using MCMC from datafile.csv. \n\
                        100000 burn in steps are used (-x)\n\
			and 100000 further iterations are sampled (-y) keeping every\n\
			1000th sample (-z).\n\
    finestructure -x 0 -y 100000 -z 1000 datafile.csv out.mcmc.xml out.mcmc.longer.xml\n\
                        continue the previous run for 100K additional steps, treating the original run as burnin.\n\
    finestructure -m T -x 0 datafile.csv out.mcmc.xml out.tree.xml\n\
			Infers a tree, using the best state seen in out.mcmc.xml as the initial state.\n\
    finestructure -m T -k 2 -T 1 datafile.csv out.mcmc.xml out.tree.xml\n\
			Infers a tree, using (-T 1) the maximum concordance state over out.mcmc.xml as the initial state. This is reported with full likelihood ordering (-k 2), useful for cutting at a given number of ppulations K (but may look bad in the GUI).\n\
    ";


static const char * fullhelp=
  "\
ADVANCED USAGE:\n\
SPECIFYING PARAMETERS:\n\
    -a <num>		Set alpha, the prior of the number of parameters\n\
			(default: 1.0).\n\
    -c <num>		Set the likelihood correction factor: L_{used}=L^{1/<corfactor>}.\n\
			(default: 1.0)\n\
NON-SQUARE MATRICES:\n\
    -R                  Specify that the data file is painted against a reference dataset, and therefore\n\
                        is not square. columns are not clustered, and currenty only an equipartition prior can be used. \n\
    -r <filename>	Specify a reference chunk lengths datafile in addition to the main datafile.\n\
                        Any -i,-X,-Y,-c,-B,-b options\n\
			*preciding* this file will affect counts; those following affect refernce.\n\
FIXING THE NUMBER OF POPULATIONS:\n\
    -K                  Fix the number of populations to whatever you started with.\n\
                        This would be set by '-I' or by an initial state file.\n\
INITIAL VALUES:\n\
    -I <x>          	Initial number of populations.  <x> is either a number\n\
			or \"n\" for the number of individuals, or \"l\" for label detected \n\
			populations.  Default is 1.\n\
    -s <s>		Sets the RNG seed to s (>0)\n\
    -i <i>		Ignores the first i lines of the input file\n\
TREE OPTIONS:\n\
    -k <num>		Change the tree building algorithm.\n\
			0:	Discard all ordering and likelihood information (default).\n\
			1:	Maintain ordering.\n\
			2:	Maintain ordering and likelihood.\n\
    -T <type>		When using a merge tree, initialisation can be set to the following:\n\
			1:	Use the \"Maximum concordance State\" as used by Leslie et al 2015 (PoBI paper, Nature).\n\
			2:	Perform merging to get to best posterior state.\n\
			3:	Perform full range of moves to to get to best posterior state.\n\
				This is the default.  Set number of attempts with -x <num>.\n\
			4:	As 1, but don't flatten maximum copy rates for the main tree.\n\
			5:	As 2, but don't flatten maximum copy rates for the main tree.\n\
			6:	As 3, but don't flatten maximum copy rates for the main tree.\n\
			7:	As 1, but maximise hyperparameters between merges.\n\
			8:	As 2, but maximise hyperparameters between merges.\n\
			9:	As 3, but maximise hyperparameters between merges.\n\
    -m <method>		Method to use.  Default: oMCMC.\n\
			<method> is either MCMCwithTree, oMCMC (MCMC without tree), \n\
			Tree, or a contraction of any.\n\
    -O <name>		File containing a state to use for ordering, if not the main file.\n\
    -C <val>            Maximum number of concordance iterations (default: 500)\n\
CHANGING THE MODEL:\n\
    -l <filename>	Specify the average copy length datafile.  -i,-X,-Y options\n\
			*preciding* this file will affect this read; you can set different\n\
			options for the copy rate datafile by specifying these -i,-X,-Y again\n\
			after the -l option.\n\
    -B <model>		Choose a model for beta:\n\
			1/e/E:	Equipartition model of Pella and Masuda.\n\
			2/c/C:	Constant model.\n\
			4/o/O:  F model of Falush et al 2003 with a single parameter\n\
				for all populations (default).\n\
    -b <num>(,<num>,..)	Hyperparameters for the current model.  \n\
                        COUNTS:\n\
			For model 1, there are no parameters (equipartition model is Dir(beta) with beta=1/K).\n\
			For model 2, one parameter, Dir(beta) with all beta_i=<num>).\n\
			(default: 1.0).\n\
			For model 4, set the hyperprior of the distribution of\n\
			delta and F. Parameters are \n\
			(k_f,k_delta,theta_f,theta_delta) for the parameters of the\n\
			gamma distribution F~Gamma(k_f,theta_f), \n\
			and delta~Gamma(k_delta,theta_delta)\n\
			(default: -b 2,2,0.01,0.01).\n\
			REFERENCE COUNTS: As COUNTS. Exception: the F model (4) has only 2 parameters (k_f,theta_f):\n\
			gamma distribution F~Gamma(k_f,theta_f), \n\
                        as there can be no increased copying within a population.\n\
FULL SET OF EXTRACTION OPTIONS:\n\
    -e <name>		Extract details from a state; can be (a unique contraction of):\n\
			X2: the normalised copying matrix (specify <chunkcounts.out> <optional mcmc.xml> <output>)\n\
			X: the copying data matrix for populations (specify <chunkcounts.out> <mcmc.xml> <output>)\n\
			meancoincidence: the mean coincidence matrix (specify <chunkcounts.out> <mcmc.xml> <output>)\n\
			maxstate: maximum observed posterior probability state (specify <chunkcounts.out> <mcmc.xml> <output>)\n\
			popidfile: extract a ChromoPainterv2 idfile containing population information from the best population state. This consists of names\\tpopulation tags\\t1 (the final 1 is for inclusion) (specify <chunkcounts.out> <tree.xml> <output>)\n\
			beta: the parameter matrix (specify <chunkcount.out> <optional mcmc.xml> <output>)\n\
			merge<:value><:split>: create a merge(or split) \n\
			  population from the mean coincidence.  (specify <chunkcount.out> <optional mcmc.xml> <output>)\n\
  			range:<from>:<to> gets the iterations in the specified range.   (specify <chunkcount.out> <optional mcmc.xml> <output.xml>)\n\
			thin:<step>: thins the output by step.  (specify <chunkcount.out> <optional mcmc.xml> <output.xml>)\n\
			probability: get the posterior probability of the data\n\
			given the conditions of the outputfile.\n\
			likelihood: samples the likelihood of the data given the conditions\n\
			in the outputfile.\n\
			tree: extract the tree in newick format and print it to a FOURTH file  (specify <chunkcounts.out> <mcmc.xml> <tree.xml> <output>)\n\
FINE CONTROL OVER THE DATA:\n\
    -X			Specifies that there are row names in the data (not necessary for \n\
			ChromoPainter or ChromoCombine style files.)\n\
    -Y			Specifies that there are column names in the data file (as -X, not necessary.)\n\
";
  
static const char * undochelp=
  "\
   Undocumented usage: These features are EXPERIMENTAL and UNEXPLORED in papers. They may not work\n\
    -M <modeltype>	Specify the type of inference model for chunk counts.  \n\
			<modeltype> accept contractions and lower case, and can be:\n\
			  1 or Finestructure: standard finestructure model (default).\n\
			  2 or Normalised: Normalise data row and columns within a population.\n\
			  3 or MergeOnly: As 2, but only compare populations being merged or split.\n\
			  4 or Individual: Prior is placed on individual rows instead of \n\
					  population rows. (slowest model).\n\
                        LENGTHS: 8 parameters:\n\
			(k_alpha0,k_beta0,k_alpha,k_beta,beta_alpha0,beta_beta0,beta_alpha,beta_beta)\n\
			MEANS: 6 parameters:\n\
			(k_betamu, k_alphamu, k_kappa, beta_alphamu,beta_betamu,beta_kappa)\n\
			Set K parameters negative for fixed =|k|\n\
			e.g. when finding a tree given the mean parameters.\n\
    -u <datatype>	Use a data inference method; one of :\n\
			counts: use only the copy counts data. (default if -l not specified)\n\
			lengths: use only the copy length data (still needs valid counts data!)\n\
			totallengths: use the mean length of chunk sizes \n\
			all: use all data (careful: this may not be statistically valid).\n\
			default: use counts and totallengths (default with -l specified).\n\
";

/*    -p <num>		Use the PCA enhanced merge/split proposals.  <num> is between zero and one\n\
			and controls the probability with which PCA proposals are used.\n\
			PCA proposals are much faster but may increase autocorrelation time.\n\
			In most cases, the increased speed is worth it, but it may cause problems\n\
			locating the posterior mode if there are a lot of individuals.  Default: 0.\n\*/
// ********************** PCA MOVE CURRENTLY DOES NOTHING
    
std::vector<double> &splitasdouble(const std::string &s, char delim, std::vector<double> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
    	elems.push_back(atof(item.c_str()));
    }
    return elems;
}
    
string getVersion(){
	string ret;
	//#ifdef PACKAGE_STRING // Omitted due to version matching issues
	//ret.append(PACKAGE_STRING);
	//#else  
	ret.append("finestructure");
	//#endif
	ret.append(" build date "); 
	ret.append(__DATE__);
	ret.append(" at ");
	ret.append(__TIME__);
	return(ret);
}

void printVersion(){
	cout<<getVersion()<<endl;
}


  State setState(my_rng * rng,string fs2,Data * d,Prior *prior,Data *dlength=NULL,Data *dref=NULL,int datainference=INFDATA_COUNTS,int modeltype=MODELTYPE_FINESTRUCTURE)
{
  //	if(opt.verbose) cout<<"Reordering on "<<fs2<<"."<<endl;
	FsXml infile(fs2);
	InfExtract2 iext2r(rng,d,&infile,prior,false,dlength,dref,datainference,modeltype);
	return(iext2r.getState());
}

int getModelType(string modelarg){
  std::transform(modelarg.begin(), modelarg.end(),modelarg.begin(), ::toupper);
  if(modelarg.substr(0,1).compare("F")==0 || modelarg.substr(0,1).compare("1")==0) return(MODELTYPE_FINESTRUCTURE);
  if(modelarg.substr(0,1).compare("N")==0 || modelarg.substr(0,1).compare("2")==0) return(MODELTYPE_NORMALISED);
  if(modelarg.substr(0,1).compare("M")==0 || modelarg.substr(0,1).compare("3")==0) return(MODELTYPE_NORMALISEDMERGEONLY);
  if(modelarg.substr(0,1).compare("I")==0 || modelarg.substr(0,1).compare("4")==0) return(MODELTYPE_INDIVIDUAL);
  return(1);
}

  int getBetaModel(string tmp, string &betapriorstring){
    int betamodel=0;
    char *tmpchar;
     char * tmpc = new char[tmp.size() + 1];
    std::copy(tmp.begin(), tmp.end(), tmpc);
    tmpc[tmp.size()] = '\0';  // Ensure null termination
    if(tmp.compare("1")==0 || tmp.substr(0,1).compare("e")==0 || tmp.substr(0,1).compare("E")==0) {betamodel=BETAMOD_EQUI;
    }else if(tmp.compare("2")==0 || tmp.substr(0,1).compare("c")==0 || tmp.substr(0,1).compare("C")==0) {betamodel=BETAMOD_CONST;
    }else if(tmp.compare("3")==0 || tmp.substr(0,1).compare("f")==0 || tmp.substr(0,1).compare("F")==0) {betamodel=BETAMOD_F;
    }else if(tmp.compare("4")==0 || tmp.substr(0,1).compare("o")==0 || tmp.substr(0,1).compare("O")==0) {betamodel=BETAMOD_F2;
    }else if(tmp.substr(0,1).compare("5")==0 || tmp.substr(0,1).compare("x")==0 || tmp.substr(0,1).compare("X")==0) {betamodel=BETAMOD_COPYMAT;tmpchar=strtok(tmpc,":"); 
      tmpchar=strtok(NULL,":");if(tmpchar==NULL){ cerr<<"Must specify a datafile for beta matrix with this prior!"<<endl;exit(0);}; betapriorstring=string(tmpchar);
    }else if(tmp.substr(0,1).compare("6")==0 || tmp.substr(0,1).compare("p")==0 || tmp.substr(0,1).compare("P")==0){
	betamodel=BETAMOD_F2_COPYMAT;tmpchar=strtok(tmpc,":"); 
      tmpchar=strtok(NULL,":");if(tmpchar==NULL){ cerr<<"Must specify a datafile for beta matrix with this prior!"<<endl;exit(0);}; betapriorstring=string(tmpchar);
    }
    return(betamodel);
  }

int finestructure(int argc, char *argv[])
{
  ProgramOptions opt;
    string comment="Command line: ";
    for(int c1=0;c1< argc;c1++) {comment.append(argv[c1]);comment.append(" ");}
    comment.append("\nVersion: ");
    comment.append(getVersion());
    std::stringstream ss;
    std::string tmp;
    string fs;// The input file
    string fs2;// The input file for ordering
    string fstree;// The *input* file for the tree (only for -e TREE)
    string betapriorstring;// the input file for the beta prior, if betamodel=5
    Prior *prior=new Prior();
    Data *dsafe=NULL; // The data set that is guaranteed to have been read in which may not be a chunkcounts file)
    Data *dcount=NULL, *dlength=NULL,*dref=NULL; // data for passing through to the state creation for counts, lengths and reference
    double corfactor=-1.0;
    bool setcorfactor=false;// Whether we have manually set the corfactor, to override the value in an input xml file
    bool setcorfactorref=false;
    int initpop=-1;
    int maxconcordance=500;
    int ignorelines=0;
    bool xhead=false,yhead=false;
    unsigned long seed=0;
    bool extract=false;string ext;
    bool havefullxmlinput=false;
    int datainference=INFDATA_ALLNOTLENGTHS;
    int treetype=TREETYPE_USEHILLCLIMBSTATE;
    int modeltype=MODELTYPE_FINESTRUCTURE;
    int treemodification=TREEMOD_FLATTEN;
    bool initpopset=false;
    string fixfile;
    string priorstring; // The string provided to us about the prior, which is retained since we may not know whether to apply it to the ref or copy matrix
    std::vector<std::string> args(argv, argv+argc) ;
    args.erase(args.begin());

    unsigned int argon=0;
    //(c = getopt (argc, argv, "x:y:z:s:l:r:u:I:i:m:p:T:XYa:b:M:F:B:e:t:o:O:c:vKk:hVS")) != -1
    while (argon<args.size())  {
  if(args[argon].compare("-x")==0 &&(argon+1<args.size())) {
    opt.burnin=atoi(args[argon+1].c_str());      
    args.erase(args.begin()+argon,args.begin()+argon+2);
  }else if(args[argon].compare("-y")==0 &&(argon+1<args.size())) {
    opt.additional=atoi(args[argon+1].c_str());
    args.erase(args.begin()+argon,args.begin()+argon+2);
  }else if(args[argon].compare("-z")==0 &&(argon+1<args.size())) {
    opt.thinin=atoi(args[argon+1].c_str());
    args.erase(args.begin()+argon,args.begin()+argon+2);
  }else if(args[argon].compare("-s")==0 &&(argon+1<args.size())) {
    seed=strtoul(args[argon+1].c_str(),NULL,10);
    args.erase(args.begin()+argon,args.begin()+argon+2);
  }else if(args[argon].compare("-l")==0 &&(argon+1<args.size())) {
    dlength=new Data(string(args[argon+1]),ignorelines,xhead,yhead,1-opt.userefonly);
    args.erase(args.begin()+argon,args.begin()+argon+2);
  }else if(args[argon].compare("-R")==0) {
    opt.userefonly=1;
    if(dref!=NULL) {cerr<<"ERROR: Cannot specify both -r and -R."<<endl; exit(1);}
    prior->setPriorRef(prior->getBetaModel());
    args.erase(args.begin()+argon);
  }else if(args[argon].compare("-r")==0 &&(argon+1<args.size())) {
    if(opt.userefonly) {cerr<<"ERROR: Cannot specify both -r and -R."<<endl; exit(1);}
    dref=new Data(string(args[argon+1]),ignorelines,xhead,yhead,false);
    if(priorstring.compare("")!=0){ // We've already had a previous prior string set
      // We need to apply it
      prior->setPriorCounts(priorstring);
      priorstring="";
    }
    if(corfactor>=0) {
      setcorfactorref=true;
      dref->setCfactor(corfactor);
      corfactor=-1.0;}
    args.erase(args.begin()+argon,args.begin()+argon+2);
  }else if(args[argon].compare("-u")==0 &&(argon+1<args.size())) {
    opt.usedata=string(args[argon+1]);
    args.erase(args.begin()+argon,args.begin()+argon+2);
  }else if(args[argon].compare("-I")==0 &&(argon+1<args.size())) {
    if(args[argon+1][0]!='n' && args[argon+1][0]!='N') {
      if(args[argon+1][0] || args[argon+1][0]=='L') {initpop=-2;
      }else initpop=atoi(args[argon+1].c_str());
    }else initpop=-1;
    initpopset=true;
    args.erase(args.begin()+argon,args.begin()+argon+2);
  }else if(args[argon].compare("-i")==0 &&(argon+1<args.size())) {
    ignorelines=atoi(args[argon+1].c_str()); break;
    args.erase(args.begin()+argon,args.begin()+argon+2);
  }else if(args[argon].compare("-e")==0 &&(argon+1<args.size())) {
    extract=true; 
    if(!initpopset) initpop=-1;
    ext=string(args[argon+1]);opt.method.assign("");
    args.erase(args.begin()+argon,args.begin()+argon+2);
  }else if(args[argon].compare("-O")==0 &&(argon+1<args.size())) {
    fs2=string(args[argon+1]);
    args.erase(args.begin()+argon,args.begin()+argon+2);
  }else if(args[argon].compare("-p")==0 &&(argon+1<args.size())) {
    opt.pcaprob=atof(args[argon+1].c_str());
    args.erase(args.begin()+argon,args.begin()+argon+2);
  }else if(args[argon].compare("-m")==0 &&(argon+1<args.size())) {
    opt.method.assign(args[argon+1]);	
    std::transform(opt.method.begin(), opt.method.end(),opt.method.begin(), ::toupper);
    args.erase(args.begin()+argon,args.begin()+argon+2);
  }else if(args[argon].compare("-M")==0 &&(argon+1<args.size())) {
    modeltype=getModelType(args[argon+1]);
    args.erase(args.begin()+argon,args.begin()+argon+2);
  }else if(args[argon].compare("-a")==0 &&(argon+1<args.size())) {
    double alpha=atof(args[argon+1].c_str());
    prior->setAlpha(alpha);
    args.erase(args.begin()+argon,args.begin()+argon+2);
  }else if(args[argon].compare("-k")==0 &&(argon+1<args.size())) {
    opt.treescale=atoi(args[argon+1].c_str());
    args.erase(args.begin()+argon,args.begin()+argon+2);
  }else if(args[argon].compare("-C")==0 &&(argon+1<args.size())) {
    maxconcordance=atoi(args[argon+1].c_str());
    args.erase(args.begin()+argon,args.begin()+argon+2);
  }else if(args[argon].compare("-B")==0 &&(argon+1<args.size())) {
    tmp=args[argon+1];
    string bms;
    int bm=getBetaModel(tmp,bms);
    if(opt.userefonly || dref!=NULL){
      prior->setPriorRef(bm);
      // Don't allow complicated priors for reference beta model
    }else{
      prior->setPriorCounts(bm);
      betapriorstring=bms;
    }
    // 
    args.erase(args.begin()+argon,args.begin()+argon+2);
  }else if(args[argon].compare("-b")==0 &&(argon+1<args.size())) {
    //
    priorstring=args[argon+1];
    args.erase(args.begin()+argon,args.begin()+argon+2);
  }else if(args[argon].compare("-T")==0 &&(argon+1<args.size())) {
    treetype=atoi(args[argon+1].c_str());
    if(treetype>3 && treetype<7) {treemodification=TREEMOD_NOFLATTEN;treetype-=3;}
    if(treetype>=7) {treemodification=TREEMOD_FLATTENHILLCLIMB;treetype-=6;}
    args.erase(args.begin()+argon,args.begin()+argon+2);
  }else if(args[argon].compare("-c")==0 &&(argon+1<args.size())) {
    corfactor=atof(args[argon+1].c_str());
    args.erase(args.begin()+argon,args.begin()+argon+2);
  }else if(args[argon].compare("-F")==0 &&(argon+1<args.size())) {
    fixfile=args[argon+1];
    args.erase(args.begin()+argon,args.begin()+argon+2);
  }else if(args[argon].compare("-t")==0 &&(argon+1<args.size())) {
	opt.test_max=atoi(args[argon+1].c_str());
    args.erase(args.begin()+argon,args.begin()+argon+2);
  }else if(args[argon].compare("-X")==0) {
    xhead=true;
    args.erase(args.begin()+argon);
  }else if(args[argon].compare("-Y")==0) {
    yhead=true;
    args.erase(args.begin()+argon);
  }else if(args[argon].compare("-K")==0) {
    opt.fixK=true;
    args.erase(args.begin()+argon);
  }else if(args[argon].compare("-v")==0) {
    opt.verbose=1;
    args.erase(args.begin()+argon);
  }else if(args[argon].compare("-S")==0) {
    opt.verbose=0;opt.silent=1;
    args.erase(args.begin()+argon);
  }else if(args[argon].compare("-h")==0) {
    cout<<help<<endl;return 0;
  }else if(args[argon].compare("-H")==0) {
    cout<<help<<endl<<fullhelp<<endl;return 0;
  }else if(args[argon].compare("-V")==0) {
    printVersion();  return 0;
  }else if(args[argon][0]=='-') {
    cout<<"Wrong arguments: did not recognise "<<args[argon]<<" (does i require an argument?)"<<endl<<help<<endl;return 1;
  }else{
    argon++;
  }
}

    if(priorstring.compare("")!=0){// We have a prior string to apply to the active dataset
      if(dref!=NULL){
	prior->setPriorRef(priorstring);
      }else{
	prior->setPriorCounts(priorstring);
      }
      priorstring="";
    }


  if (args.size()<1) {cout<<help<<endl;return 0;}
  my_rng * rng = seedrng(seed,opt.verbose);// 0 means use /dev/random or clock.
  comment.append("\nSeed: ");
  ss<<seed;
  comment.append(ss.str());

// process the main arguments
// read data
  if (args.size()>=2) {
    if(opt.verbose) {
      cout<<"Opening data file: "<<args[0];
      if(opt.userefonly) cout<<" in reference mode."<<endl;
      else cout<<" in standard all-vs-all mode."<<endl;
    }
    try{
      dsafe=new Data(args[0],ignorelines,xhead,yhead,1-opt.userefonly);
      if(opt.userefonly==0) {dcount=dsafe;
      }else dref=dsafe;
    }catch(std::string x){cout<<x<<endl;exit(0);}
  }else {cout<<"Need data and an output file!"<<endl<<help<<endl; return 0;}

  if(fixfile.size()>0) {
    if(opt.verbose) cout<<"Using fixed population file "<<fixfile<<endl;
    try{
      if(dcount!=NULL) dcount->makeSuperFromFile(fixfile);
      if(dlength!=NULL) dlength->makeSuperFromFile(fixfile);
      if(dref!=NULL) dref->makeSuperFromFile(fixfile);
    }catch(std::string x){cout<<x<<endl;exit(0);}
  }
    
// assign the corfactor
    if(corfactor>0 && dsafe->getCfactor()<0) {cerr<<"WARNING: You have specified 'C' and provided a datafile that does not contain it."<<endl<<"If you know what you are doing you can ignore this warning."<<endl<<"Otherwise you are advised to use the 'chromocombine' tool to estimate it."<<endl<<"See www.paintmychromosomes.com under 'ChromoCombine' for details."<<endl;
	}else if(dsafe->getCfactor()<0) {cerr<<"WARNING: 'c' NOT READ FROM DATA INPUT FILE."<<endl<<"This means that fineSTRUCTURE may expect to see the wrong variance"<<endl<<"in the data and you will probably experience poor clustering."<<endl<<"You are advised to use the 'chromocombine' tool to estimate 'C'."<<endl<<"See www.paintmychromosomes.com under 'ChromoCombine' for details."<<endl;
	}
    if(corfactor>0) {
      if(dcount==NULL) setcorfactorref=true;
      else setcorfactor=true;
      dsafe->setCfactor(corfactor);
    }
    if(dcount!=NULL) dcount->ensureCfactor();
    if(dlength!=NULL) dlength->ensureCfactor();
    if(dref!=NULL) dref->ensureCfactor();

    
// decide which data to use
/*    if(dlength==NULL) datainference=INFDATA_COUNTS;// no data for the lengths
    else {*/
	string sdat=opt.usedata;
	std::transform(sdat.begin(), sdat.end(),sdat.begin(), ::toupper);
	if(sdat.substr(0,1).compare("L")==0) datainference=INFDATA_LENGTHS;
	else if(sdat.substr(0,1).compare("C")==0) datainference=INFDATA_COUNTS;
	else if(sdat.substr(0,1).compare("T")==0) datainference=INFDATA_TOTALLENGTHS;
	else if(sdat.substr(0,1).compare("A")==0) datainference=INFDATA_ALL;
	else if(sdat.substr(0,1).compare("D")==0) datainference=INFDATA_ALLNOTLENGTHS;
	else datainference=INFDATA_ALLNOTLENGTHS;
	// otherwise we use all INFDATA_ALLNOTLENGTHS, the default
//    }
// assign the bvec hyperprior for betamodel (lengths are added later)
//	bvec=getBvec(betamodel,betamodelref,datainference,bvec,betapriorstring,ignorelines,xhead,yhead);
    
    State * state;
// read state if appropriate, or create a new one
  if (args.size()>2) {
    fs=string(args[1]);
    if(opt.verbose) cout<<"Reading state file: "<<fs.c_str()<<endl;
    FsXml infile(fs);
    try{
      if(infile.gotoLineContaining("<Iteration>",true)<0) {// is a population file
	state=new State(rng,dcount,fs,prior,true,dlength,dref,datainference,modeltype);
      }else{// is an xml output file
	state=new State(rng,dcount,&infile,prior,false,dlength,dref,datainference,modeltype);
	havefullxmlinput=true;
      }
    }catch(std::string x){cout<<"ERROR: "<<x<<endl;exit(0);}

  }else {
	if(opt.verbose) cout<<"Creating state from data."<<endl;
	try{state=new State(rng,dcount,initpop,prior,dlength,dref,datainference,modeltype);}catch(std::string x){cout<<"ERROR: "<<x<<endl;exit(0);}
  }
  cout<<"Processing."<<endl; 
//
    if(extract && (ext.substr(0,2).compare("TR")==0 )){ // we must have an extra argument, the third line being the tree we READ IN, the 4th being the output file
      if (args.size()>2)  {
	fstree=args[1];
	cout<<"Reading tree file "<<fstree<<endl;
      }else{
	cerr<<"For extract -e TREE, need 4 files in the order data mcmc tree output"<<endl;
	exit(0);
      }
    }
// create output file
    filebuf fb;
    string oname=args[args.size()-1];
    try{
      fb.open (oname.c_str(),ios::out);
    }catch(std::string x){
	cerr<<"Error opening file!"<<endl<<x<<endl; return 0;}
    ostream os (&fb);
// Extract only
    if(extract){
	std::transform(ext.begin(), ext.end(),ext.begin(), ::toupper);
	if(ext.substr(0,2).compare("X2")==0) {// normalised copying matrix
	  if(args.size()<2) {
	    cout<<"Error: -e X2 needs <chunkcounts.out> <optional mcmc.xml> <output>"<<endl; return 0;
	  }
	  if(fs2.size()>0) setState(rng,fs2,dsafe,prior,dlength,dref,datainference,modeltype).printX(&os,true);
		else state->printX(&os,true);
	}else if(ext.substr(0,1).compare("X")==0) {// unnormalised copying matrix
	  if(args.size()<3) {
	    cout<<"Error: -e X needs <chunkcounts.out> <mcmc.xml> <output>"<<endl; return 0;
	  }
	  if(fs2.size()>0) setState(rng,fs2,dsafe,prior,dlength,dref,datainference,modeltype).printX(&os,false);
		else state->printX(&os,false);
	}else if(ext.substr(0,1).compare("B")==0) {// beta matrix
	  if(args.size()<2) {
	    cout<<"Error: -e Beta needs <chunkcounts file> <optional mcmc xml file> <outputfile>"<<endl; return 0;
	  }
	  if(fs2.size()>0) setState(rng,fs2,dsafe,prior,dlength,dref,datainference,modeltype).printBeta(&os);
		else state->printBeta(&os);
	}else if(ext.substr(0,3).compare("MEA")==0) {// meancoincidence (pairwise) matrix
	  if(args.size()<3) {
	    cout<<"Error: -e meancoincidence needs <chunkcounts.out> <mcmc.xml> <output>"<<endl; return 0;
	  }
	FsXml *infile=new FsXml(fs);
	InfExtract iext(rng,dsafe,infile,opt.verbose);
	delete(infile);
	infile=new FsXml(fs);
	InfExtract2 iext2(rng,dsafe,infile,prior,false,dlength,dref,datainference,modeltype);
	if(fs2.size()>0){
	  State ords=setState(rng,fs2,dsafe,prior,dlength,dref,datainference,modeltype);
		try{iext.reorder(ords.allIndInOrder());
		}catch(string x){cout<<"ERROR:reorder: "<<x<<endl;exit(0);}
	}
//	iext.getMeanX();
	iext.printMeanX(&os);
	}else if(ext.substr(0,2).compare("MI")==0){// minimum distance state
	  if(args.size()<3) {
	    cout<<"Error: -e minimumdistance needs <chunkcounts.out> <mcmc.xml> <output>"<<endl; return 0;
	  }
		FsXml infile(fs);
		size_t firstcolon=ext.find_first_of(":");
		double pen=1.0;
		if(firstcolon!=string::npos) { pen=atof(ext.substr(firstcolon+1).c_str());}
		InfExtract iext(rng,dsafe,&infile,opt.verbose);
		try{state=new State(rng,dsafe,iext.getMeanX(),true, 0.0,prior,dlength,dref,datainference,modeltype);}catch(std::string x){cout<<"ERROR: "<<x<<endl;exit(0);}
		iext.makeMinSquaresState(pen,state);//iext.getState()
		iext.getState()->setprint(&os);	
	}else if(ext.substr(0,1).compare("P")==0) {// donor file
	  if(args.size()<3) {
	    cout<<"Error: -e popidfile needs <chunkcounts.out> <mcmc.xml or tree.xml> <output>"<<endl; return 0;
	  }
	  FsXml infile(fs);
	  InfExtract2 iext(rng,dsafe,&infile,prior,opt.verbose,dlength,dref,datainference,modeltype);
	  state=iext.getState();
	  InfExtractDonor iextd(state,opt.verbose);
	  iextd.printDonorFile(&os);
	}else if(ext.substr(0,3).compare("MER")==0) {
	  if(args.size()<3) {
	    cout<<"Error: -e merge needs <chunkcounts.out> <mcmc.xml> <output>"<<endl; return 0;
	  }
		double mergeval=0.95;
		bool mergerule=true;
		size_t firstcolon=ext.find_first_of(":");
		string opts;
		if(firstcolon!=string::npos) {
			opts=ext.substr(firstcolon+1);
			size_t secondcolon=opts.find_first_of(":");
			mergeval=atof(opts.substr(0,secondcolon).c_str());
			if(secondcolon!=string::npos){
			opts=opts.substr(secondcolon+1);
			if(opts.find_first_of("sS0t")!=string::npos) mergerule=false;
			}
		}
		FsXml infile(fs);
		InfExtract iext(rng,dsafe,&infile,opt.verbose);
		try{state=new State(rng,dsafe,iext.getMeanX(),mergerule, mergeval,prior,dlength,dref,datainference,modeltype);}catch(std::string x){cout<<"ERROR: "<<x<<endl;exit(0);}
		state->setprint(&os);
	}else if(ext.substr(0,3).compare("MAX")==0) {//maximum posterior state
	  if(args.size()<3) {
	    cout<<"Error: -e maxstate needs <chunkcounts.out> <mcmc.xml> <output>"<<endl; return 0;
	  }
	  FsXml infile(fs);
	  InfExtract2 iext(rng,dsafe,&infile,prior,opt.verbose,dlength,dref,datainference,modeltype);
	  state=iext.getState();
	  state->setprint(&os);
	}else if(ext.substr(0,1).compare("A")==0 ){// extract admixture matrix
		State * adstate;
		State tmpstate=setState(rng,fs2,dsafe,prior,dlength,dref,datainference,modeltype);
		if(fs2.size()>0) adstate= new State(&tmpstate);
		else adstate=new State(state);
		InfAdmixture infad(rng,dsafe,adstate,1.0,opt.verbose);
		infad.printQs(&os);
	}else if(ext.substr(0,2).compare("PM")==0 ){// extract P matrix
	  if(args.size()<3) {
	    cout<<"Error: -e pmatrix needs <chunkcounts.out> <mcmc.xml> <output>"<<endl; return 0;
	  }
		State * adstate;
		State tmpstate=setState(rng,fs2,dsafe,prior,dlength,dref,datainference,modeltype);
		if(fs2.size()>0) adstate= new State(&tmpstate);
		else adstate=new State(state);
		InfAdmixture infad(rng,dsafe,adstate,1.0,opt.verbose);
		infad.printPs(&os);
	}else if(ext.substr(0,2).compare("TH")==0 ){//thin
	  if(args.size()<3) {
	    cout<<"Error: -e thin:<step>: needs <chunkcounts.out> <mcmc.xml> <output.xml>"<<endl; return 0;
	  }
		size_t firstcolon=ext.find_first_of(":");
		string opts;
		long thin=1;
		if(firstcolon!=string::npos) {
			thin=atoi(ext.substr(firstcolon+1).c_str());
		}else{
	    cout<<"Error: -e thin:<step>: needs <chunkcounts.out> <mcmc.xml> <output.xml>"<<endl; return 0;
                }

		if(opt.verbose)cout<<"Thinning of "<<fs<<" by "<<thin<<endl;
		FsXml infile(fs);
		infile.thin(thin,&os);
	}else if(ext.substr(0,1).compare("R")==0 ){//extract range
	  if(args.size()<3) {
	    cout<<"Error: -e range:<from>:to needs <chunkcounts.out> <mcmc.xml> <output.xml>"<<endl; return 0;
	  }
		size_t firstcolon=ext.find_first_of(":");
		string opts;
		long from=0,to=0;
		if(firstcolon!=string::npos) {
			opts=ext.substr(firstcolon+1);
			size_t secondcolon=opts.find_first_of(":");
			from=atof(opts.substr(0,secondcolon).c_str());
			if(secondcolon!=string::npos){
			to=atoi(opts.substr(secondcolon+1,string::npos).c_str());
			}else{
	    cout<<"Error: -e range:<from>:to needs <chunkcounts.out> <mcmc.xml> <output.xml>"<<endl; return 0;
                }
		}else{
	    cout<<"Error: -e range:<from>:to needs <chunkcounts.out> <mcmc.xml> <output.xml>"<<endl; return 0;
                }
		if(opt.verbose)cout<<"Range extraction of "<<fs<<" from "<<from<<" to "<<to<<endl;
		FsXml infile(fs);
		infile.rangeextract(from,to,&os);
	}else if(ext.substr(0,1).compare("L")==0 ) {// Likelihood sample
	  if(args.size()<3) {
	    cout<<"Error: -e likelihood needs <chunkcounts.out> <mcmc.xml> <output>"<<endl; return 0;
	  }
		FsXml infile(fs);
		try{InfExtract5 iext5(rng,dsafe,&infile,prior,1,false);
		vector<double> liks=iext5.getLikelihoods();
		for(unsigned long c1=0;c1<liks.size()-1;c1++) os<<setprecision(12)<<liks[c1]<<", ";
 		if(liks.size()>0)os<<setprecision(12)<<liks[liks.size()-1]<<endl;	
		}catch(string x){cerr<<"Error in extraction of likelihoods:"<<x<<endl;exit(0);}
	}else if(ext.substr(0,2).compare("PR")==0 ){// Posterior Probability
	  if(args.size()<3) {
	    cout<<"Error: -e probability needs <chunkcounts.out> <mcmc.xml> <output>"<<endl; return 0;
	  }
		FsXml infile(fs);
		try{InfExtract4 iext4(rng,dsafe,&infile,prior,false);
		vector<double> posteriors=iext4.getPosteriors();
		for(unsigned int c1=0;c1<posteriors.size()-1;c1++) os<<setprecision(12)<<posteriors[c1]<<", ";
 		if(posteriors.size()>0)os<<setprecision(12)<<posteriors[posteriors.size()-1]<<endl;	
		}catch(string x){cerr<<"Error in extraction of probabilities:"<<x<<endl;exit(0);}
	}else if(ext.substr(0,2).compare("ST")==0 ){// Stats
		FsXml infile(fs);
		try{InfExtract4 iext4(rng,dsafe,&infile,prior,false);
		vector<double> posteriors=iext4.getPosteriors();
		for(unsigned int c1=0;c1<posteriors.size()-1;c1++) os<<setprecision(12)<<posteriors[c1]<<", ";
 		if(posteriors.size()>0)os<<setprecision(12)<<posteriors[posteriors.size()-1]<<endl;	
		}catch(string x){cerr<<"Error in extraction of probabilities:"<<x<<endl;exit(0);}
	}else if(ext.substr(0,2).compare("TR")==0 ){// tree in newick format
	  if(args.size()<3) {
	    cout<<"Error: -e tree needs <chunkcounts.out> <mcmc.xml> <tree.xml> <output>"<<endl; return 0;
	  }
	  Inf1 *tmptree = GlobalReadTree(rng,dsafe,fstree,prior,dlength,dref,
						  datainference,modeltype,opt.verbose+opt.silent*2);
	      tmptree->printTree(&os,false);
	}else {cerr<<"Error: invalid extraction."<<endl;}
	// end of extract functions
/// MCMC MODEL
    }else if(opt.method.compare(0,1,"M")==0|| opt.method.compare(0,1,"O")==0) {
      if(opt.verbose) cout<<"Starting MCMC calculation."<<endl;
      InfMCMC * infMCMC=GlobalRunMCMC(rng,state,&os,opt.burnin,opt.additional,opt.thinin,comment,datainference,opt.pcaprob,opt.fixK,opt.verbose+opt.silent*2);
/*
	try{InfMCMC infMCMC(d,state,dlength,datainference,opt.verbose+opt.silent*2);
	if(opt.verbose) cout<<"BURN IN PHASE"<<endl;
	infMCMC.metropolis(0,opt.burnin);
	if(opt.verbose) cout<<"MCMC PHASE"<<endl;
	infMCMC.resetCounters();
	infMCMC.metropolis(opt.burnin,opt.additional,opt.thinin,&os);
*/
	try{
	State * state2=new State(infMCMC->getState());
	if(opt.method.compare(0,1,"M")==0 ) {
	  Inf1 inf2(rng,state2,datainference,opt.verbose+opt.silent*2,opt.test_max,opt.treescale);
	  if(opt.verbose) cout<<"TREE CREATION PHASE"<<endl;
	  try{inf2.mergeHillClimb(&os,false,treemodification);}catch(std::string x){cout<<x<<endl;exit(0);}
	}
	infMCMC->exportXmlTail(&os);
	delete(infMCMC);
	delete(state2);
	}catch(std::string x){cout<<"Tree creation error: "<<x<<endl;}
/// (SEMI) DETERMINISTIC MODELS
    }else if (opt.method.compare(0,1,"S")==0) {// split tree
      cerr<<"ERROR: This method has been temporarily disabled. If you think that it would be useful for your research, contact dan.lawson@bristol.ac.uk detailing your needs and we will reimplement it properly."<<endl;
      exit(0);
      // //      Inf1 inf1(rng,initpop,datainference,modeltype,opt.verbose+opt.silent*2,opt.test_max);
      // 	inf1.exportXmlHead(&os,fs,string("SplitTree"),opt.burnin);
      // 	inf1.exportXmlComment(&os,comment);
	
      // 	if(initpop>0) {
      // 	  if(opt.verbose) cout<<"SPLIT PHASE"<<endl;
      // 	  try{inf1.splitHillClimb(true);}catch(std::string x){cout<<x<<endl;exit(0);}
      // 	}

      // 	State * state2=new State(inf1.getState());

      	// Inf1 inf2(rng,state2,datainference,opt.verbose+opt.silent*2,opt.test_max);
      	// if(opt.verbose) cout<<"MERGE PHASE"<<endl;
      	// try{inf2.mergeHillClimb(&os,false,treemodification);}catch(std::string x){cout<<x<<endl;exit(0);}
      	// inf1.exportXmlTail(&os);
/// MAIN TREE METHOD
    }else if (opt.method.compare(0,1,"T")==0) {// merge tree
      /*
	if(opt.verbose) cout<<"MERGE PHASE"<<endl;
	FsXml *infile=new FsXml(fs);
	InfExtract iext(d,infile,opt.verbose);
	Inf1 * inf1_i=NULL;
	delete(infile);
	infile=new FsXml(fs);
	InfExtract2 * iext2=NULL;
	InfMCMC * infHillClimb=NULL;
	State * state2;
	if(havefullxmlinput){
		iext2=new InfExtract2(d,infile,prior,opt.verbose,dlength,datainference);
		if(treetype==1){
			state2=new State(iext2->getState());
			inf1_i=new Inf1(d,iext2->getState(),dlength,datainference,opt.verbose,opt.test_max);
		}else if(treetype==2) {
			inf1_i=new Inf1(d,iext2->getState(),dlength,datainference,opt.verbose,opt.test_max);
			try{inf1_i->mergeHillClimb(NULL,true,false);}catch(std::string x){cout<<x<<endl;exit(0);}
			state2=new State(inf1_i->getState());

		}else if(treetype==3) {
			try{infHillClimb=new InfMCMC(d,iext2->getState(),dlength,datainference,opt.verbose);
			//infHillClimb->hillClimb(0,opt.burnin);
			infHillClimb->metropolis(0,opt.burnin);
			state2=new State(infHillClimb->getState());
			}catch(std::string x){cout<<x<<endl;exit(0);}
		}else {cerr<<"Invalid tree type (-T option)."<<endl;exit(0);}
		delete(infile);
	}else {
		if(treetype==1){
			inf1_i=new Inf1(d,state,dlength,datainference,opt.verbose,opt.test_max);
			state2=new State(state);
		}else if(treetype==2) {
			inf1_i=new Inf1(d,state,dlength,datainference,opt.verbose,opt.test_max);
			try{inf1_i->mergeHillClimb(NULL,true,false);}catch(std::string x){cout<<x<<endl;exit(0);}
			state2=new State(inf1_i->getState());

		}else if(treetype==3){
			try{infHillClimb=new InfMCMC(d,state,dlength,datainference,opt.verbose);
			infHillClimb->hillClimb(0,opt.burnin);
			state2=new State(infHillClimb->getState());
			}catch(std::string x){cout<<x<<endl;exit(0);}
		}else {cerr<<"Invalid tree type (-T option)."<<endl;exit(0);}
	}
	*/
	long newx=-1,newy=-1,newz=-1;
	string olddatafile;
	try{
	  double newc=-1.0;
	  double newcref=-1.0; //
	  getXmlHeader(fs,newc,newcref,newx,newy,newz,olddatafile);
	  if(!setcorfactor && newc>=0 && dcount!=NULL) dcount->setCfactor(newc);
	  if(!setcorfactorref && newcref>=0 &&dref!=NULL) dref->setCfactor(newcref); 
	}catch(string x){
	  cerr<<x<<endl;
	}
	compareDataFiles(olddatafile, dsafe->getFileName()); // produces warnings if there are different file names
//cout<<"treetype="<<treetype<<" (data) (fs) treetestmax="<<opt.test_max<<" hcs="<<opt.burnin<<";
//cout<<" di="<<datainference<<" mt="<<modeltype<<" (state) (dlength) hfxml="<<havefullxmlinput<<" fixK="<<opt.fixK<<" ts="<<opt.treescale<<endl;
	Inf1* inf1=mergeTree(rng,treetype,dcount, fs,opt.test_max,opt.burnin,prior,
			    datainference,modeltype, state, dlength,dref, havefullxmlinput,opt.fixK,opt.treescale,oname,maxconcordance,opt.verbose+opt.silent*2);
//	Inf1 inf1(state2,dlength,datainference,opt.verbose,opt.test_max);
	inf1->exportXmlHead(&os,fs,string("MergeTree"),opt.burnin);
	inf1->exportXmlComment(&os,comment);
	InfMCMC* tmcmc=new InfMCMC(rng,inf1->getState(),INFDATA_COUNTS,0,false);//Used to print the state
	tmcmc->exportXmlIter(&os,0); // export the iteration
//	inf1->getState()->iterPrint(&os);
	try{inf1->mergeHillClimb(NULL,false,treemodification);}catch(std::string x){cout<<x<<endl;exit(0);}
	if(opt.verbose) cout<<"Assigning certainty"<<endl;
	FsXml *infile=new FsXml(fs);
	InfExtract3 iext3(rng,dsafe,infile,inf1->getNodes(),opt.verbose);
	delete(infile);
	/*
	// NOTE: Diagonalise disabled due to errors with force files.
	if(opt.verbose) cout<<"Diagonalise tree"<<endl;
	inf1.diagonaliseOnM(d->getMatrix(),false);
*/
	if(opt.verbose) cout<<"Finish up"<<endl;

/*	if(inf1_i!=NULL) {
		inf1.reorderState(inf1_i->getState());
		inf1_i->getState()->iterPrint(&os);
	}else if(infHillClimb!=NULL){
		inf1.reorderState(infHillClimb->getState());
		infHillClimb->getState()->iterPrint(&os);
	}
*/
	inf1->printTree(&os);
	inf1->exportXmlTail(&os);
	delete(inf1);
    }else if(opt.method.compare(0,1,"A")==0) {// admixture model
	bool atest=false;
	if(opt.method.compare(0,10,"ADMIXTURET")==0) atest=true;
	if(opt.verbose &&!atest) cout<<"Admixture model..."<<endl;
	else if(opt.verbose &&atest) cout<<"Admixture model test only..."<<endl;
	//state->setprint(&cout);
	InfAdmixture infad(rng,dsafe,state,2,opt.verbose,atest);
//	infad.printPs(&cout);
	try{
	if(opt.verbose) cout<<"BURN IN PHASE"<<endl<<flush;
	infad.metropolis(0,opt.burnin);
	infad.exportXmlHead(&os);
	infad.exportXmlComment(&os,comment);
	if(opt.verbose) cout<<"MCMC PHASE"<<endl<<flush;
	infad.resetCounters();
	infad.metropolis(opt.burnin,opt.additional,opt.thinin,&os);
	infad.exportXmlTail(&os);
	}catch(std::string x){cerr<<"Error in admixture:"<<endl<<x<<endl;}
    }else {
	cerr<<"Invalid method."<<endl<<help<<endl;
    }
    freerng(rng);
    if(dsafe!=NULL) delete(dsafe);
    delete(prior);
    fb.close();
    if(!opt.silent)cout<<fssuccesstext<<endl;
    return 0;
}


} // end namespace fines
