#include <string>
#include <cstring>
#include <fstream>
#include <sstream>

#include "cp/ChromoPainterMutEM.h"
#include "finestructure/fines.h"
#include "chromocombine/ChromoCombine.h"
#include "fsproject.h"
#include "fsconstants.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


using namespace std;

namespace fines
{

/////////////////////////////////////
string getFsVersion(){
	string ret;
	ret.append("mixPainter v1.0.0");
	ret.append("\n");
	#ifdef PACKAGE_STRING // Omitted due to version matching issues
	ret.append("Report all bugs to ");
	ret.append(PACKAGE_BUGREPORT);
	ret.append("\n");
        #endif
	ret.append("mixPainter");
	ret.append(" build date "); 
	ret.append(__DATE__);
	ret.append(" at ");
	ret.append(__TIME__);
	return(ret);
}

  void badpainterhelp(){
    cout<<"mixPainter help"<<endl;
    cout<<"Usage: "<<endl<<" mixPainter <directory name> -p <phase files> -i <id file> [-r <recomb files> -H -f <fraction> -c <max cores> -v -h]"<<endl;
    cout<<"<directory name>: The location of the directory where mixPainter can store intermediate files."<<endl;
    cout<<"<phase files>: A list of CHROMOPAINTER phase format files. See README for details."<<endl;
    cout<<"<id file>: A file listing two columns, the \"ind_name population_name\" of each individual. Also accepts standard fs format id file, which includes a third usage column containing 0/1, 1 for individuals that are included in this analysis."<<endl;
    cout<<"OPTIONAL ARGUMENTS:"<<endl;
    cout<<"-r <recombination files>: A list of CHROMOPAINTER recomb files. See README for details. Can be generated with makeuniformrecfile.pl or convertrecfile.pl. If this is not specified, mixPainter will use unlinked painting, which is not recommended."<<endl;
    cout<<"-H: Specify haploid mode. Default: Diploid mode (meaning 1 individual has two rows in the phase file)."<<endl;
    cout<<"-f <fraction>: Use only <fraction> of the individuals to estimate the CHROMOPAINTER parmameters with. 0.1 is a good choice and is the default. Set to 1 for best estimation."<<endl;
    cout<<"-c: Maximum number of processing cores to use. Default: autodetect, and use all available."<<endl;
    cout<<"-h: This help."<<endl<<endl;
    cout<<"Description:"<<endl<<"mixPainter runs CHROMOPAINTER in a mode suitable for input into badMIXTURE.  It can be parallelized across multiple cores in a single computer but not across multiple machines (i.e. no HPC mode.)"<<endl;
    return;
  }

  void freeargvout(char *argvout[]){
    // free up memory
    for(int i=0;i<1000;i++)
      {
	free(argvout[i]);
      }
    free(argvout);
  }
  
} // end namespace fines


using namespace fines;

int main(int argc, char *argv[])
{
    string comment="Command line: ";
    for(int c1=0;c1< argc;c1++) {comment.append(argv[c1]);comment.append(" ");}
    comment.append("\nVersion: ");
    comment.append(getFsVersion());
    string mydirectory="";
    string mycpfile="";

    if (argc<1) {
      badpainterhelp();
      return(1);
    }

    // allocate memory for the commands we'll pass to fs
    std::vector<std::string> args(argv, argv+argc) ;
    args.erase(args.begin());

    if(args[0].substr(0,2).compare("-V")==0){ // version
      cout<<getFsVersion()<<endl;
      return(0);
    }
    if(args[0].substr(0,2).compare("-h")==0){ // help (on purpose)
      badpainterhelp();
      return(0);
    }
    
    if (argc<2) {
      badpainterhelp();
      return(1);
    }
  

    if(args[0].substr(0,1).compare("-")==0){ // help!
      badpainterhelp();
      return(1);
    }else{
      mydirectory=args[0];
      mycpfile=args[0];
      mycpfile.append(".cp");
      args.erase(args.begin());
    }

    int argcout=0;
    char **argvout;
    argvout=(char**)malloc(sizeof(char*)*1000);
    for(int i=0;i<1000;i++)
      {
	argvout[i] = (char*)malloc(1000*sizeof(char));
	strcpy( argvout[i], "" );
      }
    strcpy(argvout[argcout++],"fs");
    strcpy(argvout[argcout++],mycpfile.c_str());
    strcpy(argvout[argcout++],"-N");
    
    // things that we need to keep track of
    bool verbose=false;
    bool s6indfrac=false;
    bool idfile=false;
    bool phasefiles=false;
    bool recombfiles=false;
    // Process the remaining parameters
    unsigned int argon=0;
    while(args.size()>0) {
      argon=0;
      //      cout<<"Assessing argument "<<args[argon]<<" at location "<<argon<<endl;
      if(args[argon].compare("-h")==0 || args[0].compare("-help")==0){ // help!
	badpainterhelp();
	freeargvout(argvout);
	return 0;
	//      }else if(args[argon].compare("-v")==0) {	// verbose
	//	verbose=true;
	//	cout<<"Verbose mode"<<endl;
	//	strcpy(argvout[argcout++],"-v");
	//	argon++;
      }else if(args[argon].compare("-H")==0){
	cout<<"Enabling Haploid mode"<<endl;
	strcpy(argvout[argcout++],"-ploidy");
	strcpy(argvout[argcout++],"1");
	argon++;
      }else if(args[argon].compare("-c")==0){
	if(argon>=args.size()-1) {
	  cerr<<"Error: Invalid options to -c!"<<endl;
	  freeargvout(argvout);
	  return(1);
	}
	cout<<"Setting CPU cores to "<<args[argon+1].c_str()<<endl;
	strcpy(argvout[argcout++],"-numthreads");
	strcpy(argvout[argcout++],args[argon+1].c_str());
	argon+=2;
      }else if(args[argon].compare("-f")==0){
	s6indfrac=true;
	if(argon>=args.size()-1) {
	  cerr<<"Error: Invalid options to -f!"<<endl;
	  freeargvout(argvout);
	  return(1);
	}
	cout<<"Setting parameter estimation fraction to "<<args[argon+1].c_str()<<endl;
	strcpy(argvout[argcout++],"-s6indfrac");
	strcpy(argvout[argcout++],args[argon+1].c_str());
	argon+=2;
      }else if(args[argon].compare("-i")==0){
	idfile=true;
	if(argon>=args.size()-1) {
	  cerr<<"Error: Invalid options to -i!"<<endl;
	  freeargvout(argvout);
	  return(1);
	}
	cout<<"Setting id file to "<<args[argon+1].c_str()<<endl;
	strcpy(argvout[argcout++],"-popidfile");
	strcpy(argvout[argcout++],args[argon+1].c_str());
	argon+=2;
      }else if(args[argon].compare("-p")==0){
	phasefiles=true;
	if(argon>=args.size()-1) {
	  cerr<<"Error: Invalid options to -p!"<<endl;
	  freeargvout(argvout);
	  return(1);
	}
	cout<<"Assigning phase files:"<<endl;
	strcpy(argvout[argcout++],"-phasefiles");
	while(argon+1<args.size() && args[argon+1].substr(0,1).compare("-")!=0) {
	  cout<<".. Using phase file "<<args[argon+1].c_str()<<endl;
	  strcpy(argvout[argcout++],args[++argon].c_str());
	}
	argon++;
      }else if(args[argon].compare("-r")==0){
	recombfiles=true;
	if(argon>=args.size()-1) {
	  cerr<<"Error: Invalid options to -r!"<<endl;
	  freeargvout(argvout);
	  return(1);
	}
	cout<<"Assigning recomb files:"<<endl;
	strcpy(argvout[argcout++],"-recombfiles");
	while(argon+1<args.size() && args[argon+1].substr(0,1).compare("-")!=0) {
	  cout<<".. Using recomb file "<<args[argon+1].c_str()<<endl;
	  strcpy(argvout[argcout++],args[++argon].c_str());
	}
	argon++;
      }else{
	cout<<"Error: Unknown argument: "<<args[argon]<<endl;
	  freeargvout(argvout);
	  return(1);
      }
      //      cout<<"Deleting "<<argon<<" elements"<<endl;
      args.erase(args.begin(), args.begin() + argon);
      //      cout<<"Remaining size "<<args.size()<<endl;
    }

    if(!idfile){
      cerr<<"ERROR: -i <idfile> not specified!"<<endl;
      freeargvout(argvout);
      return(1);
    }

    if(!phasefiles){
      cerr<<"ERROR: -p <phasefile> not specified!"<<endl;
      freeargvout(argvout);
      return(1);
    }
    
    if(!s6indfrac && recombfiles){
      strcpy(argvout[argcout++],"-s6indfrac");
      strcpy(argvout[argcout++],"0.1");      
    }
    strcpy(argvout[argcout++],"-s67args:-J");

    strcpy(argvout[argcout++],"-s7chunksperregion");
    strcpy(argvout[argcout++],"1");

    strcpy(argvout[argcout++],"-indsperproc");
    strcpy(argvout[argcout++],"1");
    
    strcpy(argvout[argcout++],"-allowdep");
    strcpy(argvout[argcout++],"0");
    strcpy(argvout[argcout++],"-countdatapop");
    strcpy(argvout[argcout++],"-combines5");
    if(recombfiles){
      cout<<"Detected running in LINKED mode. This is the recommended usage."<<endl;
      strcpy(argvout[argcout++],"-makes6");
      strcpy(argvout[argcout++],"-dos6");
      strcpy(argvout[argcout++],"-combines6");
    }else{
      cout<<"Detected running in UNLINKED mode. This is not recommended where dense data are available."<<endl;
    }
    cout<<""<<endl;
    
    strcpy(argvout[argcout++],"-makes7");
    strcpy(argvout[argcout++],"-dos7");
    strcpy(argvout[argcout++],"-combines7");
    
    
    // testing
    
    // for(int i=0;i<argcout;++i){
    //   cout<<"Arg "<<i<<" : "<<argvout[i]<<endl; 
    // }
    
    int ret=fsproject(0,argcout,argvout);

    if(ret==0){
      cout<<"Successfully run CHROMOPAINTER to create a dataset suitable for mixture analysis."<<endl;
      cout<<"You can now delete the folder "<<mydirectory<<", or keep it for the records of which individuals were used in which chromopainter run."<<endl;
      cout<<"You should put the following file into badMIXTURE:"<<endl;
      if(recombfiles){
	cout<<mydirectory<<"_pop_linked.chunkcounts.out"<<endl;
      }else{
	cout<<mydirectory<<"_pop_unlinked.chunkcounts.out"<<endl;
      }
    }else{
      cout<<"Something went wrong!"<<endl;
    }
    return(ret);
}
