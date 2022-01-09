
#include <vector>
#include <string>
#include <cstring>
#include <iostream>

#include "ChromoCombineFuns.h"
#include "ChromoCombine.h"

using namespace std;

static const char * help=
    "\
    Usage: chromocombine [OPTIONS] -o <outputfileroot> <listofinputfiles>\n\
	Remember to set <outputfileroot>.\n\
    	Additionally, you must supply at least one inputfile.\n\
	Inputfiles can be specified in the following forms:\n\
	  a) ChromoPainter root file name\n\
	  b) ChromoPainter *full* file names, with any combination of endings.\n\
	     Repeated roots will be kept only once.\n\
	Note that you must not have changed the file name endings.\n\
    Options:\n\
    -o <outputfileroot>	Set the output file root name (default: output)\n\
    -d 			Specify that the input is a directory, and that *all*\n\
			valid file roots ending in <ending> are of interest.\n\
    -l			Specify that length matrices should be IGNORED.\n\
    -m			Specify that mutation matrices should be IGNORED.\n\
    -c			Specify that counts should be IGNORED.\n\
    -f <forcefile>	Process a continental forcing file.  \"c\" will be \n\
			calculated for this file separately, and stored in it.\n\
			The other output files will not be changed by this.\n\
			Therefore, fineSTRUCTURE will use the correct \"c\" when\n\
			run on the dataset without a forcefile, as well as when\n\
			one is provided. You must however use the force file\n\
			on the dataset for which \"c\" was calculated.\n\
    -F <forceoutput>	By default, the force file is *updated* with the correct\n\
			value of c.  This instead writes a *new* force file that\n\
			contains the important sections (in case you want to use the.\n\
			same force file on many datasets).\n\
    -e <ending>		Set the inputfile ending (default: \".out\")\n\
			Remember the dot!\n\
    -E <ending>		Set the outputfile ending (default: \".out\")\n\
    -C			Instead of smallish regions, use complete genomic segments for\n\
			calculation of c.  (Currently only valid if you provide a \n\
			number of input files and that each has approximately the same recombination distance)\n\
    -i <file>		The id file from ChromoPainter to create consistent columns. This is needed to provide the correct column headers if donorfiles have been used to exclude some individuals, but the individuals were themselves the populations.\n\
    -t			Test: print the files that will be processed.\n\
    			and additionally test that they exist.\n\
    -u			Unsafe mode; when summing individuals from different files, \n\
			the default is to stop if some are seen more often than others.\n\
			Use this option to override.\n\
    -v			Verbose mode.\n\
    -q                  Quiet mode. Report only terminal problems.\n\
    -h			This help message. \n\
    ";




int chromocombine(int argc, char *argv[])
{
  //  cout<<"CHROMOCOMBINE: Arguments:"<<endl;
  //  for(int c1=0;c1<argc;++c1){
  //    cout<<argv[c1]<<endl;
  //  }
  ChromoCombineWorker ccw;
  ccw.outputfileroot=string("output");
  ccw.outputfileending=string(".out");
  ccw.outputfileoutending=string(".out");
  ccw.testrun=0;
  ccw.numeffinds=0;
  ccw.wantmuts=true;
  ccw.wantcounts=true;
  ccw.wantlengths=true;
  ccw.useforce=false;
  ccw.verbose=0;
  ccw.completegenomes=false;
  bool usedir=false;
  bool saferun=true;
  char c;
  optind=1;
  while ((c = getopt (argc, argv, "o:e:E:f:F:i:vCudclmthq")) != -1) switch (c){
    case('o'):ccw.outputfileroot=string(optarg);break;
    case('e'):ccw.outputfileending=string(optarg);break;
    case('E'):ccw.outputfileoutending=string(optarg);break;
    case('f'):ccw.forcefile=string(optarg);ccw.useforce=true;
      try{ccw.makeSuperFromFile(ccw.forcefile);}
      catch(string x){cout<<x<<endl<<help<<endl;return(0);};
      break;
    case('F'):ccw.forcefileoutput=string(optarg);break;
    case('i'):if(ccw.readIds(string(optarg))==0) {
	printf("ReadIds Error: no ids found!\n");exit(1);
      }
      break;
    case('d'):usedir=true;break;
    case('c'):ccw.wantcounts=false;break;
    case('l'):ccw.wantlengths=false;break;
    case('m'):ccw.wantmuts=false;break;
    case('t'):ccw.testrun=1;break;
    case('u'):saferun=false;break;
    case('C'):ccw.completegenomes=true;break;
    case('v'):ccw.verbose=1;break;
    case('q'):ccw.verbose=-1;break;
    case('h'):cout<<help<<endl;return 0;
    default:
      cout<<"Failure to read argument: "<<c<<endl;
	    cout<<help<<endl;
            return(-1);
  }
  if (argc-optind<1) {
    cout<<"Failure to read arguments: Found "<<argc<<" formal arguments and "<<argc-optind<<" files."<<endl;
   cout<<help<<endl;return 0;
  }
  ccw.makeEndings();

  while(optind<argc){
    ccw.inputfilerootraw.push_back(string(argv[optind++]));
  }
  if(usedir) {
    try{
      ccw.getFilesFromDirectories();
    }catch(string x){
      cout<<x<<endl;
      cout<<help<<endl;
      return(0);
    }
  }
  ccw.testFilesExist();
  
  if(ccw.verbose>0){
    cout<<"Input files "<<endl;
    for(unsigned long c1=0;c1<ccw.inputfilerootraw.size();c1++) cout<<"... "<<ccw.inputfilerootraw[c1]<<endl;
  }
  
  ccw.rationaliseInputFiles();
  ccw.testFilesExist();
  
  if(ccw.testrun>0){
    if(ccw.verbose>0) cout<<"Files to be processed:"<<endl;
    for(unsigned long c1=0;c1<ccw.inputfileroot.size();c1++) {
      cout<<ccw.inputfileroot[c1];
      if(ccw.validfile[c1]==0) cout<<" ... OK!"<<endl;
      else cout<<" ... found "<<ccw.requiredendings.size()-ccw.validfile[c1]<<" of "<<ccw.requiredendings.size()<<" files!"<<endl;
    }
    return(0);
  }
  if(ccw.verbose>0)cout<<"Reading files"<<endl;
  long readfilesval=ccw.readFiles();
  if(readfilesval<=0){
    cerr<<"Error reading files!"<<endl;
    return(0);
  }
  if(ccw.verbose>0)cout<<"Read files!"<<endl;
  try{
    ccw.getForceIds();
  }catch(string x){
    cout<<x<<endl<<help<<endl;
    return(0);
  }
  long validation=ccw.validateMatrices();
  if(validation<0){
    cerr<<"Cannot continue with invalid matrix structure!"<<endl;
    return(0);
  }else if(validation>0 && saferun){
    cerr<<"Not continuing with possibly unbalanced input files.  To override, run with \"-u\" flag."<<endl;
    return(0);    
  }
  if(ccw.verbose>0)cout<<"Validated files!"<<endl;
  ccw.finaliseData();
  if(ccw.verbose>0)cout<<"Finalised data!"<<endl;
  if(ccw.wantcounts) {
    ccw.cval=ccw.calcC();
    if(ccw.verbose>0)cout<<"Calculated C!"<<endl;
  }
  ccw.writeFiles();
  cout<<"Successfully summed "<<ccw.inputfileroot.size()<<" file root(s) containing "<<ccw.indexcounts.size()<<" individuals, to new file root "<<ccw.outputfileroot<<endl;
  cout<<"Successfully written chunkcounts file "<<ccw.outputfileroot<<ccw.possiblemiddles[0]<<ccw.outputfileoutending;
  if(ccw.wantcounts)cout <<" with c value "<<ccw.cval;
  cout<<endl;
  if(ccw.useforce){
    if(ccw.verbose>0)cout<<"applying force..."<<endl;
    try{ccw.applyForceFile();
    }catch(std::string x){cout<<"Error applying force file"<<endl;exit(0);}
    double forceC=ccw.calcC(true);
    if(ccw.writeForceFile(forceC)){
      string outfile=ccw.forcefileoutput;
      if(outfile.size()==0)outfile=ccw.forcefile;
      cout<<"Successfully written forcefile "<<outfile<<" with c value "<<forceC<<endl;
    }else{
      cerr<<"Warning: Failed to write continental force file with updated c!"<<endl;
    }
  }
  if(ccw.cval>0 && ccw.wantcounts) {
    cout<<"When using these files with fineSTRUCTURE there is now no need to specify \"c\"."<<endl;
  }else if(ccw.wantcounts){
    cout<<"WARNING: \"c\" has not been correctly estimated. You can either specify it manually (not advised), or follow the advise above and in the faq at www.paintmychromosomes.com to solve the problem with its estimation."<<endl;
  }
  cout<<ccsuccesstext<<endl;
  return(0);
}
