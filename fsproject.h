#ifndef FSPROJECT_H
#define FSPROJECT_H
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <fstream>
#include <string.h>

#include "fscmds.h"
#include "fssettings.h"
#include "fsparam.h"
//#include "finestructure/fsxml.h"

using namespace std;
using namespace fines;

/**
    @brief Project mode
*/
class FsProject
{
 public:
  FsProject(std::string f,std::string d,bool verbose);///< Create from a file
  ~FsProject();///< Destructor

  void setupSections();///< Setup the file structure for the project (comments and freetext comments)
  void addHistory(std::vector<std::string> args);///< Adds the command with all arguments to the history file
  void defineParameters();///< Defines the parameters we can set
  void defineCommands();///< Defines the commands we can run
  void readFromFile();///< Read the project from the file
  bool applyVal(FsSettingsValue val);///< Read a specific value from file
  bool applyVal(std::vector<std::string> args);///<

  void writeToFile();////< Write the project to the file
  string writeToFile(string tpar);///< Return a string containing a parameter
  void deleteStage(string ending);///< Delete the stageX file (can be a number, 2a or 7a)
  void safeRemoveCommandFile(int forstage);///< Backup and remove the command file for a given stage
  void backupFileAndRemove(string filename,string backupdir="cpbackup");///< Backup a file to the cpbackup directory, then delete it
  void backupFile(string filename,string backupdir="cpbackup");///< Backup a file to the cpbackup directory
  void safeCreateFile();///< Create the project file, moving files of the same name out of the way
  void setFileName(std::string f);///< Set the filename of the project file
  void setDirectoryName(std::string d);///< Set the filename of the project file
  void setVerbose(bool verbose);///< Set the verbosity level
  void setFsmode(int fsmode);///< Whether we are in classic fs mode
  std::string getFileName();///< obtain the filename being used
  std::string getDirectoryName();///< obtain the directory being used
  std::string getExec();///< Obtain the executable being used
  
  int defaultChunksperregion();///< Default number of chunks per region
  void copyfile(std::string from, std::string to); ///< Copy a file
  void resetToStage(int newstage);///< Reset settings to a stage by removing some settings from future stages
  void createDuplicated(int newstage, string newname);///< Duplicate the appropriate parts of the file structures to enable everything that can be reused, to be reused.
  int applyFiles(string readfile,string actionifdirexists=string("merge"));///<create directories and read files as appropriate to "readfile" (either "new","read", or "detect"). If it exists we either "merge", "delete" or "stop"
  bool haveOutput(int stage,bool forcombined=false);///< check that we have the output we expect from the a stage
  string whichLinkagemode();///< Which linkage mode should we use
  void countSnps();///< Count and simply check the SNPs and recombination files
  void countData();///< Count and simply check the datasets we've been given
  void countDataPop();///< Count and simply check the population datasets we've been given (call countData if not already done)
  void makePopIds();///< make the population id files for admixture mode
    
  bool canDo(string cmd,bool vverbose=false);///< check dependencies for a command
  int parIndex(std::string partest);///< get the index of a parameter (-1 for missing)
  int cmdIndex(std::string cmdtest);///< get the index of a command (-1 for missing)
  std::string cmdInfo(std::string cmd, bool statetype=true);///< return the information/help for a command
  std::vector<std::string> getDependencies(std::vector<std::string> args);///< Gets the previous command that needs to be run
  bool checkArgs(std::vector<std::string> args,int minargs,int maxargs=-1);///< Checks that the command described by args has between minargs and maxargs arguments. maxargs=-1 means maxargs=minargs, -2 means unbounded
  int checkStage(std::string args,std::string val);/// return the stage for a *parameter*. Return 0 on a problem and 1 on OK
  void checkStage(std::string args,int maxstage);///< Check that args is allowed to run at the current stage
  void checkStage(std::vector<std::string> args,int maxstage);///< Check that args[0] is allowed to run at the current stage
  void docmd(std::string args);///< do a command line command
  void docmd(std::vector<std::string> args);///< do a command line command
  std::vector<char *> converttoargv(std::string s); ///< Convert a string into a argv
  void freeargv(std::vector<char *> ccmd); ///< Free the argv array
  
  void createIdFile(string idfile);///< Creates an ID file from the first stage file
  int numRecipients(int forstage=-1);///< Get the number of recipient individuals that are being painted for the requested stage; or by default from the id file (for standard fs) or from the donor file (otherwise)
  
  bool finishedStage(int stage);/// check if we have finished a stage

  void ensureCpEmRoot(int forstage); ///< Makes sure we have a stage 1/6 root file name
  void ensureCpRunRoot(int forstage); ///< Makes sure we have a stage 2 root file name
  void ensureCpRunARoot(int forstage); ///< Makes sure we have a stage 2a root file name
  void ensureCommandfile(int forstage); ///< Makes sure we have a stage 1 root file name

  std::string makeoutfileroot(string root,int fileon,int indstart,int indend,int forstage);///< Contruct the output file from the run details

  void matchCpInputs();///< Check that we have valid input files for cp 
  int getNhapsFromFile(int which);///< Extract the number of haps from a stage file
  int getNsnpsFromFile(unsigned int fileno);///< Extract the number of SNPs from a file
  int getNindsFromFile(bool keepall=true,string idf="");///< count the number of recipients to process (in file idf, default is the idfile)
  int getUniqueNindsFromFile(bool keepall=true,string idf="");///< count the number of unique recipients to process (in file idf, default is the idfile)
  void makeStageCpEM(int forstage);///< Make the commands to be run in stage1
  void makeStageCpRun(int forstage);///< Make the commands to be run in stage2

  void doCpStage(int stage);///< Run the commands to be run in stage1/2
  void doStagePopExtract();///< Run Stage 5, extracting populations from the FS tree
  

  bool validateStage(int forstage,bool vverbose);///< Validate a stage
  bool validateStage(int forstage,string ending,int expectedlines,string lastline,double minfrac,bool vverbose);///< Validates the output of a stage, by checking files with a given ending for an expected number of lines or a known lastline string.
  vector<int> validOutput(int forstage,string ending,int expectedlines,string lastline,bool vverbose); ///< Returns a vector of missing (-1) incomplete (0) or successful (1) for each command in the s<x>outputvec
  bool validateStageCpEm(int forstage,bool vverbose=false); ///< validate an EM stage
  bool validateStageCpRun(int forstage,bool vverbose=false); ///< Validate a painting stage

  void removeCompletedCommands(int forstage);///< Remove all commands that are completed for a given stage
  void combineCpEm(int fromstage);///< Combine the outputs of stage1 to estimate parameters
  void combineStageCpRun(int fromstage);///< Do chromocombine (all chromosomes plus each one individually)
  string combineStageCpRunSingle(int chr,string sroot, string sxcombineargs, vector<string>srootvec,int forstage);///< Do chromocombine (specified chromosome set only). Return the logfile.
  void combineStage3();///< Do MCMC validation
  void combineStage4();///< Do tree validation
  void writeHpcStage3(string cmdfile);///< Do the mcmc in an iterative manner
  
  void addEMline(string line, vector<double> *Nevec_ptr, vector<double> *muvec_ptr); ///< add the EM line to Ne and mu
  void addEMobs(string filename, vector<double> *Nevec_ptr, vector<double> *muvec_ptr);///< Extract all the final EM lines from Ne and mu
  void writeStringVectorToFile(std::vector<std::string> sv,std::string fn,std::vector<std::string> logfiles, bool addfs=true);///< write a string vector to file (optionally adding "fs " to the start)
  void writeStringToFile(std::string s,std::string fn);///< write a string to file

  ////////////////////////////////////
  void do_omp_set_num_threads(int numthreads);///< Set the number of threads, if omp is available.
  int get_omp_get_thread_num(); ///< Get the number of this thread, if running parallel  (0 else)
  int get_omp_get_num_threads(); ///< Get the number of threads available, if running parallel (0 else)

  //////////////////////////////////// 
  void continuefsmcmc();///< rerun the mcmc, continuing where we left off
  void ignoreGRfsmcmc();///< restore a previous MCMC and set the GR threshold to ignore any problems with it
  void ensurefsroot();///< Make sure wr have a finestructure root
  void ensurefsmcmc();///< Make sure we have a finestructure file name
  void ensurefstree();///< Make sure we have a finestructure tree file name
  void makefsmcmc();///< Make the fs mcmc command lines
  void makefstree();///< Make the fs tree command lines
  void dofsmcmc();///< Do the fs mcmc
  void dofstree();///< Do the fs mcmc

  ////////////////////////////////////
  void switchStdout(const char *newStream);///< switch stdout to newstream
  void revertStdout();///< Revert Stdout to console
  int getHpc();///< Whether or not we are in hpc mode
  int getIndsPerProc(int forstage=-1);///< Figure out how many inds to use per process
  string getCommandFile(int forstage=-1);///< Command file name for a specified (-1: or current) stage
  int getCommandFileCount(int forstage=-1);///< Command file count for a specified (-1: or current) stage
  int getCommandFileLength(int forstage=-1);///< Command file content (number of lines) for a specified (-1: or current) stage
  void optionsHelp();///< Default help message about automatic mode
  void outputHelp();///< Help message about the created files
  void stagesHelp();///< Help message about the various computational stages
  void getHelp(std::vector<std::string> args);///<Complex help on particular topics
  void makeExample();///< Create walkthrough

  bool readmcmctraces(int filenum);///< Read the MCMC traces for a particular run number, stored into a temporary (unsaved) object
  bool writemcmctraces(std::string filename);///< Write the MCMC traces to file
  bool mcmcConvergence(); ///< Obtain MCMC convergence diagnostics from the mcmc output, returning 1 if it is OK to proceed, 0 if we are concerned
  double mcmcGRstatistic(std::vector<std::vector<double> > *data); // Compute the Gelman Rubin potential scale reduction factor statistic for a set of runs, for a parameter
  int recommendN();///< Recommend the number of processes per node
  int recommendM();///< Recommend the number of commands per node

  // GT stuff
  void makeGt(); ///< Makes everything that is needed for a GLOBETROTTER analysis
  void ensureGtRoot();///< Ensure that we have a well-defined GT root name, for output file construction
  void ensureDonors();///< Make sure that we have donor information
  void zcatSamples(vector<string> srcfiles, string destfile); ///< Zcat all the samples into a single file, as needed by GT
  void ensureCombinedSamples();///< Combine the samples into a single file for each chromosome, as required by GT
  void cleanupSamples();///< Remove the per-individual samples files, for saving space
  void makeDataListFiles();///< Makes the list of samples and recombination files that are required by GT
  void writeGtParams(string recip);///< Write a GT parameter file
  
 protected:
  ////////////////////////////////////////////
  // Metadata 
  int nstages;///< The number of stages we store
  bool restorestage;///<Whether we are restoring the stage
  bool allowdep;///< Whether we allow dependencies to be autoresolved
  bool verbose; ///< Whether we are in verbose mode
  int fsmode; ///< Whether we are in classic fs mode and should do normal reporting (1=true if yes, 0=false if no)
  std::string filename;
  std::string dirname;
  std::string fileroot;
  std::vector<std::string> sectionnames; ///< Names of sections
  std::vector<std::string> sectioncomments; ///< Names of sections
  std::vector<FsPar> pars; ///< Names of parameters
  std::vector<int> parsize;///< Maximum index of parameters in each section
  std::vector<FsCmd> cmds; ///< Names of commands

  std::vector<std::string> freetextcomments; ///< Comments between tags
  std::string historytext; ///< history

  ////////////////////////////////////////////
  // data in computer processed form
  // processing files (not stored in the settings file)
  std::string s1commandfile;///< Where we store stage 1 commands
  std::string s2commandfile;///< Where we store stage 2 commands
  std::string s3commandfile;///< Where we store stage 3 commands
  std::string s4commandfile;///< Where we store stage 4 commands
  std::string s6commandfile;///< Where we store stage 6 commands
  std::string s7commandfile;///< Where we store stage 7 (panel painting) commands
  std::string s9commandfile;///< Where we store stage 9 (gt) commands

  std::vector<std::string> s1commands; ///< commands we construct for stage 1
  std::vector<std::string> s2commands; ///< commands we construct for stage 2
  std::vector<std::string> s3commands; ///< commands we construct for stage 3
  std::vector<std::string> s4commands; ///< commands we construct for stage 4
  std::vector<std::string> s6commands; ///< commands we construct for stage 6
  std::vector<std::string> s7commands; ///< commands we construct for stage 7
  std::vector<std::string> s9commands; ///< commands we construct for stage 9

  std::vector<std::string> s1logfiles; ///< the logfile locations to be written to file
  std::vector<std::string> s2logfiles; ///< the logfile locations to be written to file
  std::vector<std::string> s3logfiles; ///< the logfile locations to be written to file
  std::vector<std::string> s4logfiles; ///< the logfile locations to be written to file
  std::vector<std::string> s6logfiles; ///< the logfile locations to be written to file
  std::vector<std::string> s7logfiles; ///< the logfile locations to be written to file
  std::vector<std::string> s9logfiles; ///< the logfile locations to be written to file

  // for storing and restoring stdout
  int stdout_fd; 
  fpos_t stdout_pos;


  ////////////////////////////////////////////
  // universal properties
  int stage;///< the stage that we are currently processing.
  string fsfile; ///< file root used for storing directories and the fs details
  int hpc; ///< HPC mode
  int ploidy; ///< haploid or diploid mode
  int numthreads; ///< Number of threads
  string paintercmd;///< Command to use for chromopainter. The default, internal choice is "cp" which is special (called by fs cp and has --emfilesonly). Changing it uses an external command, presumably ChromoPainterv2, and forces hpc mode

  int indsperproc; ///< how we split up the individuals; default is 1, meaning everyone is processed in a different command  
  string linkagemode; ///< Whether we use linked or unlinked mode (or default: autodetect)
  bool outputlogfiles;///< Whether the redirection to log files are included in the command lists
  std::string exec; ///< the full path of this program, including the right ending (for specifying versions)
  std::vector<int> validatedoutput; ///< Whether we have validated the output of each stage
  std::vector<int> missingfiles; ///< How many files that we have counted for each stage that are not present or complete

  /////////////////////////

  
  // STAGE0 OUTPUT / STAGE 1 PROCESSING AND INPUT
  
  ////////////////////////////////////////////
  // Stage12 UNIVERSAL QUANTITIES
  std::string s12inputtype; ///< type of input files
  std::vector<std::string> phasefiles; ///< input phase files
  std::vector<std::string> recombfiles; ///< input recombination files
  std::string idfile; ///< input donor file
  std::string s12donorfile; ///< input donor file
  std::string s12args; ///< input arguments
  int ninds;///< Number of individuals for which we have data
  int nindsUsed;///< Number of individuals to be used from the data
  int nrecipsDonor;///< number of recipient individuals in the donor file
  int nsnps; ///< Number of SNPs found in total
  vector<int> nsnpsvec;///< Number of SNPs in each file

  ////////////////////////////////////////////
  // STAGE1 OUTPUT / STAGE 1a INPUT
  std::string s1args; ///< stage 1 arguments
  int s1emits;///< number of em iterations 
  std::string s1outputroot;///< Base output root
  std::vector<std::string> s1outputrootvec; ///< output files
  int s1minsnps;///< Minimum number of loci extracted for an EM block
  double s1snpfrac;///< `Fraction' of SNPs to process
  double s1indfrac;///< `Fraction' of INDS to process
  // STAGE1a OUTPUT
  
  ////////////////////////////////////////////
  // STAGE2 OUTPUT / STAGE 2a INPUT
  double Neinf; // Ne as we infer it from stage1
  double muinf; // mu as we infer it from stage1
  int s2chunksperregion;///< Number of chunks to be used to define a region
  int s2samples;///< Number of samples of the painting to obtain per recipient haplotype (default to zero, i.e. don't obtain these)

  std::string s2args; ///< input arguments
  std::string s2outputroot;///< Base output root
  std::vector<std::string> s2outputrootvec; ///< output files
  std::string s2combineargs;///< stage 2combine arguments

  ////////////////////////////////////////////
  // STAGE 3 AND 4 PARAMS 

  // STAGE2a OUTPUT
  double cval; ///< The inferred value of "c"
  std::string cproot; ///<combined output root
  std::string cpchunkcounts; ///<chromopainter final chunkcount file

  // STAGE3-4 generic properties
  string s34args;
  std::string fsroot; ///<finestructure output root

  // STAGE 3 PARAMS and OUTPUT (finestructure)
  static const int numitersDefault=100000;
  int s3iters; ///< User interface: how many iterations to do. By default, half are assigned to burnin, half to mcmc
  int s3iterssample; // -x negative for autocalculate from s3iters
  int s3itersburnin; // -y same as above
  int numskip; ///< -z same as above
  int maxretained; ///< for calculating -z
  int nummcmcruns; ///< number of mcmc runs performed

  std::string fsmcmcoutput; ///<fs name
  std::vector<std::string> fsmcmcoutputvec; ///<fs mcmc files
  std::vector<double> mcmcGR;///< The Gelman Rubin statistics computed from mcmc files
  double threshGR; ///<Threshold for rejecting convergence

  // STAGE 4 PARAMS and OUTPUT (finestructure tree)

  string s4args;
  int s4iters; // -x negative for autocalculate from s3iters
  std::string fstreeoutput; ///<fs name
  std::vector<std::string> fstreeoutputvec; ///<fs mcmc files

  // STAGE 5 OUTPUT / STAGE 6 INPUT
  std::string popidfile; ///< input reference id file, containing populations (which may be different from the idfile, e.g. constructed from the clustering)
  std::string donoridfile;///< input reference donor file (optional, can be constructed as all vs all)
  std::vector<std::string> popcmdargs; ///< the arguments that we will add to the chromopainter run
  //  std::vector<std::string>  idfiles;///< list of all the id files that we will use
  std::vector<std::string>  refidfiles; ///< list of all the reference id files that we have created for external use
  std::string refdonorfile; ///< the donor file that we have created for external use
  
  // STAGE6 OUTPUT / STAGE 6a INPUT
  int s6emits;///< number of em iterations 
  std::string s6args; ///< input arguments
  std::string s67args; ///< input arguments
  std::string s6outputroot;///< Base output root
  std::vector<std::string> s6outputrootvec; ///< output files
  int s6minsnps;///< Minimum number of loci extracted for an EM block
  double s6snpfrac;///< `Fraction' of SNPs to process
  double s6indfrac;///< `Fraction' of INDS to process

  // STAGE 6 output
  double popNeinf; // Ne as we infer it from stage6
  double popmuinf; // mu as we infer it from stage6

  // STAGE7 OUTPUT / STAGE 7a INPUT
  std::string s7combineargs;///< stage 7combine arguments
  int s7chunksperregion;///< Number of chunks to be used to define a region
  int s7samples;///< Number of samples of the painting to obtain per recipient haplotype (default to zero, i.e. don't obtain these)
  std::string s7args; ///< input arguments
  std::string s7outputroot;///< Base output root
  std::vector<std::string> s7outputrootvec; ///< output files

  // STAGE7a OUTPUT
  double popcval; ///< The inferred value of "c"
  std::string popcproot; ///<combined output root
  std::string popcpgenomelen; ///<chromopainter final genome length file


  // STAGE 9 = GT stuff 
  std::string gtroot;///< The root name given to all GT files
  std::string gtsampleslist;///< The list of all samples files
  std::string gtrecomblist;///< The list of all recombination files

  ////< Parameters passed to GT in the par file
  std::string gtbinwidth;
  std::string gtnummixingiterations;
  std::string gtbootstrapnum;
  std::string gtpropscutoff;
  
  std::vector<std::string> gtdonors; ///< The donor populations read from the donor file
  std::vector<std::string> gtrecips; ///< The recips populations read from the donor file
  std::vector<std::string> gtsamplesfiles;///< The samples files for each chromosome
  std::vector<std::string> gtparamfiles; ///< The parameter files that GT will call //********************** ADD TO EXPORT *********//
  
  // Temporary objects for calculations
  std::vector<std::vector<double> > mcmc_posterior;
  std::vector<std::vector<double> > mcmc_k;
  std::vector<std::vector<double> > mcmc_beta;
  std::vector<std::vector<double> > mcmc_delta;
  std::vector<std::vector<double> > mcmc_f;

  // Keeping track of old MCMC runs
  std::vector<std::string> old_fsmcmcoutputvec; ///<fs mcmc files
  int fscounter;
  int fsmaxattempts;
};

int fsproject(int fsmode, int argc, char *argv[]); ///< fsproject

#endif
