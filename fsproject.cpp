#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdexcept>
#include <exception>
#include <iostream>
#include <limits>
#include <string>
#include <algorithm>

/////////////////////////////
// Not cross platform?
#ifdef _OPENMP
#include <omp.h>
#endif
// Also need to think about mkdir from unistd.h

#include <sys/types.h>
#include <sys/stat.h>
/////////////////////////////

#include "cp/ChromoPainterMutEM.h"
#include "chromocombine/ChromoCombine.h"
#include "finestructure/fines.h"
#include "finestructure/fsxml.h"
#include "fssettings.h"
#include "fsutils.h"
#include "fsproject.h"
#include "fsconstants.h"
#include "fsdonor.h"
#include "alphanum.hpp"

using namespace std;


//////////////////////////////////////
FsProject::FsProject(std::string f,std::string d,bool verbose) 
{
  exec="fs";
  numthreads=0;
  nstages=9;
  allowdep=true;
  indsperproc=0;
  stage=0;
  outputlogfiles=true;

  validatedoutput=vector<int>(nstages,0);
  missingfiles=vector<int>(nstages,-1); // missing files from each stage; default: all
  
  s1commandfile=string("");
  s2commandfile=string("");
  s3commandfile=string("");
  s4commandfile=string("");
  s6commandfile=string("");
  s7commandfile=string("");
  s9commandfile=string("");
  
  // data things
  nsnps=0;
  ninds=0;
  nindsUsed=0;
  nrecipsDonor=0;  
  // chromopainter things
  ploidy=2;
  hpc=0;
  linkagemode=string("unlinked");
  s12args=string("");
  s1args=string("-in -iM --emfilesonly");
  s1minsnps=10000;// default: 10K
  s1snpfrac=0.1;// fraction of genome we will use for EM
  s1indfrac=1.0;// fraction of inds we will use for EM
  s1emits=10;// number of EM iterations
  s2args=string("");
  s2combineargs=string("");
  s34args=string("-X -Y");
  Neinf=-1; // default: we haven't got estimates
  muinf=-1; // default: we haven't got estimates
  s2chunksperregion=-1; // default: can be wrong for unlinked data!
  s2samples=0; // default: we don't want samples
  cval=-1;// inferred "c" (we don't use this except for information, we set from the files)
  paintercmd="cp";// chromopainter command. cp is the internal one, everything else must use hpc mode
  
  // finestructure things
  s3iters=numitersDefault;
  maxretained=500;
  s3iterssample=-1;
  s3itersburnin=-1;
  numskip=-1;
  nummcmcruns=2;

  fscounter=0;
  fsmaxattempts=5;
  threshGR=1.3;
  // finestructure tree things
  s4iters=numitersDefault;
  s4args=string("");

  // population - based chromopainter things
  popNeinf=-1; // default: we haven't got estimates
  popmuinf=-1; // default: we haven't got estimates
  // 
  s6minsnps=s1minsnps;// minimum number of SNPs to use
  s6snpfrac=s1snpfrac;// fraction of inds to use
  s6indfrac=s1indfrac;// fraction of inds to use
  s6emits=s1emits;// number of EM iterations
  s6args=string("-in -iM --emfilesonly");
  s67args=string("");
  s7combineargs=string("");
  s7chunksperregion=-1; // default: can be wrong for unlinked data!
  s7samples=0; // default: we don't want samples
  s7args=string("");
  //
  popcval=-1;// inferred "c" (we don't use this except for information, we set from the files)

  gtbootstrapnum="20";
  gtpropscutoff="0.001";
  gtbinwidth="0.1";
  gtnummixingiterations="5";
  
  this->verbose=verbose;
  setFileName(f);
  setDirectoryName(d);
  setupSections();

  defineParameters();
  defineCommands();
}


FsProject::~FsProject() 
{
}
///////////////////////////////

void FsProject::defineParameters()
{
  ////// Universal properties
  sectionnames.clear();
  sectioncomments.clear();
  pars.clear();

  sectionnames.push_back("0");
  sectioncomments.push_back("Data preparation and conversion. Not yet implemented");
  pars.push_back(FsPar("validatedoutput",-1,"Derived. Whether we have validated output from each stage of the analysis (0-8)",1));
  pars.push_back(FsPar("missingfiles",-1,"Derived. Number of output files that we think are incorrect from each stage of the analysis (0-8)",1));
  parsize.push_back(pars.size());

  sectionnames.push_back("1");
  sectioncomments.push_back("Universal Stage1-4 properties of the analysis");
  pars.push_back(FsPar("exec",-1,"Finestructure command line. Set this to be able to use a specific version of this software. (default: fs)"));
  pars.push_back(FsPar("hpc",-1,"THIS IS IMPORTANT FOR BIG DATASETS! Set hpc mode. 0: Commands are run 'inline' (see 'numthreads' to control how many CPU's to use). 1: Stop computation for an external batch process, creating a file containing commands to generate the results of each stage. 2: Run commands inline, but create the commands for reference. (default: 0.)"));
  pars.push_back(FsPar("numthreads",-1,"Maximum parallel threads in 'hpc=0' mode. Default: 0, meaning all available CPUs."));
  pars.push_back(FsPar("paintercmd",-1,"Command to use for chromopainter. The default, internal choice is \"cp\" which is special (called by fs cp and has --emfilesonly). Changing it uses an external command, presumably ChromoPainterv2, and forces hpc mode."));
  pars.push_back(FsPar("ploidy",-1,"Haplotypes per individual. =1 if haploid, 2 if diploid. (default: 2)"));
  pars.push_back(FsPar("linkagemode",-1,"unlinked/linked. Whether we use the linked model. default: unlinked / linked if recombination files provided."));
  pars.push_back(FsPar("indsperproc",-1,"Desired number of individuals per process (default: 0, meaning autocalculate: use 1 In HPC mode, ceiling(N/numthreads) otherwise. Try to choose it such that you get a sensible number of commands compared to the number of cores you have available."));
  pars.push_back(FsPar("outputlogfiles",-1,"1=Commands are written to file with redirection to log files. 0: no redirection. (default:1)"));
  pars.push_back(FsPar("allowdep",-1,"Whether dependency resolution is allowed. 0=no, 1=yes. Main use is for pipelining. (default:1)."));
  parsize.push_back(pars.size());

  /// Chromopainter Stage1-2 generic properties
  sectionnames.push_back("1");
  sectioncomments.push_back("ChromoPainter Stage1-2 generic properties");
  pars.push_back(FsPar("s12inputtype",1,"What type of data input (currently only \"phase\" supported)"));
  pars.push_back(FsPar("idfile",1,"IDfile location, containing the labels of each individual. REQUIRED, no default (unless -createids or -popidfile is used)."));
  pars.push_back(FsPar("s12args",1,"arguments to be passed to Chromopainter (default: empty)"));
  parsize.push_back(pars.size());

  //Quantities observed from data
  sectionnames.push_back("1");
  sectioncomments.push_back("Quantities observed from data");
  pars.push_back(FsPar("ninds",-1,"Derived. number of individuals observed in the idfile",1));
  pars.push_back(FsPar("nindsUsed",-1,"Derived. number of individuals retained for processing from the idfile",1));
  pars.push_back(FsPar("nrecipsDonor",-1,"Derived. number of individuals retained as recipients from the donor file",1));
  pars.push_back(FsPar("nsnps",-1,"Derived. number of SNPs in total, over all files",1));
  parsize.push_back(pars.size());

  //chromopainter stage1 (EM)
  sectionnames.push_back("1");
  sectioncomments.push_back("ChromoPainter Stage1 (EM) properties");
  pars.push_back(FsPar("s1args",1,"Arguments passed to stage1 (default:-in -iM --emfilesonly)"));
  pars.push_back(FsPar("s1emits",1,"Number of EM iterations (chromopainter -i <n>, default: 10)"));
  pars.push_back(FsPar("s1minsnps",1,"Minimum number of SNPs for EM estimation (for chromopainter -e, default: 10000)"));
  pars.push_back(FsPar("s1snpfrac",1,"fraction of genome to use for EM estimation. (default: 0.1)"));
  pars.push_back(FsPar("s1indfrac",1,"fraction of individuals to use for EM estimation. (default: 1.0)"));
  pars.push_back(FsPar("s1outputroot",1,"output file for stage 1 (default is autoconstructed from filename)"));
  parsize.push_back(pars.size());

  //chromopainter inferred properties for stage2 from stage1
  sectionnames.push_back("1");
  sectioncomments.push_back("ChromoPainter Stage2 properties inferred from Stage1");
  pars.push_back(FsPar("Neinf",-1,"Derived. Inferred `Effective population size Ne' (chromopainter -n).",1));
  pars.push_back(FsPar("muinf",-1,"Derived. Inferred Mutation rate mu (chromopainter -M)",1));
  parsize.push_back(pars.size());
  
  //chromopainter stage2 (painting)
  sectionnames.push_back("2");
  sectioncomments.push_back("ChromoPainter Stage2 (main run) properties");
  pars.push_back(FsPar("s2chunksperregion",2,"number of chunks in a \"region\" (-ve: use default of 100 for linked, nsnps/100 for unlinked)"));
  pars.push_back(FsPar("s2samples",2,"number of samples of the painting to obtain per recipient haplotype, for examining the details of the painting. (Populates <root>.samples.out; default 0. Warning: these file can get large)"));
  pars.push_back(FsPar("s2args",2,"Additional arguments for stage 2 (default: none, \"\")"));
  pars.push_back(FsPar("s2outputroot",2,"Output file name for stage 2 (default: autoconstructed)."));
  pars.push_back(FsPar("s2combineargs",2,"Additional arguments for stage 2 combine (fs combine; default: none, \"\")"));
  parsize.push_back(pars.size());

  // chromocombine inferred properties for stage3 from stage2
  sectionnames.push_back("2");
  sectioncomments.push_back("FineSTRUCTURE Stage3 properties inferred from Stage2");
  pars.push_back(FsPar("cval",-1,"Derived. 'c' as inferred using chromopainter. This is only used for sanity checking. See s34 args for setting it manually.",1));

  pars.push_back(FsPar("cproot",3,"The name of the final chromopainter output. (Default: <filename>, the project file name)"));
  pars.push_back(FsPar("cpchunkcounts",3,"the finestructure input file, derived name of the chunkcounts file from cproot.",1));
  parsize.push_back(pars.size());

  // finestructure stage3-4 generic properties
  sectionnames.push_back("3");
  sectioncomments.push_back("FineSTRUCTURE Stage3-4 generic properties");
  pars.push_back(FsPar("fsroot",3,"The name of the finestructure output (Default: <filename>, the project file name)."));
  pars.push_back(FsPar("s34args",3,"Additional arguments to both finestructure mcmc and tree steps. Add \"-c <val>\" to manually override 'c'."));
  parsize.push_back(pars.size());

  // finestructure stage3 MCMC inference
  sectionnames.push_back("3");
  sectioncomments.push_back("FineSTRUCTURE Stage3 MCMC inference");
  pars.push_back(FsPar("s3iters",3,"Number of TOTAL iterations to use for MCMC. By default we assign half to burnin and half to sampling. (default: 100000)"));
  pars.push_back(FsPar("s3iterssample",3,"Number of iterations to use for MCMC (default: -ve, meaning derive from s3iters)"));
  pars.push_back(FsPar("s3itersburnin",3,"Number of iterations to use for MCMC burnin (default: -ve, meaning derive from s3iters)"));
  pars.push_back(FsPar("numskip",3,"Number of mcmc iterations per retained sample; (default: -ve, meaning derive from maxretained)"));
  pars.push_back(FsPar("maxretained",3,"Maximum number of samples to retain when numskip -ve. (default: 500)"));
  pars.push_back(FsPar("nummcmcruns",3,"Number of *independent* mcmc runs. (default: 2)"));
  pars.push_back(FsPar("fsmcmcoutput",3,"Filename to use for mcmc output (default: autogenerated)"));
  pars.push_back(FsPar("mcmcGR",3,"Derived. Gelman-Rubin diagnostics obtained from combining MCMC runs, for log-posterior, K,log-beta,delta,f respectively"));
  pars.push_back(FsPar("threshGR",3,"Threshold for the Gelman-Rubin statistic to allow moving on to the tree building stage. We always move on if thresGR<0. (Default: 1.3)"));
  parsize.push_back(pars.size());

  //stage 4 finestructure tree inference
  sectionnames.push_back("4");
  sectioncomments.push_back("FineSTRUCTURE Stage4 tree inference");
  pars.push_back(FsPar("s4args",4,"Extra arguments to the tree building step. (default: none, \"\")"));
  pars.push_back(FsPar("s4iters",4,"Number of maximization steps when finding the best state from which the tree is built. (default: 100000)"));
  pars.push_back(FsPar("fstreeoutput",4,"Filename to use for finestructure tree output. (default: autogenerated)"));
  parsize.push_back(pars.size());

  sectionnames.push_back("5");
  sectioncomments.push_back("Extraction of suitable populations for admixture analysis");
  pars.push_back(FsPar("popidfile",5,"ID file to use for population definition in admixture analysis. Can be autogenerated from stage4, or specified. (Default: autogenerated)"));
  pars.push_back(FsPar("donoridfile",5,"Donor file to use for population definition in admixture analysis. Can be autogenerated to use all populations vs all populations, or specified. (Default: autogenerated)"));
  pars.push_back(FsPar("refdonorfile",5,"Donor file to use for population definition in admixture analysis when processing samples external to this dataset. (Default: autogenerated)"));
  parsize.push_back(pars.size());

  sectionnames.push_back("6");
  sectioncomments.push_back("ChromoPainter population repainting - parameter estimation");
  pars.push_back(FsPar("s6emits",6,"Number of EM iterations (chromopainter -i <n>, default: 10)"));
  pars.push_back(FsPar("s6minsnps",6,"Minimum number of SNPs for EM estimation (for chromopainter -e, default: 10000)"));
  pars.push_back(FsPar("s6snpfrac",6,"fraction of genome to use for EM estimation. (default: 0.1)"));
  pars.push_back(FsPar("s6indfrac",6,"fraction of individuals to use for EM estimation. (default: 1.0)"));
  pars.push_back(FsPar("s67args",6,"arguments to be passed to Chromopainter (default: empty)"));
  pars.push_back(FsPar("s6args",6,"Arguments passed to stage6 (default:-in -iM --emfilesonly)"));
  pars.push_back(FsPar("s6outputroot",6,"Filename to use for stage6 EM estimation. (Default: autogenerated)"));
  pars.push_back(FsPar("popNeinf",6,"Ne parameter inferred from population-level processing. (Default: autogenerated)"));
  pars.push_back(FsPar("popmuinf",6,"mu parameter inferred from population-level processing. (Default: autogenerated)"));
  parsize.push_back(pars.size());

  sectionnames.push_back("7");
  pars.push_back(FsPar("s7combineargs",7,"Additional arguments for stage 7 combine (fs combine; default: none, \"\")"));
  pars.push_back(FsPar("s7chunksperregion",7,"number of chunks in a \"region\" (-ve: use default of 100 for linked, nsnps/100 for unlinked)"));
  pars.push_back(FsPar("s7samples",7,"number of samples of the painting to obtain per recipient haplotype, for examining the details of the painting. (Populates <root>.samples.out; default 0. Warning: these file can get large)"));
  pars.push_back(FsPar("s7args",7,"Arguments passed to stage7 (default: none, \"\")"));
  pars.push_back(FsPar("s7outputroot",7,"Filename to use for stage7 EM estimation. (Default: autogenerated)"));
  pars.push_back(FsPar("popcval",-1,"Derived. 'c' as inferred using chromopainter. This is only used for sanity checking. See s34 args for setting it manually."));
  pars.push_back(FsPar("popcproot",7,"The name of the final population-based chromopainter output. (Default: <filename>_pop, the project file name)"));
  pars.push_back(FsPar("popcpgenomelen",7,"The admixture analysis input file, derived name of the chunklengths file from cproot.",1));
  sectioncomments.push_back("ChromoPainter population repainting (main run)");
  parsize.push_back(pars.size());

  sectionnames.push_back("8");
  sectioncomments.push_back("Mixture estimation");
  parsize.push_back(pars.size());

  sectionnames.push_back("9");
  sectioncomments.push_back("GLOBETROTTER");
  pars.push_back(FsPar("gtroot",9,"The location of the GLOBETROTTER output and content (Default: <root>/gt/gt)"));
  pars.push_back(FsPar("gtbootstrapnum",9,"The number of bootstrap samples in GT (Default: 20)"));
  pars.push_back(FsPar("gtpropscutoff",9,"GT cutoff for small values (Default: 0.001)"));
  pars.push_back(FsPar("gtnummixingiterations",9,"GT number of iterations of maximization (Default: 5)"));
  pars.push_back(FsPar("gtbinwidth",9,"The width of the bins (Default: 0.1)"));
  pars.push_back(FsPar("gtsampleslist",9,"A file which will contain the list of samples needed by GT (Default: <gtroot>_sampleslist.txt)"));
  pars.push_back(FsPar("gtrecomblist",9,"A file which will contain the list of recombination files needed by GT (Default: <gtroot>_recomblist.txt)"));
  parsize.push_back(pars.size());

  // Vector quantities
  sectionnames.push_back("1");
  sectioncomments.push_back("Vector quantities placed at the end for readability");
  pars.push_back(FsPar("phasefiles",1,"Comma or space separated list of all 'phase' files containing the (phased) SNP details for each haplotype. Required. Must be sorted alphanumerically to ensure chromosomes are correctly ordered. So don't use *.phase, use file{1..22}.phase. Override this with upper case -PHASEFILES."));
  pars.push_back(FsPar("recombfiles",1,"Comma or space separated list of all recombination map files containing the recombination distance between SNPs. If provided, a linked analysis is performed. Otherwise an 'unlinked' analysis is performed. Note that linkage is very important for dense markers!"));
  pars.push_back(FsPar("nsnpsvec",-1,"Derived. Comma separated list of the number of SNPs in each phase file.",1));
  pars.push_back(FsPar("s1outputrootvec",-1,"Derived. Comma separated list of the stage 1 output files names.",1));

  parsize.push_back(pars.size());

  sectionnames.push_back("2");
  sectioncomments.push_back("");
  pars.push_back(FsPar("s2outputrootvec",-1,"Derived. Comma separated list of the stage 2 output files names.",1));
  parsize.push_back(pars.size());

  sectionnames.push_back("3");
  sectioncomments.push_back("");
  pars.push_back(FsPar("fsmcmcoutputvec",-1,"Derived. Comma separated list of the stage 3 output files names.",1));
  pars.push_back(FsPar("old_fsmcmcoutputvec",-1,"Derived. Comma separated list of the stage 3 output files names, if we need to continue a too-short MCMC run.",1));

  parsize.push_back(pars.size());

  sectionnames.push_back("4");
  sectioncomments.push_back("");
  pars.push_back(FsPar("fstreeoutputvec",-1,"Derived. Comma separated list of the stage 4 output files names.",1));
  parsize.push_back(pars.size());
  
  sectionnames.push_back("5");
  sectioncomments.push_back("");
  pars.push_back(FsPar("popcmdargs",5,"Derived. Comma separated list of all individual-specific idfiles created in stage5."));
  pars.push_back(FsPar("refidfiles",5,"Derived. Comma separated list of all id files generated for use with external samples. (Default: autogenerated)"));
  parsize.push_back(pars.size());

  sectionnames.push_back("6");
  sectioncomments.push_back("");
  pars.push_back(FsPar("s6outputrootvec",-1,"Derived. Comma separated list of the stage 6 output files names.",1));
  parsize.push_back(pars.size());

  sectionnames.push_back("7");
  sectioncomments.push_back("");
  pars.push_back(FsPar("s7outputrootvec",-1,"Derived. Comma separated list of the stage 7 output files names.",1));
  parsize.push_back(pars.size());

  sectionnames.push_back("8");
  sectioncomments.push_back("");
  parsize.push_back(pars.size());

  sectionnames.push_back("9");
  sectioncomments.push_back("");
  parsize.push_back(pars.size());

  pars.push_back(FsPar("gtdonors",9,"The names of the donor populations used in the admixture analysis (Fixed by donoridfile)"));
  pars.push_back(FsPar("gtrecips",9,"The names of the recipient populations used in the admixture analysis (Fixed by donoridfile)"));
  pars.push_back(FsPar("gtsamplesfiles",9,"A vector of the (combined) per-chromosome samples files. (Fixed by data structure)"));
  pars.push_back(FsPar("gtparamfiles",9,"A vector of the GT parameter file names. (Fixed by data structure)"));

  pars.push_back(FsPar("stage",-1,"Derived. Don't mess with this! The internal measure of which stage of processing we've reached. Change it via -reset or -duplicate.",1));

}

void FsProject::defineCommands()
{
  cmds.clear();
  cmds.push_back(FsCmd("go",-1,"","Do the next things that are necessary to get a complete set of finestructure runs."));
  cmds.push_back(FsCmd("gt",-1,"","Do the preparation of the input files fora GLOBETROTTER analysis."));  
  cmds.push_back(FsCmd("import",-1,"<file>","Import some settings from an external file. If you need to set any non-trivial settings, this is the way to do it. See \"fs -hh\" for more details."));
  ///////////////
  cmds.push_back(FsCmd("createid",1,"<filename>","Create an ID file from a PROVIDED phase file. Individuals are labelled IND1-IND<N>."));

  cmds.push_back(FsCmd("remakes",-1,"<stage>","For hpc mode. Remake the current commandfile, retaining the results of any completed runs."));
  cmds.push_back(FsCmd("duplicate",-1,"<stage> <newfile.cp>","Copy the information from the current settings file, and then -reset it to <stage>."));
  cmds.push_back(FsCmd("reset",-1,"<stage>","Reset the processing to an earlier point. Further \"-go\" commands will rerun any activity from this point. Helpful for rerunning with different parameters."));
  cmds.push_back(FsCmd("configmcmc",3,"<s3itersburnin> <s3iterssample> <numskip>","Shorthand for setting the parameters of the FineSTRUCTURE MCMC. Takes arguments in the form of finestructure's -x -y -z parameters."));
  cmds.push_back(FsCmd("ignoreGR",3,"","Reset the MCMC files to a previously completed but unconverged run, allowing processing to proceed as though it were converged."));
  cmds.push_back(FsCmd("hpcs3",3,"","Generates a special HPC command for the MCMC, performing the Gelman-Rubin check and rerunning the MCMC until convergence as part of a single HPC command."));
  cmds.push_back(FsCmd("cando",-1,"<command>","Checks whether a command can be executed, without executing it."));

  ///////////////
  cmds.push_back(FsCmd("haploid",1,"","Shorthand for the parameter `ploidy:1'"));
  cmds.push_back(FsCmd("countdata",1,"","Ends stage0. Performs checks on the data and confirms that we have valid data."));
  cmds.push_back(FsCmd("countdatapop",5,"","Starts stage5. Performs checks on the -popidfile."));

  for(int i=1;i<=9;i++) {
    stringstream sscmd,sscomment;
    sscmd<<"makes"<<i;
    sscomment<<"Make the stage"<<i<<" commands.";
    cmds.push_back(FsCmd(sscmd.str(),i,"",sscomment.str()));
  }  
  for(int i=1;i<=8;i++) {
    stringstream sscmd,sscomment;
    sscmd<<"dos"<<i;
    sscomment<<"Do the stage"<<i<<" commands. This we should only be doing in single machine mode; we use -writes1 in HPC mode.";
    cmds.push_back(FsCmd(sscmd.str(),i,"",sscomment.str()));
  }  
  for(int i=1;i<=9;i++) {
    stringstream sscmd,sscomment;
    sscmd<<"writes"<<i;
    sscomment<<"Write the stage"<<i<<" commands to file, which we only need in HPC mode.";
    if(i<9) sscomment<<" In single machine mode we can instead use -dos"<<i<<".";
    cmds.push_back(FsCmd(sscmd.str(),i,"",sscomment.str()));
  }  
  cmds.push_back(FsCmd("combines1",1,"","Ends stage1 by combining the output of the stage1 commands. This means estimating the parameters mu and Ne from the output of stage1."));
  cmds.push_back(FsCmd("combines2",2,"","Ends stage2 by combining the output of the stage2 commands. This means estimating 'c' and creating the genome-wide chromopainter output for all individuals."));
  cmds.push_back(FsCmd("combines3",3,"","Ends stage3 by checking the output of the stage3 commands. This means checking that the MCMC has converged."));
  cmds.push_back(FsCmd("combines4",4,"","Ends stage4 by checking the output of the stage4 commands. This means constructing a tree of the best observed state, to form a clustering."));
  cmds.push_back(FsCmd("combines5",5,"","Ends stage5 by making the population IDs from the donor file."));
  cmds.push_back(FsCmd("combines6",6,"","Ends stage6 by combining the output of the stage6 commands. This means estimating the parameters mu and Ne for the population painting from the output of stage6."));
  cmds.push_back(FsCmd("combines7",7,"","Ends stage7 by combining the output of the stage7 commands. This means estimating 'c' and creating the genome-wide chromopainter output for the population painting."));
  cmds.push_back(FsCmd("combines8",8,"",""));

  cmds.push_back(FsCmd("validates",-1,"<stage>","Validate the output files of the stage specified."));
  /*  for(int i=1;i<=9;i++) {
    stringstream sscmd,sscomment;
    sscmd<<"validates"<<i;
    sscomment<<"Validate the output files of the stage specified.";
    cmds.push_back(FsCmd(sscmd.str(),i,"",sscomment.str()));
    } */ 

  ///
  cmds.push_back(FsCmd("paintercmd",-1,"<cmd>","Set the chromopainter version to an external command, which should be ChromoPainterv2. This enables -hpc 1 mode and makes other changes to the processing."));
  //cpname
  //fsmcmcname
}

///////////////////////////////
void FsProject::addHistory(std::vector<std::string> args)
{
  std::ostringstream hist;
  hist<<"CMDLINE: ";
  for(unsigned int c1=0;c1<args.size();++c1){
    hist<<args[c1]<<" ";
  }
  hist<<endl;
  historytext.append(hist.str());

}

void FsProject::switchStdout(const char *newStream)
{
  fflush(stdout);
  fgetpos(stdout, &stdout_pos);
  stdout_fd = dup(fileno(stdout));
  FILE * ftest = freopen(newStream, "w", stdout);
  if(ftest==NULL) throw(runtime_error("Could not switch log file!"));
}

void FsProject::revertStdout()
{
  fflush(stdout);
  dup2(stdout_fd, fileno(stdout));
  close(stdout_fd);
  clearerr(stdout);
  fsetpos(stdout, &stdout_pos);
}

bool FsProject::finishedStage(int teststage)
{
  if(stage>teststage) return(true);
  return(false);
}

int FsProject::getHpc(){
  return(hpc);
}

void FsProject::do_omp_set_num_threads(int numthreads){
#ifdef _OPENMP
  omp_set_num_threads(numthreads);
#endif
}

int FsProject::get_omp_get_thread_num(){
#ifdef _OPENMP
  return(omp_get_thread_num());
#else
  return(0);
#endif
}

int FsProject::get_omp_get_num_threads(){
#ifdef _OPENMP
  return(omp_get_num_threads());
#else
  return(1);
#endif
}

int FsProject::getIndsPerProc(int forstage){
  if(forstage==7) return(1);   // fs2.1: We need to to this individual by individual now! ***
  if(indsperproc>0) return(indsperproc);
  if(hpc) return(1);
  int numshere=numRecipients(stage);
  if(numthreads>0) {
    indsperproc=ceil(((double) numshere)/numthreads);
    return(indsperproc);
  }
  
  // *************** FIXME!

  int tnumthreads = 1,th_id;
  #pragma omp parallel private(th_id)
  {
    th_id = get_omp_get_thread_num();
    if ( th_id == 0 ) {
      tnumthreads  = get_omp_get_num_threads();
    }
  }
  indsperproc=ceil(((double) numshere)/tnumthreads);
  return(indsperproc);
}

int FsProject::defaultChunksperregion(){
  if(linkagemode.compare("linked")==0) return(100);
  else if(linkagemode.compare("unlinked")==0) {
    if(nsnps<0) throw(logic_error("s2chunksperregion"));
    return(max(10,(int)(nsnps/100)));
  }else{
    throw(logic_error("s2chunksperregion"));
  }
}

string FsProject::getCommandFile(int forstage){
  if(forstage<0) forstage=stage;
  switch(forstage){
  case 1: return(s1commandfile);
  case 2: return(s2commandfile);
  case 3: return(s3commandfile);
  case 4: return(s4commandfile);
    //  case 5: return(s5commandfile);
  case 6: return(s6commandfile);
  case 7: return(s7commandfile);
    // case 8: return(s8commandfile);
  case 9: return(s9commandfile);
  }
  return(string(""));
}

int FsProject::getCommandFileCount(int forstage){
  if(forstage<0) forstage=stage;
  switch(forstage){
  case 1: return(s1commands.size());
  case 2: return(s2commands.size());
  case 3: return(s3commands.size());
  case 4: return(s4commands.size());
  case 6: return(s6commands.size());
  case 7: return(s7commands.size());
  case 9: return(s9commands.size());
  }
  return(-1);
}

string FsProject::whichLinkagemode()
{
  if(recombfiles.size()>0) return("linked");
  return("unlinked");
}


//////////////////////////////////
bool FsProject::applyVal(std::vector<std::string> args)
{
  /*for(unsigned int c1=0;c1<args.size();c1++){
    cout<<"DEBUG2 :\""<<args[c1]<<"\""<<endl;
    }*/
  FsSettingsValue val(args);
  //cout<<"DEBUG2 :\""<<val.getName()<<"\" \""<<val.getVal()<<"\""<<endl;
  if(val.success) return(applyVal(val));
  else{
    cerr<<"ERROR: Tried to set a parameter, but this failed. Is this a malformed command line?"<<endl;
    cerr<<"HELP for "<<cmdInfo(args[0],true)<<endl;
    throw(runtime_error("Invalid parameter setting"));
  }
  return(false);
}

bool FsProject::applyVal(FsSettingsValue val)
{
  string name=val.getName();
  //  cout<<"DEBUG :\""<<name<<"\" \""<<val.getVal()<<"\" stage "<<stage<<endl;
  bool found=false;

  //////////// COMMENT 
  if(name.compare("CMD")==0){
    historytext.append(val.getEndedLine()); found=true;
  }else if(name.compare("CMDLINE")==0){
    historytext.append(val.getEndedLine()); found=true;
  }else if(name.compare("stage")==0){
    stage=val.getValAsInt();found=true;
  }

  if(found){
    if(verbose) cout<<"Successfully read "<<val.getName()<<":\""<<val.getVal()<<"\""<<endl;
    return(found);
  }
  //////////// STAGE1

  if(name.compare("exec")==0){
    exec=val.getVal(); found=true;
  }else if(name.compare("fsfile")==0){
    fsfile=val.getVal(); found=true;
  }else if(name.compare("allowdep")==0){
    allowdep=val.getValAsInt(); found=true;
  }else if(name.compare("hpc")==0){
    hpc=val.getValAsInt(); found=true;
  }else if (name.compare("ploidy")==0){
    ploidy=val.getValAsInt(); found=true;
  }else if (name.compare("paintercmd")==0){
    paintercmd=val.getVal(); found=true;
  }else if (name.compare("linkagemode")==0){
    linkagemode=val.getVal(); found=true;
  }else if (name.compare("indsperproc")==0) {
    indsperproc=val.getValAsInt(); found=true;
  }else if (name.compare("outputlogfiles")==0) {
    outputlogfiles=val.getValAsInt(); found=true;
  }else if (name.compare("allowdep")==0) {
    allowdep=val.getValAsInt(); found=true;
  }else if (name.compare("validatedoutput")==0) {
    validatedoutput=val.getValAsIntVec(); found=true;
  }else if (name.compare("missingfiles")==0) {
    missingfiles=val.getValAsIntVec(); found=true;
  }

  // STAGE1-2 generics
  if(name.compare("s12inputtype")==0){
    s12inputtype=val.getVal(); found=true;
  }else if ((name.compare("phasefiles")==0)||(name.compare("PHASEFILES")==0)) {
    vector<string> oldphase=phasefiles;
    //    if(oldphase.size()>0) {cout<<"WARNING: You have specified -phasefiles already, and have added some more. This is often accidental. Check}
    phasefiles=val.getValAsStringVec(); found=true;
    s12inputtype=string("phase");
    vector<string> tphasefiles=phasefiles;
    sort(tphasefiles.begin(),tphasefiles.end());
    if(unique(tphasefiles).size()!=tphasefiles.size()){// duplicates
      cerr<<"ERROR: duplicated phase files provided! You only need to specify each phase file once, and you should not use -phasefiles in future commands."<<endl;
      phasefiles=oldphase;
      throw(runtime_error("duplicate phase files"));
    }
    if(name.compare("PHASEFILES")==0){
      name=string("phasefiles");
    }else if(stage==0 && !isNumericallySorted(phasefiles)) {
      cerr<<"ERROR: Phase files are not lexicographically sorted, so you probably will get confusing assignments of chromosomes to file indices. Rerun with -PHASEFILES instead of -phasefiles to override this check, or rerun with sorted files (e.g. in bash: -phasefiles file{1..22}.phase)."<<endl;
      phasefiles=oldphase;
      throw(runtime_error("unsorted phase files"));
    }
  }else if (name.compare("recombfiles")==0) {
    vector<string> oldrec=recombfiles;
    recombfiles=val.getValAsStringVec(); found=true;
    if(recombfiles.size()>0) linkagemode="linked";
    if(unique(recombfiles).size()!=recombfiles.size()){// duplicates
      recombfiles=oldrec;
      cerr<<"ERROR: duplicated recomb files provided! You only need to specify each recombination file once, and you should not use --recombfiles in future commands."<<endl;
      throw(runtime_error("duplicate recombination files"));
    }
  }else if (name.compare("idfile")==0) {
    idfile=val.getVal(); found=true;
  }else if (name.compare("s12args")==0) {
    s12args=val.getVal(); found=true;
  }else if (name.compare("ninds")==0) {
    ninds=val.getValAsInt(); found=true;
  }else if (name.compare("nindsUsed")==0) {
    nindsUsed=val.getValAsInt(); found=true;
  }else if (name.compare("nrecipsDonor")==0) {
    nrecipsDonor=val.getValAsInt(); found=true;
  }else if (name.compare("nsnps")==0) {
    nsnps=val.getValAsInt(); found=true;
  }else if (name.compare("nsnpsvec")==0) {
    nsnpsvec=val.getValAsIntVec(); found=true;
  }else if (name.compare("numthreads")==0) {
    numthreads=val.getValAsInt(); found=true;
  }
  // STAGE1 properties
  if (name.compare("s1args")==0) {
    s1args=val.getVal(); found=true;
  }else if (name.compare("s1emits")==0) {
    s1emits=val.getValAsInt(); found=true;
  }else if (name.compare("s1outputroot")==0) {
    s1outputroot=val.getVal(); found=true;
  }else if (name.compare("s1outputrootvec")==0) {
    s1outputrootvec=val.getValAsStringVec(); found=true;
  }else if (name.compare("s1minsnps")==0) {
    s1minsnps=val.getValAsInt(); found=true;
  }else if (name.compare("s1snpfrac")==0) {
    s1snpfrac=val.getValAsDouble(); found=true;
  }else if (name.compare("s1indfrac")==0) {
    s1indfrac=val.getValAsDouble(); found=true;
  }

  //////////// STAGE1 POST
  if (name.compare("Neinf")==0) {
    Neinf=val.getValAsDouble(); found=true;
  }else if (name.compare("muinf")==0) {
    muinf=val.getValAsDouble(); found=true;
  }

  if(found){
    return(checkStage(name,val.getVal()));
  }
  //////////// STAGE2 PRE
  if (name.compare("s2chunksperregion")==0) {
    s2chunksperregion=val.getValAsInt(); found=true;
  }else if (name.compare("s2samples")==0) {
    s2samples=val.getValAsInt(); found=true;
  }else if (name.compare("s2args")==0) {
    s2args=val.getVal(); found=true;
  }else if (name.compare("s2outputroot")==0) {
    s2outputroot=val.getVal(); found=true;
  }else if (name.compare("s2outputrootvec")==0) {
    s2outputrootvec=val.getValAsStringVec(); found=true;
  }else if (name.compare("s2combineargs")==0) {
    s2combineargs=val.getVal(); found=true;
  }

  if(found){
    return(checkStage(name,val.getVal()));
  }

  ////////////// STAGE 3
  //stage3-4 generics
  if (name.compare("cval")==0) {
    cval=val.getValAsDouble(); found=true;
  }else if (name.compare("cproot")==0) {
    cproot=val.getVal(); found=true;
    /*    if(cpchunkcounts.compare("")==0) {
      cpchunkcounts=cproot;
      cpchunkcounts.append(".chunkcounts.out");
      }*/
  }else if (name.compare("cpchunkcounts")==0) {
    cpchunkcounts=val.getVal(); found=true;
  }
  //stage 3 parameters
  if (name.compare("s34args")==0) {
    s34args=val.getVal(); found=true;
  }else if (name.compare("fsroot")==0) {
    fsroot=val.getVal(); found=true;
  }
  if (name.compare("s3iters")==0) {
    s3iters=val.getValAsInt(); found=true;
  }else if (name.compare("s3iterssample")==0) {
    s3iterssample=val.getValAsInt(); found=true;
  }else if (name.compare("s3itersburnin")==0) {
    s3itersburnin=val.getValAsInt(); found=true;
  }else if (name.compare("numskip")==0) {
    numskip=val.getValAsInt(); found=true;
  }else if (name.compare("maxretained")==0) {
    maxretained=val.getValAsInt(); found=true;
  }else if (name.compare("nummcmcruns")==0) {
    nummcmcruns=val.getValAsInt(); found=true;
  }else if (name.compare("fsmcmcoutput")==0) {
    fsmcmcoutput=val.getVal(); found=true;
  }else if (name.compare("fsmcmcoutputvec")==0) {
    fsmcmcoutputvec=val.getValAsStringVec(); found=true;
  }else if (name.compare("old_fsmcmcoutputvec")==0) {
    old_fsmcmcoutputvec=val.getValAsStringVec(); found=true;
  }else if (name.compare("mcmcGR")==0) {
    mcmcGR=val.getValAsDoubleVec(); found=true;
  }else if (name.compare("threshGR")==0) {
    threshGR=val.getValAsDouble(); found=true;
  }

  if(found){
    return(checkStage(name,val.getVal()));
  }

  ///////////////// STAGE 4
  if (name.compare("s4iters")==0) {
    s4iters=val.getValAsInt(); found=true;
  }else if (name.compare("s4args")==0) {
    s4args=val.getVal(); found=true;
  }else if (name.compare("fstreeoutput")==0) {
    fstreeoutput=val.getVal(); found=true;
  }else if (name.compare("fstreeoutputvec")==0) {
    fstreeoutputvec=val.getValAsStringVec(); found=true;
  }
  if(found){
    return(checkStage(name,val.getVal()));
  }

  ///////////////// STAGE 4combine
  if (name.compare("s4iters")==0) {
    s4iters=val.getValAsInt(); found=true;
  }else if (name.compare("fstreeoutput")==0) {
    fstreeoutput=val.getVal(); found=true;
  }else if (name.compare("fstreeoutputvec")==0) {
    fstreeoutputvec=val.getValAsStringVec(); found=true;
  }
  if(found){
    return(checkStage(name,val.getVal()));
  }

    ///////////////// STAGE 5
  if (name.compare("popidfile")==0) {
    popidfile=val.getVal(); found=true;
  }else if (name.compare("donoridfile")==0) {
    donoridfile=val.getVal(); found=true;
  }else if (name.compare("popcmdargs")==0) {
    popcmdargs=val.getValAsStringVec(','); found=true;
  }else if (name.compare("refdonorfile")==0) {
    refdonorfile=val.getVal(); found=true;
  }else if (name.compare("refidfiles")==0) {
    refidfiles=val.getValAsStringVec(','); found=true;
  }
  if(found){
    return(checkStage(name,val.getVal()));
  }

  ///////////////// STAGE 6
  if (name.compare("s6args")==0) {
    s6args=val.getVal(); found=true;
  }else if (name.compare("s67args")==0) {
    s67args=val.getVal(); found=true;
  }else if (name.compare("s6emits")==0) {
    s6emits=val.getValAsInt(); found=true;
  }else if (name.compare("s6outputroot")==0) {
    s6outputroot=val.getVal(); found=true;
  }else if (name.compare("popNeinf")==0) {
    popNeinf=val.getValAsDouble(); found=true;
  }else if (name.compare("popmuinf")==0) {
    popmuinf=val.getValAsDouble(); found=true;
  }else if (name.compare("s6outputrootvec")==0) {
    s6outputrootvec=val.getValAsStringVec(); found=true;
  }else if (name.compare("s6minsnps")==0) {
    s6minsnps=val.getValAsInt(); found=true;
  }else if (name.compare("s6snpfrac")==0) {
    s6snpfrac=val.getValAsDouble(); found=true;
  }else if (name.compare("s6indfrac")==0) {
    s6indfrac=val.getValAsDouble(); found=true;
  }
  if(found){
    return(checkStage(name,val.getVal()));
  }

  ///////////////// STAGE 7
  if (name.compare("s7chunksperregion")==0) {
    s7chunksperregion=val.getValAsInt(); found=true;
  }else if (name.compare("s7samples")==0) {
    s7samples=val.getValAsInt(); found=true;
  }else if (name.compare("s7args")==0) {
    s7args=val.getVal(); found=true;
  }else if (name.compare("popcval")==0) {
    popcval=val.getValAsDouble(); found=true;
  }else if (name.compare("popcproot")==0) {
    popcproot=val.getVal(); found=true;
  }else if (name.compare("popcpgenomelen")==0) {
    popcpgenomelen=val.getVal(); found=true;
  }else  if (name.compare("s7outputroot")==0) {
    s7outputroot=val.getVal(); found=true;
  }else if (name.compare("s7outputrootvec")==0) {
    s7outputrootvec=val.getValAsStringVec(); found=true;
  }else if (name.compare("s7combineargs")==0) {
    s7combineargs=val.getVal(); found=true;
  }

  ///////////////// STAGE 9 (GT)
  if (name.compare("gtroot")==0) {
    gtroot=val.getVal(); found=true;
  }else if(name.compare("gtsampleslist")==0) {
    gtsampleslist=val.getVal(); found=true;
  }else if(name.compare("gtbootstrapnum")==0) {
    gtbootstrapnum=val.getVal(); found=true;
  }else if(name.compare("gtpropscutoff")==0) {
    gtpropscutoff=val.getVal(); found=true;
  }else if(name.compare("gtbinwidth")==0) {
    gtbinwidth=val.getVal(); found=true;
  }else if(name.compare("gtnummixingiterations")==0) {
    gtnummixingiterations=val.getVal(); found=true;
  }
  
  if(found){
    return(checkStage(name,val.getVal()));
  }

  cerr<<"WARNING: Setting "<<val.getName()<<" not found!"<<endl;
  return(found);

}

void FsProject::resetToStage(int newstage)
{
  if(newstage>stage){
    cerr<<"Attempted to reset to stage "<<newstage<<" but we are only in stage "<<stage<<". Can only reset to the current or past stages!"<<endl;
    throw(runtime_error("reset stage error"));
  }
  for(int i=newstage;i<nstages;i++) {
    validatedoutput[i]=0;
    missingfiles[i]=-1;
  }
  if(newstage<=0){ //reset the  stage0 stuff
    if(verbose) cout<<"Removing phase/recomb input information..."<<endl;
    phasefiles.clear();
    recombfiles.clear();
    nsnpsvec.clear();
    nsnps=0;
    ninds=0;
    nindsUsed=0;
    deleteStage("0");
  }
  if(newstage<=1){ //reset the stage1 stuff
    if(verbose) cout<<"Removing EM results..."<<endl;
    s1commands.clear();
    s1logfiles.clear();
    s1outputrootvec.clear();
    s1outputroot="";
    Neinf=-1; 
    muinf=-1;
    safeRemoveCommandFile(1);
    deleteStage("1");
  }
  if(newstage<=2){
    if(verbose) cout<<"Removing chromopainter main run results..."<<endl;
    s2commands.clear();
    s2logfiles.clear();
    s2outputrootvec.clear();
    s2outputroot="";
    s2combineargs="";
    cval=-1;
    cproot="";
    cpchunkcounts="";
    safeRemoveCommandFile(2);
    deleteStage("2");
    deleteStage("2a");
  }
  if(newstage<=3){
    if(verbose) cout<<"Removing finestructure mcmc run results..."<<endl;
    s3commands.clear();
    s3logfiles.clear();
    fsroot="";
    fsmcmcoutput="";
    fsmcmcoutputvec.clear();
    old_fsmcmcoutputvec.clear();
    mcmcGR.clear();
    safeRemoveCommandFile(3);
    deleteStage("3");
  }
  if(newstage<=4){
    if(verbose) cout<<"Removing finestructure tree run results..."<<endl;
    s4commands.clear();
    s4logfiles.clear();
    fsroot="";
    s4args="";
    fstreeoutput="";
    fstreeoutputvec.clear();
    safeRemoveCommandFile(4);
    deleteStage("4");
  }
  if(newstage<=5){
    if(verbose) cout<<"Removing admixture panel construction results..."<<endl;
    popidfile="";
    donoridfile="";
    refdonorfile="";
    refidfiles.clear();
    popcmdargs.clear();
    nrecipsDonor=0;
    deleteStage("5");
  }  
  if(newstage<=6){
    if(verbose) cout<<"Removing admixture EM results..."<<endl;
    s6commands.clear();
    s6logfiles.clear();
    s6outputrootvec.clear();
    s6outputroot="";
    popNeinf=-1; 
    popmuinf=-1;
    safeRemoveCommandFile(6);
    deleteStage("6");
  }
  if(newstage<=7){
    if(verbose) cout<<"Removing chromopainter admixture run results..."<<endl;
    s7combineargs="";
    s7commands.clear();
    s7logfiles.clear();
    s7outputrootvec.clear();
    s7outputroot="";
    popcval=-1;
    popcproot="";
    popcpgenomelen="";
    safeRemoveCommandFile(7);
    deleteStage("7");
    deleteStage("7a");
  }
  if(newstage<=9){
    if(verbose) cout<<"Removing globetrotter run results..."<<endl;
    gtroot="";
    gtsampleslist="";
    gtrecomblist="";
    gtrecips.clear();
    gtdonors.clear();
    gtsamplesfiles.clear();
    gtparamfiles.clear();
    safeRemoveCommandFile(9);
    //    deleteStage("gt");
  }
  stage=newstage;
}

void FsProject::createDuplicated(int newstage,string newname)
{
  dirname=projectroot(newname);
  newname=projectfull(dirname);
  setFileName(newname);
  dirname=projectroot(newname);
  if((newstage>2) & (newstage<5)){ // keeping the chromopainter results, redoing fs
    ensureCpRunARoot(2);
  }else{
    if(applyFiles("new",string("stop"))!=0) {
      throw(runtime_error("duplicate error"));
    }
  }
  ///////////*****************
  // Need to decide what we do when duplicating stages 6-8
}

void FsProject::readFromFile()
{
  if(filename.compare("")==0) {
    throw(logic_error("fsproject: attempting to read project before file name is provided!"));
  }
  FsSettingsReader infile(filename,verbose);
  FsSettingsValue val;
  bool done=false;
  while(!done){
    val=infile.getNext();
    //    cout<<"DEBUG2 :\""<<val.getName()<<"\" \""<<val.getVal()<<"\""<<endl;
    if(val.success){
      applyVal(val);
    }else{
      done=true;
    }
  }
}

string FsProject::writeToFile(string tpar)
{
  ostringstream os;
      if(tpar.compare("exec")==0)os<<exec;
      else if(tpar.compare("allowdep")==0)  os<<allowdep;
      else if(tpar.compare("hpc")==0)os<<hpc;
      else if(tpar.compare("numthreads")==0) os<<numthreads;
      else if(tpar.compare("ploidy")==0)  os<<ploidy;
      else if(tpar.compare("paintercmd")==0)  os<<paintercmd;
      else if(tpar.compare("linkagemode")==0)  os<<linkagemode;
      else if(tpar.compare("indsperproc")==0)  os<<indsperproc;
      else if(tpar.compare("outputlogfiles")==0)  os<<outputlogfiles;
      else if(tpar.compare("allowdep")==0)  os<<allowdep;
      else if(tpar.compare("s12inputtype")==0)  os<<s12inputtype;
      else if(tpar.compare("idfile")==0)  os<<idfile;
      else if(tpar.compare("s12args")==0)  os<<s12args;
      else if(tpar.compare("ninds")==0)  os<<ninds;
      else if(tpar.compare("nindsUsed")==0)  os<<nindsUsed;
      else if(tpar.compare("nrecipsDonor")==0)  os<<nrecipsDonor;
      else if(tpar.compare("nsnps")==0)  os<<nsnps;
      else if(tpar.compare("s1args")==0)  os<<s1args;
      else if(tpar.compare("s1emits")==0)  os<<s1emits;
      else if(tpar.compare("s1minsnps")==0)  os<<s1minsnps;
      else if(tpar.compare("s1snpfrac")==0)  os<<s1snpfrac;
      else if(tpar.compare("s1indfrac")==0)  os<<s1indfrac;
      else if(tpar.compare("s1outputroot")==0)  os<<s1outputroot;
      else if(tpar.compare("Neinf")==0)  os<<Neinf;
      else if(tpar.compare("muinf")==0)  os<<muinf;
      else if(tpar.compare("s2chunksperregion")==0)  os<<s2chunksperregion;
      else if(tpar.compare("s2samples")==0)  os<<s2samples;
      else if(tpar.compare("s2args")==0)  os<<s2args;
      else if(tpar.compare("s2outputroot")==0)  os<<s2outputroot;
      else if(tpar.compare("s2combineargs")==0)  os<<s2combineargs;
      else if(tpar.compare("cval")==0)  os<<cval;
      else if(tpar.compare("cproot")==0)  os<<cproot;
      else if(tpar.compare("cpchunkcounts")==0)  os<<cpchunkcounts;
      else if(tpar.compare("s34args")==0)  os<<s34args;
      else if(tpar.compare("fsroot")==0)  os<<fsroot;
      else if(tpar.compare("s3iters")==0)  os<<s3iters;
      else if(tpar.compare("s3iterssample")==0)  os<<s3iterssample;
      else if(tpar.compare("s3itersburnin")==0)  os<<s3itersburnin;
      else if(tpar.compare("numskip")==0)  os<<numskip;
      else if(tpar.compare("maxretained")==0)  os<<maxretained;
      else if(tpar.compare("nummcmcruns")==0)  os<<nummcmcruns;
      else if(tpar.compare("fsmcmcoutput")==0)  os<<fsmcmcoutput;
      else if(tpar.compare("mcmcGR")==0)  os<<ssvec(mcmcGR);
      else if(tpar.compare("threshGR")==0)  os<<threshGR;
      else if(tpar.compare("s4iters")==0)  os<<s4iters;
      else if(tpar.compare("s4args")==0)  os<<s4args;
      else if(tpar.compare("fstreeoutput")==0)  os<<fstreeoutput;
      
      else if(tpar.compare("popidfile")==0)  os<<popidfile;
      else if(tpar.compare("donoridfile")==0)  os<<donoridfile;
      else if(tpar.compare("refdonorfile")==0)  os<<refdonorfile;
      else if(tpar.compare("s67args")==0)  os<<s67args;
      else if(tpar.compare("s6args")==0)  os<<s6args;
      else if(tpar.compare("s6emits")==0)  os<<s6emits;
      else if(tpar.compare("s6minsnps")==0)  os<<s6minsnps;
      else if(tpar.compare("s6snpfrac")==0)  os<<s6snpfrac;
      else if(tpar.compare("s6indfrac")==0)  os<<s6indfrac;
      else if(tpar.compare("s6outputroot")==0)  os<<s6outputroot;
      else if(tpar.compare("popNeinf")==0)  os<<popNeinf;
      else if(tpar.compare("popmuinf")==0)  os<<popmuinf;

      else if(tpar.compare("s7combineargs")==0)  os<<s7combineargs;
      else if(tpar.compare("s7chunksperregion")==0)  os<<s7chunksperregion;
      else if(tpar.compare("s7samples")==0)  os<<s7samples;
      else if(tpar.compare("s7args")==0)  os<<s7args;
      else if(tpar.compare("s7outputroot")==0)  os<<s7outputroot;
      else if(tpar.compare("popcval")==0)  os<<popcval;
      else if(tpar.compare("popcproot")==0)  os<<popcproot;
      else if(tpar.compare("popcpgenomelen")==0)  os<<popcpgenomelen;

      else if(tpar.compare("gtroot")==0)  os<<gtroot;
      else if(tpar.compare("gtsampleslist")==0)  os<<gtsampleslist;
      else if(tpar.compare("gtrecomblist")==0)  os<<gtrecomblist;
      else if(tpar.compare("gtbootstrapnum")==0)  os<<gtbootstrapnum;
      else if(tpar.compare("gtpropscutoff")==0)  os<<gtpropscutoff;
      else if(tpar.compare("gtbinwidth")==0)  os<<gtbinwidth;
      else if(tpar.compare("gtnummixingiterations")==0)  os<<gtnummixingiterations;
      
      else if(tpar.compare("phasefiles")==0)  os<<ssvec(phasefiles);
      else if(tpar.compare("recombfiles")==0)  os<<ssvec(recombfiles);
      else if(tpar.compare("nsnpsvec")==0)  os<<ssvec(nsnpsvec);
      else if(tpar.compare("s1outputrootvec")==0)  os<<ssvec(s1outputrootvec);
      else if(tpar.compare("s2outputrootvec")==0)  os<<ssvec(s2outputrootvec);
      else if(tpar.compare("s6outputrootvec")==0)  os<<ssvec(s6outputrootvec);
      else if(tpar.compare("s7outputrootvec")==0)  os<<ssvec(s7outputrootvec);
      else if(tpar.compare("fsmcmcoutputvec")==0)  os<<ssvec(fsmcmcoutputvec);
      else if(tpar.compare("old_fsmcmcoutputvec")==0)  os<<ssvec(old_fsmcmcoutputvec);
      else if(tpar.compare("fstreeoutputvec")==0)  os<<ssvec(fstreeoutputvec);
      else if(tpar.compare("validatedoutput")==0)  os<<ssvec(validatedoutput);
      else if(tpar.compare("missingfiles")==0)  os<<ssvec(missingfiles);
      else if(tpar.compare("popcmdargs")==0)  os<<ssvec(popcmdargs);
      else if(tpar.compare("refidfiles")==0)  os<<ssvec(refidfiles);
      else if(tpar.compare("gtdonors")==0)  os<<ssvec(gtdonors);
      else if(tpar.compare("gtrecips")==0)  os<<ssvec(gtrecips);
      else if(tpar.compare("gtsamplesfiles")==0)  os<<ssvec(gtsamplesfiles);
      else if(tpar.compare("gtparamfiles")==0)  os<<ssvec(gtparamfiles);
      else {
	stringstream ss;
	ss<<"Trying to write a command ("<<tpar<<")"<<" that doesn't exist!";
	  throw(runtime_error("fsproject: parameter doesn't exist!"));
      }
      return(os.str());
}

void FsProject::writeToFile()
{
  if(filename.compare("")==0) {
    throw(logic_error("fsproject: attempting to write project before file name is provided!"));
  }
  filebuf fb;
  try{
    fb.open (filename.c_str(),ios::out);
  }catch(exception &x){
    cerr<<"Error opening file!"<<endl<<x.what()<<endl; 
    throw(runtime_error("fsproject: cannot write to file!"));
  }
  ostream os (&fb);
  //////////////////////////////////////////////////
  os<<"section:fsproject"<<endl;
  os<<freetextcomments[0];
  os<<"section:history"<<endl;
  os<<historytext<<endl;
  os<<"section:parameters"<<endl;

  int paron=0;
  string lastsection="",thissection;
  for(int secon=0;secon<(int)sectionnames.size();++secon) {
    thissection=sectionnames[secon];
 if(thissection.compare(lastsection)!=0) {
      os<<endl<<"###################"<<endl<<"stage:"<<thissection<<endl;
      lastsection=thissection;
    }
    os<<"## "<<sectioncomments[secon]<<endl<<endl;
    while(paron<parsize[secon]){
      string tpar=pars[paron].getName();
      os<<tpar<<":";
      try{
	os<<writeToFile(tpar);
      }catch(exception &x){
	// Parameter doesn't exist
	throw(x);
      }
      os<<"  # "<<pars[paron].getHelp()<<endl;
      ++paron;
    }
  }
  os<<endl<<"section:fsprojectend"<<endl;
  os<<"stage:"<<stage<<"  # "<<pars[paron].getHelp()<<endl;
  fb.close();
}
  
void FsProject::setupSections()
{
  string comments=constcommentspre;
  comments.append(constname0);
  comments.append(constname1);
  comments.append(constname2);
  comments.append(constname3);
  comments.append(constname4);
  comments.append(constname5);
  comments.append(constname6);
  comments.append(constname7);
  comments.append(constname8);
  comments.append(constcommentspost);
  freetextcomments.resize(0);
  freetextcomments.push_back(comments);
  freetextcomments.push_back(constname0);
  freetextcomments.push_back(constname1);
  freetextcomments.push_back(constname2);
  freetextcomments.push_back(constname3);
  freetextcomments.push_back(constname4);
  freetextcomments.push_back(constname5);
  freetextcomments.push_back(constname6);
  freetextcomments.push_back(constname7);
  freetextcomments.push_back(constname8);
}

void FsProject::setFileName(std::string f)
{
  filename=f;
  fileroot=f.substr(0,f.find_last_of('.'));
}

std::string FsProject::getFileName()
{
  return(filename);
}

std::string FsProject::getExec()
{
  return(exec);
}

void FsProject::setDirectoryName(std::string d)
{
  dirname=d;
}

std::string FsProject::getDirectoryName()
{
  return(dirname);
}

void FsProject::setVerbose(bool verbose) 
{
  this->verbose=verbose;
}

void FsProject::setFsmode(int fsmode) 
{
  this->fsmode=fsmode;
}

void FsProject::countSnps(){
  // Counts the number of SNPs in the recombination and phase files to check that the match

  // Count from the phase files
  nsnpsvec.clear();
  nsnps=0;
  for(unsigned int c1=0;c1<phasefiles.size();c1++){
    int tnsnps=getNsnpsFromFile(c1);
    nsnpsvec.push_back(tnsnps);
    nsnps+=tnsnps;
  }
  if(recombfiles.size()==0){
    cout<<"Not using recombination map, and therefore using unlinked mode."<<endl;
  }else{
    if(unique(recombfiles).size()!=recombfiles.size()){// duplicates
      cerr<<"ERROR: duplicated recomb files provided! You only need to specify each recombination file once, and you should not use --recombfiles in future commands."<<endl;
      throw(runtime_error("fsproject: duplicate recombination files"));
    }
    if(phasefiles.size()!=recombfiles.size()) throw(runtime_error("fsproject: Different number of phase to recombination files!"));
    for(unsigned int c1=0;c1<recombfiles.size();++c1){
      if(countLines(recombfiles[c1])!=nsnpsvec[c1]+1){
	cerr<<"SNP difference between phase file ("<<phasefiles[c1]<<") with "<<nsnpsvec[c1]<<" SNPs and recombination map ("<<recombfiles[c1]<<") with "<<countLines(recombfiles[c1])-1<<" SNPs. Did you specify the right recombination map? Do you need to recreate the map with convertrecfile.pl? Did you include the header line? Did you check for whether each line has a newline character?"<<endl;
	throw(runtime_error("fsproject: SNP conflict between phase and recombination files!"));
      }
    }
    
  }
  
  if(verbose){
    cout<<"Using "<<nsnps<<" total SNPs over "<<phasefiles.size()<<" phase files."<<endl;
  }
}

void FsProject::countData(){
  // Counts the data and checks that we have sane 
  if(verbose) cout<<"Counting data"<<endl;
  ninds=getNindsFromFile(true);
  int nindsUnique=getUniqueNindsFromFile(true);
  if(ninds!=nindsUnique) {
  cerr<<"ERROR: ID file "<<idfile<<" appears to contain multiple IDs of the same name. Recreate this file without duplicates, or test with -createid <filename> to create a valid (uninformative) ID file."<<endl;
      throw(runtime_error("data error"));
  } 
  nindsUsed=getNindsFromFile(false);
  for(unsigned int c1=0;c1<phasefiles.size();c1++) {
    if(ninds==0 || getNhapsFromFile(c1)==0){
      cerr<<"ERROR: Failed to process a file. Found "<<ninds<<" individuals in the ID file called "<<idfile<<", and in phase file called "<<phasefiles[c1]<<" with "<<getNhapsFromFile(c1)<<" inds. This is a data problem."<<endl;
      throw(runtime_error("data error"));
    }
    int thaps=getNhapsFromFile(c1);
    if(thaps!=ninds*ploidy){
      cerr<<"ERROR: Mismatch between number of individuals in the ID file called "<<idfile<<" with "<<ninds<<" inds, and in phase file called "<<phasefiles[c1]<<" with "<<getNhapsFromFile(c1)<<" haplotypes. This can be a data or a ploidy problem."<<endl;
      double tploidyrem=thaps % ninds;
      double tploidy=thaps / ninds;
      if(tploidyrem==0) {
         cerr<<"INFORMATION: You appear to have "<<tploidy<<" haplotype";
	 if(tploidy>1)cerr<<"s";
	 cerr<<" per individual";
	 if(tploidy==1) cerr<<", i.e. haploid data";
	 cerr<<". Suggest rerunning with \"-ploidy "<<tploidy<<"\"";
	 if(tploidy==2) cerr<<" (or simply omitting the ploidy flag)";
	 cerr<<"."<<endl;
      }
      throw(runtime_error("data error"));
    }
  }
  if(verbose){
    cout<<"Counted "<<ninds<<" Individuals making up "<<ninds*ploidy<<" haplotypes, using "<<nindsUsed<<" of them."<<endl;
  }
  countSnps();
  stage=1;
  validatedoutput[0]=1;
}

void FsProject::countDataPop(){
  // This is where we check that we have a valid popidfile for comparison t0 the idfile. We also run -countdata if it has not already been run.
  if(verbose) cout<<"Counting data from popidfile"<<endl;
  if(stage<1 && validatedoutput[0]==0) countData();

  // Checks from here
  std::vector<std::string> tnames=getIdsFromFile(idfile,true);
  std::vector<std::string> pnames=getIdsFromFile(popidfile,true);
  if(tnames.size() != pnames.size()) {
    cerr<<"ERROR: -idfile has "<<tnames.size()<<" total individuals, whilst -popidfile has "<<pnames.size()<<".  These should refer to the same set of individuals."<<endl;
    throw(runtime_error("countDataPop id mismatch"));
  }
  std::vector<std::string> tnames2=getIdsFromFile(idfile,false);
  std::vector<std::string> pnames2=getIdsFromFile(popidfile,false);
  if(tnames2.size() != pnames2.size()) {
    cerr<<"ERROR: -idfile has "<<tnames2.size()<<" retained individuals, whilst -popidfile has "<<pnames2.size()<<".  These should refer to the same set of individuals."<<endl;
    throw(runtime_error("countDataPop id mismatch"));
  }
  for(unsigned int c1=0;c1<tnames.size();++c1){
    if(tnames[c1]!=pnames[c1]){
      cerr<<"ERROR: -idfile has name "<<tnames[c1]<<" at row "<<c1<<" wheras -popidfile has name "<<pnames[c1]<<".  These should refer to the same set of individuals."<<endl;
      throw(runtime_error("countDataPop id mismatch"));
    }
  }
  countSnps();
  //
  stage=5;
}

void FsProject::copyfile(std::string from, std::string to)
{
  std::ifstream  src(from.c_str(), std::ios::binary);
  std::ofstream  dst(to.c_str(),   std::ios::binary);
  dst << src.rdbuf();
}

void FsProject::deleteStage(string ending){
  stringstream ss;
  ss<<dirname<<"/stage"<<ending<<"/";
  if(directoryExists(ss.str())) {
    deleteFolderTree(ss.str());
  }
}

void FsProject::safeRemoveCommandFile(int forstage){
  ensureCommandfile(forstage);
  switch(forstage){
  case 1: backupFileAndRemove(s1commandfile,"commandfiles");s1commandfile="";break;
  case 2: backupFileAndRemove(s2commandfile,"commandfiles");s2commandfile="";break;
  case 3: backupFileAndRemove(s3commandfile,"commandfiles");s3commandfile="";break;
  case 4: backupFileAndRemove(s4commandfile,"commandfiles");s4commandfile="";break;
  case 6: backupFileAndRemove(s6commandfile,"commandfiles");s6commandfile="";break;
  case 7: backupFileAndRemove(s7commandfile,"commandfiles");s7commandfile="";break;
  case 9: backupFileAndRemove(s9commandfile,"commandfiles");s9commandfile="";break;
  default:break;
  }
}

void FsProject::backupFileAndRemove(string filename,string backupdir){
  if( access( filename.c_str(), F_OK ) != -1 ) { // file exists
    if(verbose) cout<<"FsProject: backing up and removing file: "<<filename<<endl;
    backupFile(filename);
    unlink(filename.c_str());
  }
}

void FsProject::backupFile(string filename,string backupdir){

  if( access( filename.c_str(), F_OK ) != -1 ) { // file exists
    if(verbose) cout<<"File: "<<filename<<" already exists; backing it up..."<<endl;

    std::ostringstream backupfilenamess;
    int i=0;
    do{
      i++; // first backup is called <x>1.bak
      backupfilenamess.str("");
      backupfilenamess << dirname<<"/"<<backupdir<<"/";
      ensureDirectory(backupfilenamess.str());
      backupfilenamess<<baseName(filename) << i << ".bak";
      if(i>10000) {
	throw(logic_error("FsProject::safeCreateFile concern: found > 10000 backups, can't find a safe name to rename old project file?"));
      }
    }while(access( backupfilenamess.str().c_str(), F_OK ) != -1 ); // backupfile exists

    if(verbose) cout<<"Backing up "<<filename<<" to "<<backupfilenamess.str()<<endl;
    copyfile(filename,backupfilenamess.str());
  } 
  
}

void FsProject::safeCreateFile()
{
  backupFile(filename);
  writeToFile();
}

void FsProject::writeStringToFile(std::string s,std::string fn) 
{
  filebuf fb;
  try{
    fb.open (fn.c_str(),ios::out);
  }catch(exception &x){
    cerr<<"Error opening file!"<<endl<<x.what()<<endl; 
    throw(runtime_error("fsproject: cannot write to file!"));
  }
  ostream os (&fb);
  os<<s;
  fb.close();
}

void FsProject::writeStringVectorToFile(std::vector<std::string> sv,std::string fn,std::vector<std::string> logfiles,bool addfs) 
{
  filebuf fb;
  try{
    fb.open (fn.c_str(),ios::out);
  }catch(exception &x){
    cerr<<"Error opening file!"<<endl<<x.what()<<endl; 
    throw(runtime_error("fsproject: cannot write to file!"));
  }
  ostream os (&fb);
  for(unsigned int i=0;i<sv.size();i++){
    if(addfs) os<<exec<<" ";
    os<<sv[i];
    if(logfiles.size()==sv.size()) os<<" > "<<logfiles[i]<<" 2>&1";
    os<<endl;
  }
  fb.close();
}

std::vector<char *> FsProject::converttoargv(std::string s)
{
  std::string::iterator new_end = std::unique(s.begin(), s.end(), BothAreSpaces);
  s.erase(new_end, s.end()); 
  std::vector<std::string> mycmd=split(s,' ');
  std::vector<char *> ccmd(mycmd.size());
  for (unsigned int j = 0; j < mycmd.size(); ++j) {
    ccmd[j] = (char *)malloc(mycmd[j].size()+1);
    strcpy(ccmd[j],mycmd[j].c_str());
  }
  return(ccmd);
}

void FsProject::freeargv(std::vector<char *> ccmd)
{
  for (unsigned int j = 0; j < ccmd.size(); ++j) {
    free(ccmd[j]);
  }
}

void FsProject::doCpStage(int stage) {
#ifndef _OPENMP
  printf("WARNING: You have compiled this code without OpenMP support. Parallel processing is not possible without using HPC mode (\"-hpc 1\") and external parallelization. Consider recompiling after reconfiguring openMP.\n");
#endif
  if(fsmode==0) {
    if(stage==6){
      cout<<"Performing PARAMETER ESTIMATION, which is followed by CHROMOSOME PAINTING. They should take similar amounts of time using the default settings."<<endl;
    }else{
      cout<<"Performing CHROMOSOME PAINTING."<<endl;
    }
  }
  std::vector<std::string> sv;
  if(stage==1) {sv=s1commands;
  }else if(stage==2) {sv=s2commands;
  }else if(stage==6) {sv=s6commands;
  }else if(stage==7) {sv=s7commands;
  }else{
    throw(logic_error("fsproject: calling an invalid chromopainter stage (should be 1-2 for finestructure and 6-7 for admixture)"));
  }

  int cmdon=0;
  if(numthreads>0) {
    do_omp_set_num_threads(numthreads);
  }
  int allok=1;
#pragma omp parallel for
  for(unsigned int i=0;i<sv.size();i++)  {
    //    int thread_number = omp_get_thread_num();
    //    if(!allok)continue;
    string tsv=sv[i];
    string logfile;
    std::ostringstream ss;
    if(stage==1) logfile=s1outputrootvec[i];
    else if(stage==2) logfile=s2outputrootvec[i];
    else if(stage==6) logfile=s6outputrootvec[i];
    else if(stage==7) logfile=s7outputrootvec[i];
    logfile.append(".log");
    tsv.append(" --noexitonerrors ");
    tsv.append(logfile);

    std::vector<char *> argv=converttoargv(tsv);
    
#pragma omp atomic
    cmdon++;
#pragma omp critical
    {
    if(fsmode>0) {
      ss<<"Running stage "<<stage;
      if((stage==1) | (stage==6)) ss<<" (chromopainter parameter estimation) ";
      else ss<<" (chromopainter painting) ";
      if(stage>=6) ss<<"(admixture) ";
      ss<<"command number (~"<<cmdon<<") of "<<sv.size();
      ss<<" (logging to "<<logfile<<")";
      ss<<"\n";
    }else{
      ss<<"Running";
      if((stage==1) | (stage==6)) ss<<" chromopainter parameter estimation, ";
      else ss<<" chromopainter painting, ";
      ss<<"command number (~"<<cmdon<<") of "<<sv.size();
      ss<<"\n";
    }
    cout<<(ParallelStream()<<ss.str()).toString();
    if(verbose) cout<<(ParallelStream()<<"RUNNING S"<<stage<<" CMD:"<<tsv<<"\n").toString();
    }
    
    //switchStdout(logfile.c_str());
    int rv=chromopainter(argv.size(),argv.data());
    freeargv(argv);
    //revertStdout();
    
    // check that it ran correctly
    if((getLastLine(logfile).compare(cpsuccesstext)!=0)|| (rv>0)) {
#pragma omp critical
      {
      cerr<<"ChromoPainter Run "<<i<<" failed! Return value was "<<rv<<endl;
      string el=getLineContaining("Exiting",logfile);
      string ll=getLastLine(logfile);
      cerr<<"Error line was \""<<el<<"\""<<endl;
      if(ll.compare("")==0) {
	cerr<<"Logging failed. This usually means that chromopainter terminated abnormally."<<endl;
      }else{
	cerr<<"See log file ("<<logfile<<") for more details."<<endl;
      }
      allok=0;
      }
    }// endif
  }// end for
  if(!allok && ((stage==2)| (stage==7))) throw(runtime_error("chromopainter"));
}

void FsProject::addEMline(string line, vector<double> *Nevec_ptr, vector<double> *muvec_ptr){
  // add an observation of parameters to the vectosr Nevec_ptr and muvec_ptr
  std::vector<std::string> sline=split(line,' ');
  if(sline.size()<(unsigned int)(3+ploidy)){
    throw(runtime_error("combineEM file error"));
  }
  Nevec_ptr->push_back(atof(sline[1+ploidy].c_str()));
  muvec_ptr->push_back(atof(sline[2+ploidy].c_str()));
  if(verbose) cout<<"... Found EM values: Ne="<<sline[1+ploidy].c_str()<<" and mu="<<sline[2+ploidy].c_str()<<endl;
}

void FsProject::addEMobs(string filename, vector<double> *Nevec_ptr, vector<double> *muvec_ptr){
  // add the last line of an EM file to be the observation of the parameters
  ifstream file(filename.c_str());
  string lastline(""),line;
  int foundnew;
  while(getline(file, line)){
    foundnew=0;
    if(line.substr(0,5).compare("EMPAR")==0){
      foundnew=1;
    }
    if(foundnew==1 && lastline.compare("")!=0) addEMline(lastline,Nevec_ptr,muvec_ptr);
    lastline=line;
  }
  if(lastline.compare("")!=0) {
      addEMline(lastline,Nevec_ptr,muvec_ptr);
  }
  file.close();

}

void FsProject::combineCpEm(int fromstage)
{ // extract all EMPARAMETERS from stage1 files
  vector<double> Nevec;
  vector<double> muvec;
  vector<string> thisstageemfiles;
  vector<string> soutputrootvec;
  if(fromstage==1){
    soutputrootvec=s1outputrootvec;
  }else{
    soutputrootvec=s6outputrootvec;
  }
  // sanity checks
  if(soutputrootvec.size()<1){ // invalid; we need output files to combine
    cerr<<"Combining ChromoPainter EM files for parameter inference requires stage"<<fromstage<<" to have been defined and run."<<endl;
    stringstream ss;
    ss<<"fsproject: combines"<<fromstage<<" stage"<<fromstage<<" files undefined";
    throw(runtime_error(ss.str()));
  }
  int nerrs=0;
  for(unsigned int c1=0;c1<soutputrootvec.size();c1++){
    string tfile=soutputrootvec[c1];
    tfile.append(".EMprobs.out");
    thisstageemfiles.push_back(tfile);
    if(verbose) cout<<"Adding EM results for file "<<tfile<<endl;
    if( access( tfile.c_str(), F_OK ) == -1 ) { // file does not exist
      cerr<<"INFORMATION: EM output files from stage "<<fromstage<<" not found.";
      if(hpc==1) cerr<<" You are in hpc 1 mode, so stage "<<fromstage<<" has probably not yet been executed."<<endl;
      else cerr<<" You are not in HPC mode, which means something has gone wrong with -dos"<<fromstage<<". Investigate the logs."<<endl;
    stringstream ss;
    ss<<"fsproject: combines"<<stage<<" stage"<<stage<<" missing results";
      throw(runtime_error(ss.str()));
    }
    try{
      addEMobs(tfile,&Nevec,&muvec);
    }catch(runtime_error& e) {
	string swhat=e.what();
	if(swhat.compare("combineEM file error")==0){
	  ++nerrs;
	  cerr<<"WARNING: EMfile "<<tfile<<" could not be processed. This could mean that the parameters converged to an invalid point, or that the computation was not completed."<<endl;
	}
    }
  }

  // Check for success
  if((nerrs>=(int)soutputrootvec.size())| (Nevec.size()==0)| (muvec.size()==0)){
    cerr<<"ERROR: No EM runs contain valid parameter inference. Something went wrong."<<endl;
    cerr<<"If ChromoPainter stage"<<fromstage<<" did not complete, try \"-reset "<<fromstage<<" -go\" to rerun from that stage."<<endl;
    throw(runtime_error("EM estimation failed"));
  }else if(nerrs>0){
    cerr<<"WARNING: "<<nerrs<<" EM runs (of "<<soutputrootvec.size()<<") failed. Continuing by disregarding those estimates. This may lead to problems."<<endl;
  }

  // combine
  if(verbose) cout<<"Combining "<<Nevec.size()<<" Ne values from "<<thisstageemfiles.size()<<" files"<<endl;
  double tNeinf=0,tmuinf=0;
  tNeinf=0;
  tmuinf=0;
  for(unsigned int c1=0;c1<Nevec.size();c1++)tNeinf+=Nevec[c1];
  for(unsigned int c1=0;c1<muvec.size();c1++)tmuinf+=muvec[c1];
  tNeinf/=(double)Nevec.size();
  tmuinf/=(double)muvec.size();
  cout<<"Inferred Ne="<<tNeinf<<" and mu="<<tmuinf<<endl;
  if(fromstage==1){
    Neinf=tNeinf;
    muinf=tmuinf;
  }else{
    popNeinf=tNeinf;
    popmuinf=tmuinf;
  }
  stage=fromstage+1;
  validatedoutput[fromstage]=1;
}

bool FsProject::validateStage(int forstage,string ending,int expectedlines,string lastline,double minfrac,bool vverbose)
{
  int numexist=0;
  int numok=0;

  vector<int> isvalid = validOutput(forstage,ending,expectedlines,lastline,vverbose);
  for(unsigned int c1=0;c1<isvalid.size();c1++){
    if(isvalid[c1]>=0){ // file exists
      ++numexist;
    }
    if(isvalid[c1]==1){
      ++numok;
    }
  }
  
  if(numexist==0){
    if(vverbose) {
      cerr<<"INFORMATION:  Output from stage "<<forstage<<" not found.";
      if(hpc==1) cerr<<" You are in hpc 1 mode, so stage "<<forstage<<" has probably not yet been executed."<<endl;
      else cerr<<" You are not in HPC mode, which means something has gone wrong with -dos"<<forstage<<". Investigate the logs."<<endl;
    }
    return(false);
  }else if(numok<(int)ceil(isvalid.size()*minfrac)){
    if(vverbose ) {
      cerr<<"INFORMATION:  Output from stage "<<forstage<<" is incomplete. Found "<<numok<<" valid files, "<<numexist<<" existing files out of a total of "<<isvalid.size()<<" (tolerating a success rate of "<<minfrac<<")"<<endl;
      
      if(hpc==1) cerr<<" You are in hpc 1 mode, so stage "<<forstage<<" has probably partially completed. "<<endl;
      else cerr<<" You are not in HPC mode, which means something has gone wrong with -dos"<<forstage<<". Investigate the logs."<<endl;
    }

    stringstream ss;
    ss<<"fsproject: combines"<<forstage<<" stage"<<forstage<<" incomplete results";
    throw(runtime_error(ss.str()));
  }else if((int)(numok)!=(int)isvalid.size()){
    if(vverbose ) {
      cerr<<"INFORMATION:  Output from stage "<<forstage<<" is incomplete. Found "<<numok<<" valid files, "<<numexist<<" existing files out of a total of "<<isvalid.size()<<" (tolerating a success rate of "<<minfrac<<"). Although the output is incomplete, enough are completed for the run to continue."<<endl;
    }
  }
  return(true);    
}

vector<int> FsProject::validOutput(int forstage,string ending,int expectedlines,string lastline,bool vverbose)
{
  if(hpc) ensureCommandfile(forstage);
  vector<int> isvalid;
  vector<string> srootvec;
  switch(forstage){
  case 1:srootvec=s1outputrootvec;break;
  case 2:srootvec=s2outputrootvec;break;
    //  case 3:srootvec=s3outputrootvec;break;
    //  case 4:srootvec=s4outputrootvec;break;
  case 6:srootvec=s6outputrootvec;break;
  case 7:srootvec=s7outputrootvec;break;
    // case 9:
  default: throw(logic_error("validateStage: invalid stage"));
  }
  isvalid=vector<int>(srootvec.size(),0);
  
  for(unsigned int c1=0;c1<srootvec.size();c1++){
    string tfile=srootvec[c1];
    tfile.append(ending);
    if( access( tfile.c_str(), F_OK ) == -1 ) { // file does not exist
      if(verbose) cout<<"FILE "<<tfile<<" does not exist."<<endl;
      isvalid[c1]=-1;
    }else if(countLines(tfile) < expectedlines){ // if we've said we know the expected number of lines (+ve)
      // file exists but has the wrong number of rows
      if(verbose) cout<<"FILE "<<tfile<<" has "<<countLines(tfile)<<" lines, but expected "<<expectedlines<<endl;
    }else if(lastline.compare("")!=0){ // if we've said that we know the last line
      // check it
      if(lastline.compare(getLastLine(tfile))==0) isvalid[c1]=1;
    }else{
      isvalid[c1]=1;
    }
  }
  
  return(isvalid);
}


bool FsProject::validateStageCpEm(int forstage, bool vverbose)
{
  int emits=1;
  if(forstage==1) emits=s1emits;
  else if(forstage==6) emits=s6emits;
  return(validateStage(forstage,".EMprobs.out",(getIndsPerProc(forstage)*(emits+1)),"",0.75,vverbose));
}

bool FsProject::validateStageCpRun(int forstage,bool vverbose)
{   
  return(validateStage(forstage,".chunkcounts.out",getIndsPerProc(forstage)+1,"",1.0,vverbose));
}

bool FsProject::validateStage(int forstage,bool vverbose){
  if(forstage==1 || forstage==6){
    return(validateStageCpEm(forstage,vverbose));
  }else if(forstage==2 || forstage==7){
    return(validateStageCpRun(forstage,vverbose));
  }else{
    cout<<"WARNING: Validation of stage "<<forstage<<" not properly implemented. Returning valid."<<endl;
    return(true);
  }
}

void FsProject::removeCompletedCommands(int forstage)
{
  vector<int> isvalid;
  int emits=1;
  if(forstage==1) emits=s1emits;
  else if(forstage==6) emits=s6emits;

  if(forstage==1 || forstage==6) isvalid=validOutput(forstage,".EMprobs.out",
						     (getIndsPerProc(forstage)*(emits+1)),"",false);
  else if(forstage==2 || forstage==7) isvalid=validOutput(forstage,".chunkcounts.out",
						     (getIndsPerProc(forstage)+1),"",false);
  else throw (runtime_error("removeCompletedCommands: stage not implemented"));
  vector <string> sxoutputrootvec,sxcommands;
  switch(forstage){
  case 1: sxoutputrootvec=s1outputrootvec;sxcommands=s2commands;break;
  case 2: sxoutputrootvec=s2outputrootvec;sxcommands=s2commands;break;
    //  case 3: sxoutputrootvec=s3outputrootvec;break;
    //  case 4: sxoutputrootvec=s4outputrootvec;break;
  case 6: sxoutputrootvec=s6outputrootvec;sxcommands=s6commands;break;
  case 7: sxoutputrootvec=s7outputrootvec;sxcommands=s7commands;break;
  }
  if(isvalid.size()!=sxoutputrootvec.size()){
    cerr<<"removeCompletedCommands: have "<<isvalid.size()<<" commands to check but "<<sxoutputrootvec.size()<<" commands stored in the file!"<<endl;
    throw(logic_error("removeCompletedCommands"));
  }
  ////////////////
  int nkeep=0;
  for(int c1=isvalid.size()-1;c1>=0;--c1){
    if(isvalid[c1]==1){
      sxcommands.erase(sxcommands.begin()+c1);
      ++nkeep;
    }
  }
  switch(forstage){
  case 1: s1commands=sxcommands;break;
  case 2: s2commands=sxcommands;break;
    //  case 3: s3outputrootvec=sxoutputrootvec;break;
    //  case 4: s4outputrootvec=sxoutputrootvec;break;
  case 6: s6commands=sxcommands;break;
  case 7: s7commands=sxcommands;break;
  }
  if(verbose) cout<<"Reprocessing "<<sxcommands.size()<<" commands that were not correctly completed. Retained results for "<<nkeep<<" that were completed successfully."<<endl;
}

string FsProject::combineStageCpRunSingle(int chr,string sroot, string sxcombineargs, vector<string>srootvec,int forstage)
{// Combine stage 2/7 for either all chromosomes (if chr is negative) or just the specified one.

  string ts="combine";
  ts.append(sxcombineargs);
  if(fsmode>0){
    ts.append(" -v ");
  }else {
    ts.append(" -q ");
  }
  ts.append(" -o ");
  if(chr>=0){
    ostringstream ss;
    ss<<dirname<<"/"<<"stage"<<forstage<<"a/"<<"singlechromosome_file"<<chr+1;
    sroot=ss.str();
  }
  ts.append(sroot);
    
  unsigned int filestart=0,fileend=srootvec.size();

  if(chr>=0){
    int filesperchr=srootvec.size()/phasefiles.size();
    filestart=chr*filesperchr;
    fileend=(chr+1)*filesperchr;
  }
  
  for(unsigned int c1=filestart;c1<fileend;c1++){
    string tfile=srootvec[c1];
    tfile.append(".chunkcounts.out");
    ts.append(" ");
    ts.append(srootvec[c1]);
  }

  std::vector<char *> argv=converttoargv(ts);

  std::ostringstream ss;

  ss<<dirname<<"/stage"<<forstage<<"a/";
  ensureDirectory(ss.str());
  ss<<baseName(sroot)<<".stage"<<forstage;
  //  if(chr>=0) ss<<".file"<<chr;
  ss<<".combine.log";
  string logfile=ss.str();
  if(chr<0 && fsmode>0) cout<<"Combining stage"<<forstage<<" files to file root "<< sroot<<" (logging to "<<logfile<<")"<<endl;
  if(chr>=0  && fsmode>0) cout<<"Combining stage"<<forstage<<" (data file "<<chr+1<<") files to file root "<< sroot<<" (logging to "<<logfile<<")"<<endl;
  if(verbose) cout<<"RUNNING STAGE"<<forstage<<"a CMD:"<<ts<<endl;
  switchStdout(logfile.c_str());

  int rv=1;
  try{
    rv=chromocombine(argv.size(),argv.data());
  }catch(exception &x){
    cerr<<"Caught ChromoCombine error: "<<x.what()<<endl;
  }
  freeargv(argv);
  revertStdout();
  // Check that it ran correctly
  if((getLastLine(logfile).compare(ccsuccesstext)!=0) | (rv!=0)){
    cerr<<"ChromoCombine failed! See log file ("<<logfile<<") for details."<<endl;
    cerr<<"If ChromoPainter stage"<<forstage<<" did not complete, try \"-reset "<<forstage<<" -go\" to rerun from that stage."<<endl;
      throw(runtime_error("chromocombine"));
  }
    // checkthe chunkcount file, extract C and check it worked OK
  return(logfile);
}

void FsProject::combineStageCpRun(int forstage)
{
  vector<string> srootvec;
  string sroot;
  string sxcombineargs;
  if(forstage==2){
    srootvec=s2outputrootvec;
    sxcombineargs=s2combineargs;
    sroot=cproot;
  }else{
    srootvec=s7outputrootvec;
    sxcombineargs=s7combineargs;
    sroot=popcproot;
  }

  // Combine the chromosomes one-by-one
  for(int chr=0;chr<(int)phasefiles.size();++chr){
    combineStageCpRunSingle(chr,sroot,sxcombineargs,srootvec,forstage);
  }
  // Combine all the data
  string logfile = combineStageCpRunSingle(-1,sroot,sxcombineargs,srootvec,forstage);
  
  string chunkfile=sroot;
  chunkfile.append(".chunkcounts.out");

  ifstream file(chunkfile.c_str());
  string header;
  getline(file, header);
  file.close();
  std::vector<std::string> headervec=split(header,' ');
  double tcval=-1;
  if(headervec.size()==2){
    if(headervec[0].compare("#Cfactor")==0){
      istringstream ( headervec[1] ) >> tcval;
    }
  }
  if(tcval<0){
    cerr<<"ChromoCombine failed! The combined file "<<chunkfile<<" doesn't contain valid information about 'c'. Something went wrong; See log file ("<<logfile<<") for details."<<endl;
    throw(runtime_error("chromocombine"));
  }else if((fsmode>0) && (tcval==0)) {
    if(forstage==2) cerr<<constccerror<<endl;
    else cerr<<constccerror7<<endl;
    throw(runtime_error("chromocombine"));
  }else if(fsmode>0){
    cout<<"Successfully run ChromoCombine stage! Inferred a 'c' value of "<<tcval<<endl;
    if(linkagemode.compare("linked")==0){
      if(((tcval>2) | (tcval<0.1))) cerr<<"WARNING: in linked mode, we usually expect 'c' to be between 0.1 and 2. You are advised to examine whether there have been processing problems (for example, has the parameter inference become stuck at parameters implying effectively unlinked data?)"<<endl;
    }else{
      if((tcval>0.1)) cerr<<"WARNING: in unlinked mode, we usually expect 'c' to be less than 0.1. You are advised to examine whether there have been processing problems."<<endl;
    }
  }
  if(forstage==2){
    cpchunkcounts=chunkfile;
    cval=tcval;
  }else{
    popcpgenomelen=sroot;
    popcpgenomelen.append(".chunklengths.out");
  }
  stage=forstage+1;
  validatedoutput[forstage]=1;
}

void FsProject::combineStage3()
{ // check MCMC output
  bool haveoutput=true;
  if(fsmcmcoutputvec.size()==0){haveoutput=false;
  }else if( access( fsmcmcoutputvec[0].c_str(), F_OK ) == -1 )haveoutput=false;
  if(!haveoutput){ // file does not exist
    //    cerr<<"Combining stage3 files requires stage3 to have been run. Has stage3 been run? Have all \"*.xml*\" files been recovered from remote processing?"<<endl;
    cerr<<"If stage3 did not complete, try \"-reset 3 -go\" to rerun from stage3."<<endl;

    throw(runtime_error("fsproject: combines3 stage3 missing results"));
  }

  bool converged=0;
  try{
    converged=mcmcConvergence();
}catch(exception &x){
    cerr<<"Error opening file!"<<endl<<x.what()<<endl; 
    throw(runtime_error("-combines3: MCMC reading error. Need to rerun MCMC!"));
  }
  double maxgr=mcmcGR[0];
  for(unsigned int c1=1;c1<mcmcGR.size();c1++) maxgr=max(maxgr,mcmcGR[c1]);
  if(!converged){
    // Reject the MCMC run, double the duration and continue the run
    validatedoutput[3]=0;
    stage=3;
      cout<<"WARNING: Failed Gelman-Rubin MCMC diagnostic check with maximum potential scale reduction factor "<<maxgr<<" (threshold "<<threshGR<<")"<<endl;
  }else{
    if(verbose){
      if(threshGR>0) cout<<"Passed";
      else cout<<"Skipped";
      cout << " Gelman-Rubin MCMC diagnostic check with maximum potential scale reduction factor "<<maxgr<<" (threshold "<<threshGR<<")"<<endl;
    }
    // Move it to the root of the directory
    size_t found =fsmcmcoutputvec[0].find("/stage3/");
    stringstream ss;
    ss<<dirName(dirname)<<fsmcmcoutputvec[0].substr(found+8);
    if(verbose) cout<<"Moving "<<fsmcmcoutputvec[0]<<" to "<<ss.str()<<endl;
    rename(fsmcmcoutputvec[0].c_str(),ss.str().c_str());
    fsmcmcoutputvec[0]=ss.str();
    stage=4;
    validatedoutput[3]=1;
  }
}

void FsProject::combineStage4()
{ // check TREE output
  if(haveOutput(4)){ // THIS NEEDS TO BE MORE THOROUGH!
    // Move it to the root of the directory
    size_t found =fstreeoutputvec[0].find("/stage4/");
    stringstream ss;
    ss<<dirName(dirname)<<fstreeoutputvec[0].substr(found+8);
    if(verbose) cout<<"Moving "<<fstreeoutputvec[0]<<" to "<<ss.str()<<endl;
    rename(fstreeoutputvec[0].c_str(),ss.str().c_str());
    fstreeoutputvec[0]=ss.str();

    stage=5;
    validatedoutput[4]=1;
  }else{
    cerr<<"If stage4 did not complete, try \"-reset 4 -go\" to rerun from stage4."<<endl;
    throw(runtime_error("fsproject: combines4 stage4 missing results"));
  }
}

void FsProject::writeHpcStage3(string cmdfile)
{
  ensurefsmcmc();
  ostringstream ss;
  ss<<filename<<" -hpc 0 -allowdep 0 -makes3 -dos3 -combines3 -allowdep 1 -hpc 1 -v";
  cout<<"CREATING S3 CMD:"<<ss.str()<<endl;
  s3commands.clear();
  s3logfiles.clear();
  s3commands.push_back(ss.str());
  s3logfiles.push_back(fsmcmcoutput);
  s3logfiles[0].append(".hpcs3.log");
  writeStringVectorToFile(s3commands,cmdfile,s3logfiles,true);
}

void FsProject::createIdFile(string idfile)
{
  if(phasefiles.size()==0) {
    cerr<<"ERROR: Require phase file to create idfile!"<<endl;
    throw(runtime_error("-createid"));
  }
  ///////////// write the ID file
  int nhaps=getNhapsFromFile(0);
  int nindst=nhaps/ploidy;
  if(nindst!=(int)ceil(nhaps/ploidy)) {
    cerr<<"ERROR: phase file contains a number of haps that is not a multiple of the ploidy!"<<endl;
    throw(runtime_error("-createid"));
  }

  filebuf fb;
  try{
    fb.open (idfile.c_str(),ios::out);
  }catch(exception &x){
    cerr<<"Error opening file!"<<endl<<x.what()<<endl; 
    throw(runtime_error("-createid: cannot write to file!"));
  }
  ostream os (&fb);
  for(int i=0;i<nindst;i++){
    os<<"IND"<<i+1<<endl;
  }
  fb.close();
}

int FsProject::numRecipients(int forstage)
{
  if(forstage<0 && popidfile.compare("")==0){
    forstage=1;
  }else if(forstage<0) {
    forstage=6;
  }
  if(forstage<5){
    return (nindsUsed);
  }else return(nrecipsDonor);
}

int FsProject::recommendN()
{
  int ncmds=getCommandFileCount();
  if(ncmds==0) ncmds=getCommandFileLength();
  if(stage==3 && ncmds==1){
#ifdef _OPENMP
    return(nummcmcruns);
#endif
  }
  if(ncmds<8) return(ncmds);
  else return(8);
}

int FsProject::recommendM()
{
  // data:
  int ncmds=getCommandFileCount();
  if(ncmds==0) ncmds=getCommandFileLength();
  int n=recommendN();

  if(stage==3 && ncmds==1){
#ifdef _OPENMP
    return(nummcmcruns);
#endif
  }
  // Linear growth (c*x) until x0 and then logarithmic (alpha*log(a*x+b), matching the value and gradient at x0.
  double x0=160;// max number of commands for which we parallelize completely (approached exponentially slowly)
  double c=1.0/(double)n;// 1/number of commands per batch
  double alpha=8.0; // exponent in the growth rate of the logarithmic section of the curve

  if(ncmds<=x0) return(n);
  double a=exp(c*x0/alpha)*c/alpha;
  double b=exp(c*x0/alpha)-a*x0;
  int mrec=ceil(alpha*log(a*ncmds+b));
  
  return((int)floor(ncmds/mrec));
}

bool FsProject::haveOutput(int stage,bool forcombined)
{
  // If we previously checked the output and found it was ok:
  if(missingfiles[stage]==0 && (!forcombined)) return(true); 

  if(stage==1){
    if(s1outputrootvec.size()==0) return(false);
    return(validateStageCpEm(stage,true));
    missingfiles[stage]=0;
  }else if((stage==2) && !(forcombined)){
    if(s2outputrootvec.size()==0) return(false);
    return(validateStageCpRun(stage,true));
  }else if((stage==2) && (forcombined)){
    if(cpchunkcounts.compare("")==0) return(false);
    string f=cpchunkcounts;
    if( access( f.c_str(), F_OK ) == -1 ) { // file doesn't exist
      //      if(verbose) cout<<"Do not have stage2 combined output "<<f<<endl;
      return(false);
    }
  }else if((stage==3) && (!forcombined)){
    if(fsmcmcoutputvec.size()==0) return(false);
    for(unsigned int i=0;i<fsmcmcoutputvec.size();i++){
      string f=fsmcmcoutputvec[i];
      if( access( f.c_str(), F_OK ) == -1 ) { // file doesn't exist
	//	if(verbose) cout<<"Do not have stage3 output "<<f<<endl;
	return(false);
      }
      if (getLastLine(f).compare(string("</outputFile>"))!=0) {
	if(verbose) cout<<"Do not have complete stage3 output "<<f<<endl;
	return(false);
      }
    }
  }else if((stage==3) && (forcombined)){
    if(fsmcmcoutputvec.size()==0) return(false);
    for(unsigned int i=0;i<fsmcmcoutputvec.size();i++){
      string f=fsmcmcoutputvec[i];
      if( access( f.c_str(), F_OK ) == -1 ) { // file doesn't exist
	//	if(verbose) cout<<"Do not have stage3 output "<<f<<endl;
	return(false);
      }
      if (getLastLine(f).compare(string("</outputFile>"))!=0) {
	if(verbose) cout<<"Do not have complete stage3 output "<<f<<endl;
	return(false);
      }
    }
  }else if(stage==4){
    if(fstreeoutputvec.size()==0) return(false);
    for(unsigned int i=0;i<fstreeoutputvec.size();i++){
      string f=fstreeoutputvec[i];
      if( access( f.c_str(), F_OK ) == -1 ) { // file doesn't exist
	//	if(verbose) cout<<"Do not have stage4 output "<<f<<endl;
	return(false);
      }
      if (getLastLine(f).compare(string("</outputFile>"))!=0) {
	if(verbose) cout<<"Do not have complete stage4 output "<<f<<endl;
	return(false);
      }
    }
  }else if(stage==5){
    if(popidfile.compare("")==0) {
	if(verbose) cout<<"Do not have complete stage5 output in the form of a popidfile."<<endl;
      return(false);
    }
  }else if(stage==6){
    if(s6outputrootvec.size()==0) return(false);
    return(validateStageCpEm(stage,true));
  }else if((stage==7) && !(forcombined)){
    if(s7outputrootvec.size()==0) return(false);
    return(validateStageCpRun(stage,true));
  }else if((stage==7) && (forcombined)){
    string f=popcpgenomelen;
    if( access( f.c_str(), F_OK ) == -1 ) { // file doesn't exist
      //      if(verbose) cout<<"Do not have stage2 combined output "<<f<<endl;
      return(false);
    }
  }
  return(true);
}

void FsProject::makePopIds()
{
  FsDonor donormodel(popidfile,donoridfile,dirname,0,(int)phasefiles.size(),
		     fsmode,verbose);
  donormodel.createRequiredFiles();
  donoridfile=donormodel.getDonorFile();
  gtdonors=donormodel.donorLabels();
  gtrecips=donormodel.recipLabels();
  popcmdargs = donormodel.getCommandContent();
  nrecipsDonor = donormodel.numRecipients();
  refidfiles = donormodel.getReferenceIdFiles();
  refdonorfile = donormodel.getReferenceDonorFile();
    
  if(verbose){
    for(unsigned int c1=0;c1<popcmdargs.size();++c1){
      cout<<popcmdargs[c1]<<endl;
    }
  }
  validatedoutput[5]=1;
  stage=6;
}

bool FsProject::canDo(string cmd, bool vverbose)
{
  //////
  // TO DO : check for data for -makes1 and -makes2 
  if(cmd.compare("-makes1")==0){
    if(ninds<=0) {
      if(vverbose) cout<<"ninds is not defined."<<endl;
      return(false);
    }
  }else if(cmd.compare("-dos1")==0){
    if(s1commands.size()==0) {
      if(vverbose) cout<<"s1commands has not been defined."<<endl;
      return(false);
    }
  }else if(cmd.compare("-combines1")==0){
    if(!haveOutput(1)) {
      if(vverbose) cout<<"s1 output not detected."<<endl;
      return(false);
    }
  }else if(cmd.compare("-makes2")==0){
    if(((muinf<=0)||(Neinf<=0))&(linkagemode.compare("linked")==0)){
      if(vverbose) cout<<"s2 muinf and/or Neinf are not defined (linked mode)."<<endl;
      return(false);
    }
    if((ninds<=0)&(linkagemode.compare("unlinked")==0)) {
      if(vverbose) cout<<"ninds is not defined (unlinked mode)."<<endl;
      return(false);
    }
  }else if(cmd.compare("-dos2")==0){
    if(s2commands.size()==0) {
      if(vverbose) cout<<"s2commands has not been defined."<<endl;
      return(false);
    }
  }else if(cmd.compare("-combines2")==0){
    if(!haveOutput(2)) {
      if(vverbose) cout<<"s2 output not detected."<<endl;
      return(false);
    }
  }else if(cmd.compare("-makes3")==0){
    if(!haveOutput(2,true)) {
      if(vverbose) cout<<"s2 combine output not detected."<<endl;
      return(false); // check whether combine done
    }
  }else if(cmd.compare("-dos3")==0){
    if(s3commands.size()==0) {
      if(vverbose) cout<<"s3commands has not been defined."<<endl;
      return(false);
    }
  }else if(cmd.compare("-combines3")==0){
    if(!haveOutput(3)) {
      if(vverbose) cout<<"s3 output not detected."<<endl;
      return(false);
    }
  }else if(cmd.compare("-hpcs3")==0){ // -hpcs3 replaces -makes3 and -dos3
    if(!haveOutput(2,true)) {
      if(vverbose) cout<<"s2 combine output not detected."<<endl;
      return(false); // check whether combine done
    }
  }else if(cmd.compare("-makes4")==0){
    if(!validatedoutput[3]) {
      if(vverbose) cout<<"s3 output not validated."<<endl;
      return(false);
    }
  }else if(cmd.compare("-dos4")==0){
    if(s4commands.size()==0) {
      if(vverbose) cout<<"s4commands not defined."<<endl;
      return(false);
    }
  }else if(cmd.compare("-combines4")==0){
    if(!haveOutput(4)) {
      if(vverbose) cout<<"s4 output not detected."<<endl;
      return(false);
    }
    // fs 2.1 properties from here
  }else if(cmd.compare("-dos5")==0){
    if(!validatedoutput[4]) {
      if(vverbose) cout<<"s4 output not validated."<<endl;
      return(false);
    }
  }else if(cmd.compare("-combines5")==0){
    if(!haveOutput(5)) {
      if(vverbose) cout<<"s5 output not detected."<<endl;
      return(false);
    }
  }else if(cmd.compare("-makes6")==0){
    if(!validatedoutput[5]) {
      if(vverbose) cout<<"s5 output not validated."<<endl;
      return(false);
    }
  }else if(cmd.compare("-dos6")==0){
    if(s6commands.size()==0) {
      if(vverbose) cout<<"s6commands not defined."<<endl;
      return(false);
    }
  }else if(cmd.compare("-combines6")==0){
    if(!haveOutput(6)) {
      if(vverbose) cout<<"s6 output not detected."<<endl;
      return(false);
    }
  }else if(cmd.compare("-makes7")==0){
    //
    if(((popmuinf<=0)||(popNeinf<=0))&(linkagemode.compare("linked")==0)) {
      if(vverbose) cout<<"s7 popmuinf and/or popNeinf not defined (linked mode)."<<endl;
      return(false);
    }
    // TODO *** It would be wise to check that the popidfile counts up correctly as containing populations
    if((ninds<=0)&(!haveOutput(5))&(linkagemode.compare("unlinked")==0)) {
      if(vverbose) cout<<"s5 output not detected (unlinked mode)."<<endl;
      return(false);
    }
  }else if(cmd.compare("-dos7")==0){
    if(s7commands.size()==0) {
      if(vverbose) cout<<"s7commands not defined."<<endl;
      return(false);
    }
  }else if(cmd.compare("-combines7")==0){
    if(!haveOutput(7)) {
      if(vverbose) cout<<"s7 output not detected."<<endl;
      return(false);
    }
  }else if(cmd.compare("-makes8")==0){
    if(!haveOutput(7,true)) {
      if(vverbose) cout<<"s7 combine output not detected."<<endl;
      return(false); // check whether combine done
    }
  }else if(cmd.compare("-go")==0){
    if((popidfile.compare("")==0) & !validatedoutput[4]) {
      if(vverbose) cout<<"s4 output not validated (in standard fs mode)."<<endl;
      return(false); // standard fs mode
    }
    if((popidfile.compare("")!=0) & !validatedoutput[7]) {
      if(vverbose) cout<<"s7 output not validated (in admixture mode)."<<endl;
      return(false); // admixture mode
    }
  }else if(cmd.compare("-makes9")==0){
    if((popidfile.compare("")!=0) & !validatedoutput[7]) {
      if(vverbose) cout<<"s7 output not validated (in admixture mode)."<<endl;
      return(false); // for making gt
    }
  }else if(cmd.compare("-writes9")==0){
    if(gtroot.compare("")==0) {
      if(vverbose) cout<<"-makes9 not run."<<endl;
      return(false); // writing gt
    }
  }else if(cmd.compare("-gt")==0){
    if(gtroot.compare("")==0) {
      if(vverbose) cout<<"-gt not run."<<endl;
      return(false); // writing gt
    }
  }
  return(true);
}

std::vector<std::string> FsProject::getDependencies(std::vector<std::string> args)
{
  // The dependency tree/pipeline. It is actually trivial in fs2.0:
  // makes1->dos1->combines1->makes2->dos2->combines2->makes3->dos3->makes4->dos4
  // with complications because do becomes write for hpc mode
  // -hpcs3 complicates things further in hpc mode!
   // in fs2.1 we add:
  // dos5->combines5->makes6->dos6->combnines6
  // dos5 does not get hpc'd
 
  std::vector<std::string> ret;
  if(args[0].compare("-gt")==0){
    if (!canDo(args[0]))  ret.push_back("-writes9");
  }
  if(args[0].compare("-writes9")==0){
    if (!canDo(args[0]))  ret.push_back("-makes9");
  }
  if(args[0].compare("-makes9")==0){
    if (!canDo(args[0]))  ret.push_back("-combines7");
  }
  if(args[0].compare("-go")==0){
     if((popidfile.compare("")==0) & !canDo(args[0])) ret.push_back("-combines4");
     else if ((popidfile.compare("")!=0) & !canDo(args[0]))  ret.push_back("-combines7");
  }
  if(args[0].compare("-combines7")==0){
    if(!canDo(args[0])){
      if(hpc) ret.push_back("-writes7");
      if(hpc!=1) ret.push_back("-dos7");
    }
  }
  if((args[0].compare("-dos7")==0)|(args[0].compare("-writes7")==0)){
    if(!canDo("-dos7")) ret.push_back("-makes7");
  }
  if(args[0].compare("-makes7")==0){
    if(whichLinkagemode().compare("linked")==0){
      if(!canDo(args[0])) ret.push_back("-combines6");
    }
    if(whichLinkagemode().compare("unlinked")==0){
      if(!canDo(args[0])) ret.push_back("-combines5");
    }
  }  
  if(args[0].compare("-combines6")==0){
    if(!canDo(args[0])){
      if(hpc) ret.push_back("-writes6");
      if(hpc!=1) ret.push_back("-dos6");
    }
  }
  if((args[0].compare("-dos6")==0)|(args[0].compare("-writes6")==0)){
    if(!canDo("-dos6")) ret.push_back("-makes6");
  }
  if(args[0].compare("-makes6")==0){
    if(!canDo(args[0])) ret.push_back("-combines5");
  }
  if(args[0].compare("-combines5")==0){
    // Either we have been given a popidfile, in which case we should call countdata,
    if(haveOutput(5) && nsnpsvec.size()==0) ret.push_back("-countdatapop");
    // or we should try to construct it from an mcmc run:
    else if(!canDo(args[0])) ret.push_back("-dos5");
  }
  if(args[0].compare("-dos5")==0){
    if(!canDo(args[0])) ret.push_back("-combines4");
  }

  // Everything from fs2.0
  if(args[0].compare("-combines4")==0){
    if(!canDo(args[0])){
      if(hpc) ret.push_back("-writes4");
      if(hpc!=1) ret.push_back("-dos4");
    }
  }
  if((args[0].compare("-dos4")==0)|(args[0].compare("-writes4")==0)){
    if(!canDo("-dos4")) ret.push_back("-makes4");
  }
  if(args[0].compare("-makes4")==0){
    if(!canDo(args[0])) ret.push_back("-combines3");
  }
  if(args[0].compare("-combines3")==0){
    if(!canDo(args[0])) {
      int usehpcs3=0;
#ifdef _OPENMP
      usehpcs3=1;
#endif
      if(hpc==1 && !usehpcs3) {ret.push_back("-writes3");
      }else if(hpc==1 && usehpcs3) {ret.push_back("-hpcs3");
      }else ret.push_back("-dos3"); //!hpc mode
    }
  }
  if((args[0].compare("-dos3")==0)|(args[0].compare("-writes3")==0)){
    if(!canDo("-dos3")) ret.push_back("-makes3");
  }
  if(args[0].compare("-makes3")==0){
    if(!canDo(args[0])) ret.push_back("-combines2");
  }
  if(args[0].compare("-hpcs3")==0){
    if(!canDo(args[0])) ret.push_back("-combines2");
  }
  if(args[0].compare("-combines2")==0){
    if(!canDo(args[0])) {
      if(hpc) ret.push_back("-writes2");
      if(hpc!=1)  ret.push_back("-dos2");
    }
  }
  if((args[0].compare("-dos2")==0)|(args[0].compare("-writes2")==0)){
    if(!canDo("-dos2")) ret.push_back("-makes2");
  }
  if(args[0].compare("-makes2")==0){
    ensureCpRunARoot(2);
    if(whichLinkagemode().compare("linked")==0){
      if(!canDo(args[0])) ret.push_back("-combines1");
    }
    if(whichLinkagemode().compare("unlinked")==0){
      if(!canDo(args[0])) ret.push_back("-countdata");
    }
  }
  if(args[0].compare("-combines1")==0){
    if(!canDo(args[0])) {
      if(hpc) ret.push_back("-writes1");
      if(hpc!=1)  ret.push_back("-dos1");
    }
  }
  if((args[0].compare("-dos1")==0)|(args[0].compare("-writes1")==0)){
    if(!canDo("-dos1")) ret.push_back("-makes1");
  }
  if(args[0].compare("-makes1")==0){
    if(!canDo(args[0])) ret.push_back("-countdata");
    ensureCpEmRoot(1);
  }

  return(ret);
}

int FsProject::parIndex(std::string partest) {
  for(unsigned int c1=0;c1<pars.size();c1++){
    if(partest.compare(pars[c1].getName())==0){
      return(c1);
    }
  }
  return(-1);
}

int FsProject::cmdIndex(std::string cmdtest) {
  for(unsigned int c1=0;c1<cmds.size();c1++) {
    if(cmdtest.compare(cmds[c1].getName())==0){
      return(c1);
    }
  }
  return(-1);
}

std::string FsProject::cmdInfo(std::string cmd,bool statetype){
  std::ostringstream ss;
  while(cmd.substr(0,1).compare("-")==0){ // strip leading -
    cmd=cmd.substr(1,cmd.size());
  }
    
  // Parameters
  int parindex=parIndex(cmd);
  if(parindex>=0) {
    if(statetype) ss<<"Parameter ";
    ss<<pars[parindex].getName()<<" : "<<pars[parindex].getHelp();
    return(ss.str());
  }

  // Actions
  int cmdindex=cmdIndex(cmd);
  if(cmdindex>=0) {
    if(statetype) ss<<"Action    ";
    ss<<"-"<<cmds[cmdindex].getName();
    if(cmds[cmdindex].getShortargs().compare("")!=0)ss<<" "<<cmds[cmdindex].getShortargs();
    ss<<" : "<<cmds[cmdindex].getHelp();
    return(ss.str());
  }

  ss<<cmd<<" : "<<string("No help available. Is this a valid action or parameter?");
  return(ss.str());
}

bool FsProject::checkArgs(std::vector<std::string> args,int minargs,int maxargs)
{
  int nargs= (int)args.size()-1;
  if(maxargs==-1) maxargs=minargs;
  if(nargs<minargs || (nargs>maxargs && maxargs>=0)) {
    cerr<<"ERROR: Incorrect number of arguments: "<<args[0]<<" takes "<<minargs;
    if(maxargs== -2) {cerr<<"-INF";
    }else if(maxargs>minargs) cerr<<"-"<<maxargs;
    cerr<<" argument(s)."<<endl;
    cerr<<"Help for "<<cmdInfo(args[0])<<endl;
    throw(runtime_error(args[0]));
  }
  return(true);
}

int FsProject::checkStage(std::string args,std::string val) {
  int tpar=parIndex(args);
  if(tpar<0) {
    return(0);
  }
  try{checkStage(args,pars[tpar].getStage()); 
  }catch(exception &x) {
    return(0);
  }
  if(verbose) cout<<"Successfully read "<<args<<":\""<<val<<"\""<<endl;
  return(1);
}

void FsProject::checkStage(std::string args,int maxstage){
  //  cout<<"checkstage args "<<args<<" stage "<<stage<<" maxstage "<<maxstage<<endl;
    if(stage>maxstage && maxstage>=0) {
      cerr<<"ERROR: Tried to set "<<args<<" when in processing stage "<<stage<<", but this can only be done before stage "<<maxstage<<".  "<<constdupreset<<endl ;
      throw(runtime_error("commands out of order"));
    }
  
}

void FsProject::checkStage(std::vector<std::string> args,int maxstage){
  checkStage(args[0],maxstage);
}

void FsProject::docmd(std::string args)
{
  vector<string> fullargs;
  fullargs.push_back(args);
  docmd(fullargs);
}

void FsProject::docmd(std::vector<std::string> args)
{
  bool found=false;
  if(args.size()==0) throw(logic_error("fsproject: docmd empty command!"));
  //  if(verbose) cout<<"Processing docmd("<<args[0]<<",...)"<<endl;
  ///////////////////////////////////////
  // resolve dependencies
  std::vector<std::string> prev=getDependencies(args);
  if(prev.size()>0) {
    if(allowdep){
      if(verbose) cout<<"NOTE: Determined that "<<prev[0]<<" needs to be run first; will attempt to perform this and return to this command."<<endl;
      for(unsigned int c1=0;c1<prev.size();c1++){ // do all the dependencies, not just the first one
	std::vector<std::string> thisprev;
	thisprev.push_back(prev[c1]);
	docmd(thisprev);
      }
    }else{
      cerr<<"ERROR: Command "<<args[0]<<" requires that "<<prev[0]<<" has been run, but dependency resolution has been disabled. You need to manually resolve the dependencies."<<endl;
      throw(runtime_error("dependencies unmet"));
    }
  }
 
  ///////////////////////////////////////
  // HELPER FUNCTIONS
  if(args[0].compare("-import")==0){ // import settings from a text file
    found=checkArgs(args,1);
    string oldfile=filename;
    filename=args[1];
    try{
      readFromFile();
    }catch(exception &e){
      cerr<<"Error reading settings file!"<<endl<<e.what()<<endl; 
      filename=oldfile;
      throw(runtime_error("Command line failure"));
    }
    filename=oldfile;
  }else if(args[0].compare("-duplicate")==0){ // duplicate a project to a new settings file
    found=checkArgs(args,2);
    int newstage=0;
    try{newstage=stringToInt(args[1]);
    }catch(exception &x){
      cerr<<"-duplicate argument error: first argument should be an integer."<<endl;
      cerr<<"Help for "<<cmdInfo(args[0])<<endl;
      throw(runtime_error("-duplcate arguments error"));
    }
    resetToStage(newstage);
    createDuplicated(newstage,args[2]);
  }else if(args[0].compare("-reset")==0){ // duplicate a project to a new settings file
    found=checkArgs(args,1);
    int newstage=0;
    istringstream ( args[1] ) >> newstage;
    resetToStage(newstage);
  }else if(args[0].compare("-paintercmd")==0){// change the chromopainter command to an external call
    found=checkArgs(args,1);
    cerr<<"IMPORTANT: Specifying EXTERNAL command "<<args[1]<<" to be used for the ChromoPainter steps. This enables -hpc 1 mode, and therefore requires manually running the stages. It disables subsampling snps, which slows things down; use -s6indfrac 0.1 or similar. also sets sets -s1args and -s6args to \"-in -iM\" (removing --emfilesonly, which is fs specific). This increases the amount of data retained from the EM step - you can remove all files but the EMprobs.out after EM parameter estimation is performed."<<endl;
    hpc=1;
    s1args=string("-in -iM");
    s6args=string("-in -iM");
    s1minsnps=-1;
    s1snpfrac=-1;
    s6minsnps=-1;
    s6snpfrac=-1;
    paintercmd=args[1];
  }else if(args[0].compare("-cando")==0){
    found=checkArgs(args,1);
    string cmd=args[1];
    if(fsmode>0) cout<<"Checking command: "<<cmd<<"... will exit without modifying command files."<<endl;    
    bool found=false;
    bool res=canDo(cmd,true);
    int tpar=-1;
    int maxstage=-1;
    bool cando=true;
    
    if(cmdIndex(cmd)>=0) found=true;
    if(found) {
      if(fsmode>0) cout<<cmd<<" is a valid command."<<endl;
      tpar=cmdIndex(cmd);
      maxstage=cmds[tpar].getStage();
    }else {
      if(parIndex(cmd)>=0) found=true;
      if(found) {
	tpar=parIndex(cmd);
	maxstage=pars[tpar].getStage();
	if(fsmode>0) cout<<cmd<<" is a valid parameter."<<endl;
      }
    }
    
    if(!found) {
      if(fsmode>0) cout<<cmd<<" does not appear to be a command or parameter"<<endl;
      cando=false;
    }
    if(cando) {
      if(fsmode>0) cout<<cmdInfo(cmd,true)<<endl;
      FsSettingsValue val(args);
      if(stage>maxstage && maxstage>=0){
	if(fsmode>0) cout<<"Using "<<cmd<<" must be done in stage "<<maxstage<<" or less, but we are in stage "<<stage<<endl;
	cando=false;
      }
    }
    if(cando){
      args.erase(args.begin());
      args[0]=string("-").append(cmd);
      
      vector<string> dep=getDependencies(args);
      for(unsigned int i=0;i<dep.size();i++){
	cando=false;
	if(fsmode>0) cout<<"Command "<<cmd<<" has unmet dependency: "<<dep[i]<<" in the current mode. Only the first dependency is given here."<<endl;
      }
    }
    if(fsmode>0) {
      if(cando) cout<<"-cando "<<cmd<<": true"<<endl;
      else cout<<"-cando "<<cmd<<": false"<<endl;
    }
    return;
  }

  ///////////////////////////////////////
  // DATA IO

  if(args[0].compare("-createid")==0){ // ID FILE (Important)
    found=checkArgs(args,1);
    checkStage(args,1);
    if(verbose) cout<<"CREATING ID FILE: "<<args[1]<<endl;
    cerr<<"WARNING: When creating the ID file, we must have the correct -phasefiles. The ploidy must also be set (diploid by default, -haploid for haploids)."<<endl;
    createIdFile(args[1]);
    idfile=args[1];
  }

  /*if(args[0].compare("-recombfiles")==0){ // RECOMBINATION FILES (Important)
    found=checkArgs(args,1,-2);
    checkStage(args,1);
    vector<string> oldrec=recombfiles;
    for(unsigned int i=1;i<args.size();i++){
      recombfiles.push_back(args[i]);
    }
    if(unique(recombfiles).size()!=recombfiles.size()){// duplicates
      recombfiles=oldrec;
      cerr<<"ERROR: duplicated recomb files provided! You only need to specify each recombination file once, and you should not use --recombfiles in future commands."<<endl;
      throw(runtime_error("duplicate stage files"));
   }

    s12inputtype=string("phase");
    linkagemode=string("linked");
    if(verbose) cout<<"RECOMB FILES: read "<<recombfiles.size()<<" starting with "<<recombfiles[0]<<endl;
  }*/

/*  if(args[0].compare("-donorfile")==0){
    found=checkArgs(args,1,-2);
    checkStage(args,1);
    cerr<<"ERROR: Donor files are currently not supported in fs. Use -idfile only."<<endl;
    throw(runtime_error(args[0]));
        if(args.size()!=2) {
      cerr<<"ERROR: Incorrect number of arguments: -donorfile takes exactly 1 argument"<<endl;
      throw(runtime_error(args[0]));
    }
    s12donorfile=args[1];
    if(verbose) cout<<"DONOR FILE: "<<s12donorfile<<endl;
    found=true;
  }*/

  if(args[0].compare("-popidfile")==0){
    found=checkArgs(args,1);
    //    checkStage(args,1);
    popidfile=args[1];
    if(idfile.compare("")==0) idfile=popidfile;
    if(verbose) cout<<"REFERENCE POPULATION ID FILE: "<<popidfile<<endl;
    found=true;
  }
  
  if(args[0].compare("-popdonorfile")==0){
    found=checkArgs(args,1);
    donoridfile=args[1];
    if(verbose) cout<<"REFERENCE POPULATION DONOR FILE: "<<donoridfile<<endl;
    found=true;
  }
 

  if(args[0].compare("-datalist")==0){ // DATA LIST (Not implemented)
    found=checkArgs(args,1);
    checkStage(args,1);
    cout<<"NOT IMPLEMENTED: But this should set data from a file"<<endl;
  } 

  //////////////////////////
  if(args[0].compare("-countdata")==0){ /// Count and check the data
    found=checkArgs(args,0);
    checkStage(args,1);
    countData();
  }

  if(args[0].compare("-countdatapop")==0){ /// Count and check the data
    found=checkArgs(args,0);
    checkStage(args,5);
    countDataPop();
  }
  ///////////////////////////////////
  // Validate commands
  if(args[0].substr(0,10).compare("-validates")==0){ /// Check MAKE commands
    if(args[0].size()<=10 && args.size()==1){
      cerr<<"ERROR: -validates provided without any stage specified"<<endl;
      throw(runtime_error(args[0]));
    }
    int forstage=atoi(args[0].substr(10,1).c_str());
    if(args.size()>1) forstage=atoi(args[1].c_str());
    found=checkArgs(args,0,1);
    checkStage(args,forstage);
    if(fsmode>0) cout<<"Validating stage: "<<forstage<<"... will exit without modifying command files."<<endl;
    try{
      validateStage(forstage,true);
    }catch(exception &x){
      cout<<"VALIDATION FAIL"<<endl;
      return;
    }
    cout<<"VALIDATION SUCCESSFUL"<<endl;
    return;
  }

  ///////////////////////////////////
  // Making commands
  if(args[0].substr(0,6).compare("-makes")==0){ /// Check MAKE commands
    if(args[0].size()<=6){
      cerr<<"ERROR: -makes provided without any stage specified"<<endl;
      throw(runtime_error(args[0]));
    }
    int forstage=atoi(args[0].substr(6,1).c_str());
    found=checkArgs(args,0);
    checkStage(args,forstage);
    // if(getCommandFileLength(forstage)>0 && hpc>0){
    //   std::ostringstream ss;
    //   ss<<"-writep exists stage"<<stage;
    //   throw(runtime_error(ss.str()));
    // }
  }
  
  if(args[0].compare("-makes1")==0){ /// Create STAGE 2 COMMANDS
    if(recombfiles.size()==0) {
      cerr<<"ERROR: stage1 not required for unlinked model, but no recombination files provided. Did you forget to provide the recombination map?"<<endl;
      throw(runtime_error(args[0]));
    }
    makeStageCpEM(1);
  }

  if(args[0].compare("-makes2")==0){ /// Create STAGE 2 COMMANDS
    if(recombfiles.size()>0) linkagemode=string("linked");
    makeStageCpRun(2);
  }

  if(args[0].compare("-makes3")==0) { // finestructure MCMC
  //
    if(args.size()>1) istringstream ( args[1] ) >> numskip;
    if(args.size()>2) istringstream ( args[2] ) >> s3iters;
    makefsmcmc();
  }

  if(args[0].compare("-makes4")==0) { // finestructure TREE
  //
    if(verbose) cout<<"Making tree in stage "<<stage<<endl;
    if(args.size()==2) istringstream ( args[1] ) >> s4iters;
    makefstree();
  }

  if(args[0].compare("-makes6")==0){ /// Create STAGE 6 COMMANDS
   if(recombfiles.size()==0) {
      cerr<<"ERROR: stage6 not required for unlinked model, but no recombination files provided. Did you forget to provide the recombination map?"<<endl;
      throw(runtime_error(args[0]));
    }
    //    cerr<<"DEBUG: making cpem(6) found ="<<found<<endl;
    makeStageCpEM(6);
  }

  if(args[0].compare("-makes7")==0){ /// Create STAGE 7 COMMANDS
    if(recombfiles.size()>0) linkagemode=string("linked");
    makeStageCpRun(7);
  }

  if(args[0].compare("-makes9")==0){ // Make the GLOBETROTTER output files
    found=checkArgs(args,0);
    makeGt();
    stage=9;
  }

  ///////////////////////////////////
  // Doing commands

  if(args[0].compare("-gt")==0){ // do everything for GT
    // This is where we explain what to run and stuff
    found=checkArgs(args,0);
    if(hpc==0) hpc=1;
    cout<<"GLOBETROTTER: Explanation here."<<endl;
    std::ostringstream ss;
    ss<<"-writes missing stage"<<9;
    throw(runtime_error(ss.str()));
  }

  
  if(args[0].compare("-go")==0){ // do everything! Resolve by dependencies
    if(popidfile.compare("")==0) {
      found=checkArgs(args,0);
      cout<<"Finestructure complete!"<<endl;
      cout<<"Get started by running the GUI with:"<<endl<<"\"finegui -c " <<cpchunkcounts<<" -m "<<fsmcmcoutputvec[0]<<" -t "<<fstreeoutputvec[0];
      if(fsmcmcoutput.size()>1) cout<<" -m2 "<<fsmcmcoutputvec[1]<<" -t2 "<<fstreeoutputvec[1];
      cout<<"\""<<endl;
      cout<<"Then click \"File->Open\", then \"Read Data File\", \"Read Pairwise Coincidence\" and \"Read Tree\". Then you can explore the results."<<endl;
      if(fsmcmcoutput.size()>1) cout<<"Check convergence results by click \"File->Manage Second Dataset\", then \"Read Data File\", \"Read Pairwise Coincidence\" and \"Read Tree\". Then close the window and \"View->Pairwise Coincidence\", then \"Second view->Enable Alternative Diagonal View\" and \"Second view->Use second dataset\", then finally \"Second view->Pairwise Coincidence\". The top right diagonal shows the second MCMC run."<<endl;
      cout<<"You may also perform steps for an admixture analysis by either:\n\
     a) Specifying an idfile containing population information with -popidfile <file> -go;\n\
     b) Creating this automatically from the FineSTRUCTURE run (currently incomplete, see fs fs -h for -popidfile for how to extract this manually).\n"<<endl;
    }else{
      found=checkArgs(args,0);
      cout<<"ChromoPainter Admixture analysis ready. This has generated valid palettes from which admixture estimation via a mixture model is valid. Stage8  will become this admixture analysis in the future."<<endl;
    }
  }
  if(args[0].compare("-dos1")==0){
    found=checkArgs(args,0);
    checkStage(args,1);
    if(verbose) cout<<"RUNNING "<<s1commands.size()<<" stage1 command(s)!"<<endl;
    doCpStage(1);
  }
  if(args[0].compare("-dos2")==0){ 
    found=checkArgs(args,0);
    checkStage(args,2);
    if(verbose) cout<<"RUNNING "<<s2commands.size()<<" stage2 command(s)!"<<endl;
    doCpStage(2);
  }


  ///////////////////////////////////
  // Combining commands
  if(args[0].substr(0,9).compare("-combines")==0){ /// Check MAKE commands
    if(args[0].size()<=9){
      cerr<<"ERROR: -combines provided without any stage specified"<<endl;
      throw(runtime_error(args[0]));
    }
    int forstage=atoi(args[0].substr(9,1).c_str());
    found=checkArgs(args,0,1);
    checkStage(args,forstage);
  }
  
  if(args[0].compare("-combines1")==0){ // EM parameter estimates combination
    combineCpEm(1);
  }

  if(args[0].compare("-combines2")==0){ // chromopainter main run combination
    if(args.size()>1) cproot=args[1];
    combineStageCpRun(2);
  }

  if(args[0].compare("-combines3")==0) { // finestructure MCMC validation
    combineStage3();
    if(!validatedoutput[3]){
      if(++fscounter <= fsmaxattempts){
	cout<<"WARNING: Stage 3 convergence criterion failed! Use \"-ignoreGR\" to continue regardless. (Set the parameter \"-threshGR:-1\" to ignore the GR statistic in future.)  Re-running MCMC for longer: this is attempt "<<fscounter<<" of "<<fsmaxattempts<<"."<<endl;
	continuefsmcmc();
      }else{
	cout<<"WARNING: Stage 3 convergence criterion failed! Set -threshGR -1 to ignore the GR statistic. Continuing to tree inference due to exceeding maximum number of attempts."<<endl;
      }
    }
  }

  if(args[0].compare("-combines4")==0) { // finestructure TREE validation
    combineStage4();
  }

  if(args[0].compare("-combines5")==0){
    found=checkArgs(args,0);
    makePopIds();
    found=true;
  }
  
  if(args[0].compare("-combines6")==0){ // EM parameter estimates combination
    combineCpEm(6);
  }

  if(args[0].compare("-combines7")==0){ // chromopainter main run combination
    if(args.size()>1) popcproot=args[1];
    else ensureCpRunARoot(7);
    combineStageCpRun(7);
  }

  ///////////////////////////////////
  // Combining commands
  if(args[0].substr(0,8).compare("-remakes")==0){ /// Check MAKE commands
    int forstage;
    if(args.size()>1 && args[0].size()<=8){
      forstage=atoi(args[1].c_str());
    }else if(args[0].size()>8){
      forstage=atoi(args[0].substr(8,1).c_str());
    }else{
      cerr<<"ERROR: -remakes provided without any stage specified"<<endl;
      throw(runtime_error(args[0]));
    }

    if(verbose) cout<<"Processing -remakes"<<forstage<<endl;
    found=checkArgs(args,0,1);
    checkStage(args,forstage);
    if(args.size()==2) {
      forstage=atoi(args[1].c_str());
    }
    std::stringstream makecmd, writecmd;
    makecmd<<"-makes"<<forstage;
    if(verbose) cout<<"Remaking commands ("<<makecmd.str()<<")"<<endl;
    docmd(makecmd.str());
    if(verbose) cout<<"Removing incomplete commands."<<endl;
    removeCompletedCommands(forstage);
    safeRemoveCommandFile(forstage);
    if(hpc==1){
      writecmd<<"-writes"<<forstage;
      if(verbose) cout<<"Rewriting commands ("<<writecmd.str()<<")"<<endl;
      docmd(writecmd.str());
    }else{
      writecmd<<"-dos"<<forstage;
      if(verbose) cout<<"Rewriting commands ("<<writecmd.str()<<")"<<endl;
      docmd(writecmd.str());
    }
  }

  
  ///////////////////////////////////
  // Writing command for hpc commands
  if(args[0].substr(0,7).compare("-writes")==0){ /// Check MAKE commands
    if(args[0].size()<=7){
      cerr<<"ERROR: -writes provided without any stage specified"<<endl;
      throw(runtime_error(args[0]));
    }
    int forstage=atoi(args[0].substr(7,1).c_str());
    found=checkArgs(args,0,1);
    checkStage(args,forstage);
    ensureCommandfile(forstage);
    if(getCommandFileLength(forstage)>0 && hpc>0){
      std::ostringstream ss;
      ss<<"-writes exists stage"<<stage;
      throw(runtime_error(ss.str()));
    }
  }

  if(args[0].compare("-writes1")==0){
    if(args.size()==2) s1commandfile=args[1];
    if(verbose) cout<<"Writing "<<s1commands.size()<<" stage 1 commands to "<<s1commandfile<<endl;
    bool addfs=!paintercmd.compare("cp");
    writeStringVectorToFile(s1commands,s1commandfile,s1logfiles,addfs);
  }

  if(args[0].compare("-writes2")==0){
    if(args.size()==2) s2commandfile=args[1];
    if(verbose) cout<<"Writing "<<s2commands.size()<<" stage 2 commands to "<<s2commandfile<<endl;
    bool addfs=!paintercmd.compare("cp");
    writeStringVectorToFile(s2commands,s2commandfile,s2logfiles,addfs);
  }

  if(args[0].compare("-writes3")==0){
    if(args.size()==2) s3commandfile=args[1];
    if(verbose) cout<<"Writing "<<s3commands.size()<<" stage 3 commands to "<<s3commandfile<<endl;
    writeStringVectorToFile(s3commands,s3commandfile,s3logfiles);
  }

  if(args[0].compare("-hpcs3")==0){
    found=checkArgs(args,0,1);
    checkStage(args,3);
    if(args.size()==2) s3commandfile=args[1];
    else ensureCommandfile(3);
    if(verbose) cout<<"Writing a single, parallel stage 3 command to "<<s3commandfile<<endl;
    writeHpcStage3(s3commandfile);
    writeStringVectorToFile(s3commands,s3commandfile,s3logfiles);
  }

  if(args[0].compare("-writes4")==0){
    if(args.size()==2) s4commandfile=args[1];
    if(verbose) cout<<"Writing "<<s4commands.size()<<" stage 4 commands to "<<s4commandfile<<endl;
    writeStringVectorToFile(s4commands,s4commandfile,s4logfiles);
  }

  // fs2.1
  if(args[0].compare("-writes6")==0){
    if(args.size()==2) s6commandfile=args[1];
    if(verbose) cout<<"Writing "<<s6commands.size()<<" stage 6 commands to "<<s6commandfile<<endl;
    bool addfs=!paintercmd.compare("cp");
    writeStringVectorToFile(s6commands,s6commandfile,s6logfiles,addfs);
  }

  if(args[0].compare("-writes7")==0){
    if(args.size()==2) s7commandfile=args[1];
    if(verbose) cout<<"Writing "<<s7commands.size()<<" stage 7 commands to "<<s7commandfile<<endl;
    bool addfs=!paintercmd.compare("cp");
    writeStringVectorToFile(s7commands,s7commandfile,s7logfiles,addfs);
  }

  if(args[0].compare("-writes9")==0){ // Write the GLOBETROTTER output files
    if(args.size()==2) s9commandfile=args[1];
    if(verbose) cout<<"Writing "<<s9commands.size()<<" stage 9 (GT) commands to "<<s9commandfile<<endl;
    writeStringVectorToFile(s9commands,s9commandfile,s9logfiles,false);
  }

  ///////////////////////////////////
  // Properties of the data
  if(args[0].compare("-haploid")==0){
    found=checkArgs(args,0);
    checkStage(args,1);
    ploidy=1;
    if(verbose) cout<<"Setting Haploid mode"<<endl;
  }

  ///////////////////////////////////
  // FINESTRUCTURE commands
  if(args[0].compare("-configmcmc")==0) { // set mcmc parameters
  //
    found=checkArgs(args,3);
    checkStage(args,3);
    s3iters=-1;
    istringstream ( args[1] ) >> s3itersburnin;
    istringstream ( args[2] ) >> s3iterssample;
    istringstream ( args[3] ) >> numskip;
  }
  if(args[0].compare("-ignoreGR")==0) { // set mcmc parameters
  //
    found=checkArgs(args,0);
    checkStage(args,3);
    ignoreGRfsmcmc();
  }
  if(args[0].compare("-dos3")==0) { // finestructure MCMC
    found=checkArgs(args,0);
    checkStage(args,3);
    dofsmcmc();
  }
  if(args[0].compare("-dos4")==0) { // finestructure TREE
    found=checkArgs(args,0);
    checkStage(args,4);
    dofstree();
  }
  ////////////////////////////////////////////
  // fs 2.1
  if(args[0].compare("-dos5")==0) { // extract populations from FineSTRUCTURE
    found=checkArgs(args,0);
    checkStage(args,5);
    doStagePopExtract();
  }
  if(args[0].compare("-dos6")==0){
    found=checkArgs(args,0);
    checkStage(args,6);
    if(verbose) cout<<"RUNNING "<<s6commands.size()<<" stage6 command(s)!"<<endl;
    doCpStage(6);
  }
  if(args[0].compare("-dos7")==0){ 
    found=checkArgs(args,0);
    checkStage(args,7);
    if(verbose) cout<<"RUNNING "<<s7commands.size()<<" stage7 command(s)!"<<endl;
    doCpStage(7);
  }

  
  ///////////////////////////////////
  // Catch all others

  if(!found){
    try{
      //      if(verbose) cout<<"DEBUG: Searching for parameter "<<args[0]<<endl;
      found=applyVal(args);
      //      cerr<<"DEBUG: found was "<<found<<endl;
    }catch (exception& e)  {
      throw(e);
    }
    if(!found){
      //      cerr<<"DEBUG: ERROR: Command or parameter "<<args[0]<<" not recognised"<<endl;
      throw(runtime_error(args[0]));
    }
  }

  
  /// Update the history
  std::ostringstream newhist;
  newhist<<"CMD:"<<args[0];
  for(unsigned int i=1;i<args.size();i++){
    newhist<<" "<<args[i];
  }
  newhist<<endl;
  historytext.append(newhist.str());

  // Stop if we need things to be run
  char * tstr_c=(char*) args[0].data();
  if(hpc==1 && strstr (tstr_c,"-writes")!=NULL){
    std::ostringstream ss;
    ss<<"-writes missing stage"<<stage;
    throw(runtime_error(ss.str()));
  }

  // Update the project file
  if(fsmode>0){
    safeCreateFile();
    if(verbose) cout<<"Writing project file..."<<endl;
    writeToFile();
  }
}

void FsProject::matchCpInputs(){
  if(phasefiles.size()==0) throw(runtime_error("fsproject: No phase files defined!"));

  if(linkagemode.compare("linked")!=0) {
    // we have the unlinked model and we don't need recombination data
    return;
  }
  // else linked mode
  if(phasefiles.size()!=recombfiles.size()) throw(runtime_error("fsproject: Different number of phase to recombination files!"));

  return;
}

int FsProject::getCommandFileLength(int forstage){
  if(forstage<0) forstage=stage;
  switch(forstage){
  case 1: return(countLines(s1commandfile));
  case 2: return(countLines(s2commandfile));
  case 3: return(countLines(s3commandfile));
  case 4: return(countLines(s4commandfile));
  case 6: return(countLines(s6commandfile));
  case 7: return(countLines(s7commandfile));
  case 9: return(countLines(s9commandfile));
  default:
    cerr<<"ERROR: getCommandfileLength received an invalid option!"<<endl;
    throw(logic_error("haveCommandfileContent"));
  }
}


void FsProject::ensureCommandfile(int forstage)
{
  string ts;
  switch(forstage){
  case 1: ts=s1commandfile;break;
  case 2: ts=s2commandfile;break;
  case 3: ts=s3commandfile;break;
  case 4: ts=s4commandfile;break;
  case 6: ts=s6commandfile;break;
  case 7: ts=s7commandfile;break;
  case 9: ts=s9commandfile;break;
  default: 
    cerr<<"ERROR: ensureCommandfile received an invalid option!"<<endl;
    throw(logic_error("ensureCommandfile"));
  }
  if(ts.compare("")==0){
    ostringstream ss;
    ss<<dirname<<"/commandfiles";
    ensureDirectory(ss.str());
    ss<<"/commandfile"<<forstage<<".txt";
    switch(forstage){
    case 1: s1commandfile=ss.str();break;
    case 2: s2commandfile=ss.str();break;
    case 3: s3commandfile=ss.str();break;
    case 4: s4commandfile=ss.str();break;
    case 6: s6commandfile=ss.str();break;
    case 7: s7commandfile=ss.str();break;
    case 9: s9commandfile=ss.str();break;
    default: 
      throw(logic_error("ensureCommandfile"));
    }
  }
}

void FsProject::ensureCpEmRoot(int forstage){
  string thisroot;
  
  if(forstage==1) thisroot=s1outputroot;
  else thisroot=s6outputroot;
  
  if(thisroot.compare("")==0){
    ostringstream ss;
    ss<<fileroot<<"_stage"<<forstage<<"_tmp_EM_linked";
    if(ploidy==1)  ss<<"_haploid";
    thisroot=string(ss.str());
  }

  if(forstage==1) s1outputroot=thisroot;
  else s6outputroot=thisroot;
}

void FsProject::ensureCpRunRoot(int forstage){
  string thisroot;
  
  if(forstage==2) thisroot=s2outputroot;
  else thisroot=s7outputroot;

  if(thisroot.compare("")==0) {
    ostringstream ss;
    ss<<fileroot<<"_stage"<<forstage<<"_tmp_mainrun."<<linkagemode;
    if(ploidy==1)  ss<<"_haploid";
    thisroot=string(ss.str());
  }

  if(forstage==2) s2outputroot=thisroot;
  else s7outputroot=thisroot;
}

void FsProject::ensureCpRunARoot(int forstage){
  string thisroot;
  int chunks;
  chunks=s2chunksperregion;
  
  if(forstage==2) {
    thisroot=cproot;
  }else {
    thisroot=popcproot;
  }

  if(thisroot.compare("")==0){
    ostringstream ss;
    ss<<fileroot;
    if(forstage==7) ss<<"_pop";
    ss<<"_"<<linkagemode;
    if(chunks>0)  ss<<"_cpr"<<chunks;
    if(ploidy==1)  ss<<"_hap";
    thisroot=ss.str();
  }

  if(forstage==2) cproot=thisroot;
  else popcproot=thisroot;
}

void FsProject::ensurefsroot(){
  if(fsroot.compare("")==0){
    ostringstream ss;
    ss<<fileroot<<"_"<<linkagemode;
    if(s2chunksperregion>0)  ss<<"_cpr"<<s2chunksperregion;
    if(ploidy==1)  ss<<"_hap";
    fsroot=ss.str();
  }
}

void FsProject::ensurefsmcmc(){
  ensurefsroot();
  if(fsmcmcoutput.compare("")==0){
    ostringstream ss;
    ss<<dirname<<"/stage3";
    ensureDirectory(ss.str());
    ss<<"/"<<baseName(fsroot)<<"_mcmc";
    fsmcmcoutput=ss.str();
  }
  fsmcmcoutputvec.clear();
  for(int c1=0;c1<nummcmcruns;c1++){
    ostringstream ss;
    ss<<fsmcmcoutput;
    if(c1>0) ss<<"_run"<<c1;
    ss<<".xml";
    fsmcmcoutputvec.push_back(ss.str());
  }
}


void FsProject::ensurefstree(){
  ensurefsroot();
  if(fstreeoutput.compare("")==0){
    ostringstream ss;
    ss<<dirname<<"/stage4";
    ensureDirectory(ss.str());
    ss<<"/"<<baseName(fsroot)<<"_tree";
    fstreeoutput=ss.str();
  }
  fstreeoutputvec.clear();
  for(int c1=0;c1<nummcmcruns;c1++){
    ostringstream ss;
    ss<<fstreeoutput;
    if(c1>0) ss<<"_run"<<c1;
    ss<<".xml";
    fstreeoutputvec.push_back(ss.str());
  }
}

std::string FsProject::makeoutfileroot(string root,int fileon,int indstart,int indend,int forstage){
  std::ostringstream ss;
  ss<<dirname<<"/"<<"stage"<<forstage<<"/";
  ensureDirectory(ss.str());

  if((indstart==0) & (indend==0)) ss<<baseName(root)<<"_file"<<fileon+1<<"_allinds";
  else if (indstart==indend){
    ss<<baseName(root)<<"_file"<<fileon+1<<"_ind"<<indstart;
  }else{
    ss<<baseName(root)<<"_file"<<fileon+1<<"_ind"<<indstart<<"-"<<indend;
  }
  return(ss.str());
}

int FsProject::getNsnpsFromFile(unsigned int fileno){
  if(fileno>=phasefiles.size()) {
    cerr<<"ERROR: Need to have phase file "<< fileno<<" to count the number of SNPs in it!"<<endl;
    throw(logic_error("Asked to count SNPs of a file we don't have"));
  }
  vector<int> cpv=getChromoPainterHeaderInfo(phasefiles[fileno],ploidy);
  return(cpv[1]);
}

int FsProject::getNhapsFromFile(int which){

  // Otherwise get it from the phase file
  if((int)phasefiles.size()<which){
    cerr<<"ERROR: Need to have phase file to count the number of individiuals!"<<endl;
    throw(runtime_error("data missing"));
  }
  vector<int> cpv=getChromoPainterHeaderInfo(phasefiles[which],ploidy);
  return(cpv[0]);
}

bool FsProject::readmcmctraces(int filenum){
  string file=fsmcmcoutputvec[filenum];
  FsXml *fs=new FsXml(file);
  streampos fpos=fs->gotoLineContaining("<Iteration>");
  double t_posterior,t_beta,t_delta,t_f,t_k;
  int counts=0;
  
  
  mcmc_posterior.push_back(std::vector<double>());
  mcmc_k.push_back(std::vector<double>());
  mcmc_delta.push_back(std::vector<double>());
  mcmc_beta.push_back(std::vector<double>());
  mcmc_f.push_back(std::vector<double>());

  while(!fs->eof() && fpos>=0) {
    counts++;
    if(fpos>0) {
      string s_pos=fs->getParam("Posterior",fpos);
      string s_K=fs->getParam("K",fpos);
      string s_Beta=fs->getParam("beta",fpos);
      string s_Delta=fs->getParam("delta",fpos);
      string s_F=fs->getParam("F",fpos);
      istringstream ( s_pos ) >> t_posterior;
      istringstream ( s_K ) >> t_k;
      istringstream ( s_Beta ) >> t_beta;
      istringstream ( s_Delta ) >> t_delta;
      istringstream ( s_F ) >> t_f;
      
      mcmc_posterior[filenum].push_back(t_posterior);
      mcmc_k[filenum].push_back(t_k);
      mcmc_beta[filenum].push_back(log(t_beta));
      mcmc_delta[filenum].push_back(t_delta);
      mcmc_f[filenum].push_back(t_f);
      fpos=fs->gotoNextLineContaining("<Iteration>");
    }
  }
  delete(fs);
  if(verbose) cout<<"Read "<<mcmc_posterior[filenum].size()<<" iterations."<<endl;
  return(1);
}

double FsProject::mcmcGRstatistic(std::vector<std::vector<double> > *data){
  double nruns=(double)data->size();
  vector<double> runsums;
  vector<double> runsumsquares;
  double ntot=0;
  vector<double> runn;
  double combinedsum=0;
  double combinedsumsquares=0;
  // construct the sums and sums of squares that are needed
  for(unsigned int c1=0;c1<data->size();++c1){
    runsums.push_back(0);
    runsumsquares.push_back(0);
    runn.push_back(data->at(c1).size());
    ntot+=data->at(c1).size();
    for(unsigned int c2=0;c2<data->at(c1).size();++c2){
      combinedsum+=data->at(c1)[c2];
      runsums[c1]+=data->at(c1)[c2];
      combinedsumsquares+=data->at(c1)[c2]*data->at(c1)[c2];
      runsumsquares[c1]+=data->at(c1)[c2]*data->at(c1)[c2];
    }
  }
  // compute W, the within chain variance
  double totvar= combinedsumsquares/ntot - (combinedsum/ntot)*(combinedsum/ntot);
  totvar*=ntot/(ntot-1);
  vector<double> runvarest;
  double W=0;
  for(unsigned int c1=0;c1<data->size();++c1){
    double tv=(runsumsquares[c1]/runn[c1] - (runsums[c1]/runn[c1])*runsums[c1]/runn[c1]);
    runvarest.push_back(tv * runn[c1]/(runn[c1]-1));
    W += runvarest[c1]/nruns;
  }
  //  cout<<"GR Calc W="<<W <<" totvar="<<totvar<<endl;
  // compute B, the between chain variance
  double B=0;
  for(unsigned int c1=0;c1<data->size();++c1){
    B+= runn[c1]/(nruns-1) *(runsums[c1]/runn[c1] - combinedsum/ntot)*(runsums[c1]/runn[c1] - combinedsum/ntot);
  }
  //  cout<<"GR Calc B="<<B<<endl;
  // Compute the variance estimator for the combined chains
  double nbar=ntot/nruns;
  double sigma_hatsq=(nbar-1)*W/nbar + B/nbar;
  double Rhat = sqrt(sigma_hatsq/W);
  //  cout<<"GR Calc sigmahatsq="<<sigma_hatsq<<" Rhat="<<Rhat<<endl;
  if(sigma_hatsq==0 && W==0) Rhat=1;
  return(Rhat);
}

bool FsProject::mcmcConvergence(){
  if(fsmcmcoutputvec.size()==0){
    cerr<<"WARNING: mcmcConvergence: Haven't yet run the MCMC! Shouldn't be checking it at this stage."<<endl;
    return(1);
  }else if(fsmcmcoutputvec.size()==1){
    cerr<<"WARNING: mcmcConvergence: have only run one MCMC chain. Cannot check convergence."<<endl;
    return(1);
  }

  // Read the posterior traces
  mcmc_posterior.clear();
  mcmc_k.clear();
  mcmc_beta.clear();
  mcmc_delta.clear();
  mcmc_f.clear();

  for(unsigned int i=0;i<fsmcmcoutputvec.size();++i){
    if(verbose) cout<<"Reading MCMC traces for file "<<i+1<<" of "<<fsmcmcoutputvec.size()<<endl;
    if(!readmcmctraces(i)){
      cerr<<"ERROR: Cannot read MCMC traces for file "<<fsmcmcoutputvec[i]<<endl;
    }
  }
  
  stringstream mcmctracesfile;
  mcmctracesfile<<fsmcmcoutput<<".mcmctraces.tab";
  writemcmctraces(mcmctracesfile.str());

  // Compute the Gelman Rubin diagnostics
  mcmcGR.clear();
  mcmcGR.push_back(mcmcGRstatistic(&mcmc_posterior));
  mcmcGR.push_back(mcmcGRstatistic(&mcmc_k));
  mcmcGR.push_back(mcmcGRstatistic(&mcmc_beta));
  mcmcGR.push_back(mcmcGRstatistic(&mcmc_delta));
  mcmcGR.push_back(mcmcGRstatistic(&mcmc_f));
  if(verbose)  {
    cout<<"Gelman-Rubin statistics:";
    cout<<" GR(Log Posterior) = "<<mcmcGR[0];
    cout<<" GR(Number of populations K) = "<<mcmcGR[1];
    cout<<" GR(Log Beta) = "<<mcmcGR[2];
    cout<<" GR(Delta) = "<<mcmcGR[3];
    cout<<" GR(f) = "<<mcmcGR[4]<<endl;
  }
  if(threshGR<0) return(1);
  for(unsigned int i=0;i<mcmcGR.size();++i) { 
    if(mcmcGR[i]>threshGR) return(0);
  }
  return(1);
}

bool FsProject::writemcmctraces(std::string filename){
  ofstream ofile;
  ofile.open (filename.c_str());

  if(mcmc_posterior.size()==0) return(0);
  for(unsigned int c2=0;c2<mcmc_posterior.size();c2++){
    ofile<<"LogPosterior"<<c2<<" K"<<c2<<" beta"<<c2<<" delta"<<c2<<" f"<<c2<<" ";
  }
  ofile<<endl;
  for(unsigned int c1=0;c1<mcmc_posterior[0].size();c1++){
    for(unsigned int c2=0;c2<mcmc_posterior.size();c2++){
      ofile<<mcmc_posterior[c2][c1]<<" ";
      ofile<<mcmc_k[c2][c1]<<" ";
      ofile<<mcmc_beta[c2][c1]<<" ";
      ofile<<mcmc_delta[c2][c1]<<" ";
      ofile<<mcmc_f[c2][c1]<<" ";
    }
    ofile<<endl;
  }
  ofile.close();
  return(1);
}

void FsProject::ignoreGRfsmcmc(){
  if(stage>=3){ // If we are supposed to have output
    cout<<"INFO: Restoring previous MCMC results that failed the GR statistic test."<<endl;
    s3iterssample/=2;
    s3itersburnin=s3iterssample;
    numskip/=2;

    for(unsigned int c1=0;c1<old_fsmcmcoutputvec.size();++c1){
      if(access( filename.c_str(), F_OK ) == -1 ) { // file does not exist
	cerr<<"ERROR: Cannot find original MCMC file. Was the MCMC actually run, and did we reject it on the basis on Gelman-Rubin statistics? If not, try \"-threshGR:-1 -go\" instead of \"-ignoreGR\"."<<endl;
	throw(runtime_error("-ignoreGR cannot reconstruct original MCMC filename"));
      }
      if(verbose) cout<<"Renaming "<<old_fsmcmcoutputvec[c1].c_str()<<" to "<<fsmcmcoutputvec[c1].c_str()<<endl;
      if(rename(old_fsmcmcoutputvec[c1].c_str(),fsmcmcoutputvec[c1].c_str())!=0){
	cerr<<"WARNING: -ignoreGR cannot rename "<<old_fsmcmcoutputvec[c1].c_str()<<" to "<<fsmcmcoutputvec[c1].c_str()<<". Rerunning MCMC instead!"<<endl;
      }
    }
  }else{
      cout<<"INFO: Ignoring the Gelman-Rubin statistic for MCMC."<<endl;
  }

  // Things we do regardless
  threshGR=-1;
  old_fsmcmcoutputvec.clear();
}

void FsProject::continuefsmcmc(){
  cout<<"INFO: Doubling sample time and continuing from previous run"<<endl;
  old_fsmcmcoutputvec.clear();
  for(unsigned int c1=0;c1<fsmcmcoutputvec.size();++c1){
    size_t found =fsmcmcoutputvec[c1].find(".xml");
    if(found==string::npos || access( filename.c_str(), F_OK ) == -1 ) { // invalid filename, or file does not exist
      cerr<<"ERROR: Cannot reconstruct original MCMC filename from the parameter file. Was the MCMC actually run, and did we reject it on the basis on Gelman-Rubin statistics? If not, try \"-threshGR:-1 -go\" instead of \"-ignoreGR\"."<<endl;
      throw(runtime_error("continuefsmcmc cannot reconstruct original MCMC filename"));
    }
    string tstr=fsmcmcoutputvec[c1].substr(0,found);
    ostringstream ss;
    ss<<tstr<<"_x"<<s3itersburnin<<"_y"<<s3iterssample<<"_z"<<numskip<<".xml";
    old_fsmcmcoutputvec.push_back(ss.str());
    rename(fsmcmcoutputvec[c1].c_str(),old_fsmcmcoutputvec[c1].c_str());
  }
  fsmcmcoutput.clear();
  fsmcmcoutputvec.clear();

  if(s3iters>0){
    s3itersburnin=0;
    s3iterssample=s3iters;
    s3iters=-1;
  }else{
    s3itersburnin=0;
    s3iterssample*=2;
  }
  numskip*=2;

  
  vector<string> redocmd;
  redocmd.push_back("-makes3");
  docmd(redocmd);
  redocmd.clear();
  if(hpc) redocmd.push_back("-writes3");
  if(hpc!=1) redocmd.push_back("-dos3");
  docmd(redocmd);
  if(hpc==1) throw(runtime_error("fsproject: combines3 stage3 missing results"));
  redocmd.clear();
  redocmd.push_back("-combines3");
  docmd(redocmd);
}

int FsProject::getNindsFromFile(bool keepall,string idf){
  if(idf.compare("")==0) idf=idfile;
  if(idf.compare("")==0) {
    cerr<<"ERROR: An idfile is required! You can create it with -createid or read an existing one with -idfile ."<<endl;
    throw(logic_error("idfile missing"));
  }
  std::vector<std::string> tnames=getIdsFromFile(idf,keepall);
  return((int) tnames.size());
}

int FsProject::getUniqueNindsFromFile(bool keepall,string idf){
  if(idf.compare("")==0) idf=idfile;
  if(idf.compare("")==0) {
    cerr<<"ERROR: An idfile is required! You can create it with -createid or read an existing one with -idfile ."<<endl;
    throw(logic_error("idfile missing"));
  }
  std::vector<std::string> tnames=unique(getIdsFromFile(idf,keepall));
  return((int) tnames.size());
}

void FsProject::makeStageCpEM(int forstage)
{
  if(linkagemode.compare("linked")!=0) {
    cerr<<"WARNING: Requested CP inference for an unlinked model!"<<endl;
    return; ///< Nothing to do 
  }
  if(nsnpsvec.size()!=phasefiles.size()){
    cerr<<"ERROR: Number of snps not corrected counted? Have "<<phasefiles.size()<< " phase files but only "<<nsnpsvec.size()<<" counts of SNPs"<<endl;
    throw(logic_error("makeStageCpEM nsnpsvec problem"));
  }

  matchCpInputs(); ///< checks all is well
  ensureCpEmRoot(forstage);///< construct the output file name root

  ////////////////////////////////
  vector<string> scommands,srootvec,slogfiles;
  string soutputroot;
  if(forstage==1) soutputroot=s1outputroot;
  else soutputroot=s6outputroot;

    int sXminsnps;
  double sXindfrac,sXsnpfrac;
  if(forstage<=2){
    sXindfrac=s1indfrac;
    sXsnpfrac=s1snpfrac;
    sXminsnps=s1minsnps;
  }else{
    sXindfrac=s6indfrac;
    sXsnpfrac=s6snpfrac;
    sXminsnps=s6minsnps;
  }

  /// Figure out which individuals go in which files
  int indsmax=numRecipients(forstage); 
  vector<int> indsvec,indsvecall; // a vector of individuals
  int tindsperproc=getIndsPerProc(forstage);
  int wantedinds=(int)(indsmax * sXindfrac);
  if(wantedinds>indsmax || wantedinds<=0){
    cerr<<"ERROR: Requested an invalid number ("<<wantedinds<<") of individuals due to -s"<<forstage<<"indfrac: 0<s"<<forstage<<"indfrac<=1, with at least one individual."<<endl;
    throw(runtime_error("makeStageCpEM sXindfrac problem"));
  }

  // Check if we need to drop to processing each ind separately, because we only process a subset of individuals
  if((wantedinds<indsmax) | (forstage!=1)){ // if we are only processing a subset
    tindsperproc=1;
    for(int indon=1;indon<=indsmax;indon++) indsvecall.push_back(indon);
  }
  
  // Check if we need to drop to processing each ind separately, because we only process a subset of snps
  vector<bool> usesubset;
  vector<int> wantedsnps;
  for(unsigned int fileon=0;fileon<phasefiles.size();fileon++) {
    int tsnps=nsnpsvec[fileon];
    if((sXsnpfrac>0)&&(sXminsnps>0)) tsnps=max((int)(nsnpsvec[fileon]*sXsnpfrac),sXminsnps);
    wantedsnps.push_back(tsnps);
    if(wantedsnps[fileon]<nsnpsvec[fileon]){ // we are allowed to do less processing
      tindsperproc=1;// must do each ind separately
      usesubset.push_back(true);
    }else usesubset.push_back(false);
  }

  //  cerr<<"DEBUG: Starting cmd construction"<<endl;
 // Construct the commands
  for(unsigned int fileon=0;fileon<phasefiles.size();fileon++) {
    if(indsvecall.size()>0){ // process a different random subset per file
      indsvec=sampleVec(indsvecall,wantedinds);
    }

    int sXemits;
    string sXargs,sXYargs;
    if(forstage<=2){
      sXargs=s1args;
      sXYargs=s12args;
      sXemits=s1emits;
    }else{
      sXargs=s6args;
      sXYargs=s67args;
      sXemits=s6emits;
    }
    /// Loop over individuals
    for(int tindon=1;tindon<=wantedinds;tindon+=tindsperproc){
      int indon=tindon;// tindon is the index of the individual we are processing; indon is the actual individual (in case we sample)
      std::ostringstream ss;
      int indonend=min(indon+tindsperproc-1,indsmax); // make sure we don't process past the end individual
      if(tindsperproc==0){// process all in one go
	indon=indonend=0;
      }else if(indsvecall.size()>0){ // process only a subset
	if(indon-1>=(int)indsvec.size()) throw(logic_error("makeStageCpEM: indsvec error, reached an individual that shouldn't exist"));
	indon=indonend=indsvec[indon-1];
      }
      ss<<paintercmd;
      string toutfile=makeoutfileroot(soutputroot,fileon,indon,indonend,forstage);
      if(ploidy==1) ss<<" -j";
      ss<<" -i "<<sXemits;
      if(usesubset[fileon]){
	ss<<" -e "<<wantedsnps[fileon];
      }
      if(sXargs.compare("")!=0) ss<<" "<<sXargs;
      if(sXYargs.compare("")!=0) ss<<" "<<sXYargs;
      if((idfile.compare("")!=0) && (forstage!=6)) ss<<" -t "<<idfile;
      if(forstage==6){
	ss<<" "<<popcmdargs[indon-1]<<" ";// indon starts at 1
      }else{
	ss<<" -a "<<indon<<" "<<indonend;
      }
      ss<<" -g "<<phasefiles[fileon]<<" -r "<<recombfiles[fileon]<<" -o "<<toutfile;

      srootvec.push_back(toutfile);
      scommands.push_back(ss.str());
      toutfile.append(".log");
      if(outputlogfiles) slogfiles.push_back(toutfile);

      if(verbose){
	cout<<"CREATING S"<<forstage<<" CMD:"<<ss.str()<<endl;
      }
      if(tindsperproc==0) tindon=wantedinds+1; // make sure we end in -a 0 0 mode
    }
  }

  if(forstage==1){
    s1commands=scommands;
    s1logfiles=slogfiles;
    s1outputrootvec=srootvec;
  }else{
    s6commands=scommands;
    s6logfiles=slogfiles;
    s6outputrootvec=srootvec; 
 }
}

void FsProject::makeStageCpRun(int forstage)
{
  matchCpInputs(); ///< checks all is well
  ensureCpRunRoot(forstage);///< construct the output file name root

  
  if(linkagemode.compare("linked")==0){// check for previous stage having been run
    if( (((Neinf<0)| (muinf<0))&&(forstage==2)) | (((popNeinf<0)| (popmuinf<0))&&(forstage==7))){
      cerr<<"Must have inferred Ne and mu when running in linkage mode. Has stage"<<forstage-1<<" been completed? Has -combines"<<forstage-1<<" been run?"<<endl;
      throw(runtime_error("fsproject: makeStageCpRun requires combine previous stage"));
    }
  }
  
  /// Figure out which individuals go in which files
  int tindsperproc=getIndsPerProc(forstage);

  vector<string> scommands,srootvec,slogfiles;
  string soutputroot;
  double useNe,usemu;
  if(forstage==2) {
    soutputroot=s2outputroot;
    useNe=Neinf;
    usemu=muinf;
  }else {
    soutputroot=s7outputroot;
    useNe=popNeinf;
    usemu=popmuinf;
  }

  int sXsamples,sXchunksperregion;
  string sXargs,sYXargs;
  if(forstage<=2){
    sXargs=s2args;
    sYXargs=s12args;
    sXsamples=s2samples;
    sXchunksperregion=s2chunksperregion;
  }else{
    sXargs=s7args;
    sYXargs=s67args;
    sXsamples=s7samples;
    sXchunksperregion=s7chunksperregion;
  }
  // Construct the commands
  int nindshere=numRecipients(forstage);
  if(verbose) cout<<"Preparing command lines for "<<nindshere<<" recipients."<<endl;
  for(unsigned int fileon=0;fileon<phasefiles.size();fileon++){
    for(int indon=1;indon<=nindshere;indon+=tindsperproc){
      std::ostringstream ss;
      int indonend=min(indon+tindsperproc-1,nindshere); // make sure we don't process past the end individual
      ss<<paintercmd<<" ";
      string toutfile=makeoutfileroot(soutputroot,fileon,indon,indonend,forstage);
      if(ploidy==1) ss<<" -j";
      if(sXargs.compare("")!=0) ss<<" "<<sXargs<<" "; 
      if(sYXargs.compare("")!=0) ss<<" "<<sYXargs<<" ";
      if((idfile.compare("")!=0) & (forstage!=7)) ss<<" -t "<<idfile<<" ";
      
      if(forstage==7){
	ss<<popcmdargs[indon-1]<<" "; // indon starts at 1
      }else{
	ss<<" -a "<<indon<<" "<<indonend;
      }

      if(linkagemode.compare("linked")==0) {
	ss<<" -r "<<recombfiles[fileon]<<" -n "<<useNe<<" -M "<<usemu;
      }else {
	ss<<" -u";
      }
      if(sXsamples>0) ss<<" -s "<<sXsamples;
      if(sXchunksperregion>=0) ss<<" -k "<<sXchunksperregion;
      else ss<<" -k "<<defaultChunksperregion();
      ss<<" -g "<<phasefiles[fileon];
      ss<<" -o "<<toutfile;
      srootvec.push_back(toutfile);
      scommands.push_back(ss.str());
      toutfile.append(".log");
      if(outputlogfiles) slogfiles.push_back(toutfile);
      if(verbose){
	cout<<"CREATING S"<<forstage<<" CMD:"<<ss.str()<<endl;
      }
    }
  }
  stage=forstage;
  if(forstage==2){
    s2commands=scommands;
    s2logfiles=slogfiles;
    s2outputrootvec=srootvec;
  }else{
    s7commands=scommands;
    s7logfiles=slogfiles;
    s7outputrootvec=srootvec; 
 }

}


void FsProject::makefsmcmc(){
  if(s3iters>0){
    s3iterssample=s3iters/2;
    s3itersburnin=s3iters/2;
  }
  if(numskip<0){
    numskip = max(1,(int)floor(s3iterssample/maxretained));
  }
  ensurefsmcmc();
  s3commands.clear();
  for(int runon=0;runon<nummcmcruns;runon++){
    std::ostringstream ss;
    ss<<"fs -s "<<runon*2+1;
    if(s34args.compare("")!=0) ss<<" "<<s34args;
    ss<<" -x "<<s3itersburnin<<" -y "<<s3iterssample<<" -z "<<numskip<<" "<<cpchunkcounts<<" ";
    if(old_fsmcmcoutputvec.size()==fsmcmcoutputvec.size()) ss<<old_fsmcmcoutputvec[runon]<<" ";
    ss<<fsmcmcoutputvec[runon];
    s3commands.push_back(ss.str());
    string toutfile=fsmcmcoutputvec[runon];
    toutfile.append(".log");
    if(outputlogfiles) s3logfiles.push_back(toutfile);

    if(verbose){
      cout<<"CREATING S3 CMD:"<<ss.str()<<endl;
    }
  }
  stage=3;
}

void FsProject::makefstree(){
  ensurefstree();
  s4commands.clear();
  for(int runon=0;runon<nummcmcruns;runon++){
    std::ostringstream ss;
    ss<<"fs -m T -s "<<runon*2+1;
    if(s34args.compare("")!=0) ss<<" "<<s34args;
    ss<<" -x "<<s4iters<<" "<<s4args<<" "<<cpchunkcounts<<" "<<fsmcmcoutputvec[runon]<<" "<<fstreeoutputvec[runon];
    s4commands.push_back(ss.str());
    string toutfile=fstreeoutputvec[runon];
    toutfile.append(".log");
    if(outputlogfiles) s4logfiles.push_back(toutfile);
    if(verbose){
      cout<<"CREATING S4 CMD:"<<ss.str()<<endl;
    }
  }
}

void FsProject::dofsmcmc(){
  // Some sanity checks?
  // have we run chromocombine?
  if(cval<0) {
      cerr<<"FineStructure MCMC cannot be run without a valid 'c' value. Did chromocombine run? Did it calculate 'c'?"<<endl;
      throw(runtime_error("finestructure requires -combines2 to have completed successfully"));
  }
  // 
  int cmdon=0;
  if(numthreads>0) {
    do_omp_set_num_threads(numthreads);
  }
#pragma omp parallel for
  for(int runon=0;runon<nummcmcruns;runon++){
    int thread_number = get_omp_get_thread_num();

    string tsv=s3commands[runon];
    if(thread_number>0) tsv.append(" -S");
    std::vector<char *> argv=converttoargv(tsv);
    //    string logfile=fsmcmcoutputvec[runon];
    //    logfile.append(".log");
#pragma omp atomic
    cmdon++;

#pragma omp critical
    cout<<(ParallelStream()<<"Running stage 3 (mcmc) command number (~"<<cmdon<<") of "<<nummcmcruns).toString()<<endl;
    if(verbose) cout<<"RUNNING S3 CMD:"<<tsv<<endl;


    int rv=finestructure(argv.size(),argv.data());
    freeargv(argv);

    // Check that it ran correctly
    if(rv!=0){
      cerr<<"FineStructure MCMC failed!"<<endl;// See log file ("<<logfile<<") for details."<<endl;
      throw(runtime_error("finestructure"));
      }
  }// end for

}

void FsProject::dofstree(){
  // Some sanity checks?
  
  // 
  int cmdon=0;
  if(numthreads>0) {
    do_omp_set_num_threads(numthreads);
  }
#pragma omp parallel for
  for(int runon=0;runon<nummcmcruns;runon++){
    int thread_number = get_omp_get_thread_num();

    string tsv=s4commands[runon];
    if(thread_number>0) tsv.append(" -S");
    std::vector<char *> argv=converttoargv(tsv);
    //    string logfile=fstreeoutputvec[runon];
    //    logfile.append(".log");
#pragma omp atomic
    cmdon++;

#pragma omp critical
    cout<<(ParallelStream()<<"Running stage 4 (tree) command number (~"<<cmdon<<") of "<<nummcmcruns).toString()<<endl;
    if(verbose) cout<<"RUNNING S4 CMD:"<<tsv<<endl;

    //switchStdout(logfile.c_str());

    int rv=finestructure(argv.size(),argv.data());
    freeargv(argv);
    //    revertStdout();
    // Check that it ran correctly
    if(rv!=0){
      cerr<<"FineStructure TREE failed!"<<endl;// See log file ("<<logfile<<") for details."<<endl;
      throw(runtime_error("finestructure"));
      }
}// end for

}

void FsProject::doStagePopExtract(){
  cout<<"ERROR: Unimplemented stage5"<<endl;
}

void FsProject::optionsHelp(){
  cout<<fsoptionshelp;
  cout<<"IMPORTANT PARAMETERS:"<<endl;
  cout<<cmdInfo("idfile",false)<<endl;
  cout<<cmdInfo("phasefiles",false)<<endl;
  cout<<cmdInfo("recombfiles",false)<<endl;
  cout<<"IMPORTANT ACTIONS:"<<endl;
  for(int c1=0;c1<3;c1++){
    cout<<"   ";
    cout<<cmdInfo(cmds[c1].getName(),false)<<endl;
  }
}

void FsProject::outputHelp(){
  cout<<"FILES CREATED, in order of importance."<<endl<<"IMPORTANT FILES:"<<endl;
  cout<<"\t<projectname>.cp: The finestructure parameter file, containing the state of the pipeline."<<endl;
  cout<<"\t<projectname>_<linked>.chunkcounts.out: Created by stage 2 combine: The final chromopainter painting matrix, giving the number of chunks donated to individuals in rows from individuals in columns. The first line containt the estimate of \"c\"."<<endl;
  cout<<"\t<projectname>_<linked>.chunklengths.out: Created by stage 2 combine: The final chromopainter painting matrix, giving the total recombination map distance donated to individuals in rows from individuals in columns."<<endl;
  cout<<"\t<projectname>_<linked>.mcmc.xml: Created by stage 3 combine: The main MCMC file of the clustering performed by fineSTRUCTURE."<<endl;
  cout<<"\t<projectname>_<linked>.tree.xml: Created by stage 4 combine: The main \"tree\" created from the best MCMC state by fineSTRUCTURE."<<endl;  
  cout<<"\t<projectname>: A folder containing all pipeline files."<<endl;
  cout<<"\t<projectname>/commandfiles/commandfile<X>.txt: The commands to be run to complete stage X. (-hpc 1 mode only)"<<endl;
  cout<<"USEFUL FILES:"<<endl;
  cout<<"\t<projectname>/stage<X>: folders containing all pipeline files for a stage X."<<endl;
  cout<<"\t<projectname>/cpbackup/<projectname>.cp<X>.bak: Backups of the parameter file, created after every action."<<endl;
  cout<<"\t<projectname>/stage1/*_EM_linked_file<f>_ind<i>.EMprobs.out: Created by stage 1: The chromopainter parameter estimate files (indexed f=1..<num_phase_files>, in the order given) for the individuals in the order encountered in the idfile (omitting individuals specified as such)."<<endl;
  cout<<"\t<projectname>/stage2/*_mainrun_file<f>_ind<i>.*: Created by stage 2: All chromopainter files created with the same parameters for all individuals (indexed f=1..<num_phase_files>, in the order given) for the individuals in the order encountered in the idfile (omitting individuals specified as such).  See the \"fs cp\" help for details."<<endl;
  cout<<"\t<projectname>/stage3/*_linked_mcmc_run<r>.xml: Created by stage 3: all further MCMC runs beyond the first (r=1..nummcmcruns-1)."<<endl;
  cout<<"\t<projectname>/stage3/*_mcmc.mcmctraces.tab: Created by stage 3 combine: The mcmc samples from all runs in a single file."<<endl;
  cout<<"\t<projectname>/stage4/*_linked_mcmc_run<r>.xml: Created by stage 4: all further trees beyond the first (r=1..nummcmcruns-1)."<<endl;
  cout<<"OTHER FILES:"<<endl;
  cout<<"\t<projectname>_<linked>.mutationprobs.out: Created by stage 2 combine: The final chromopainter painting matrix, giving the *expected number of SNPs donated with error* to individuals in rows from individuals in columns."<<endl;
  cout<<"\t<projectname>_<linked>.regionchunkcounts.out: Created by stage 2 combine: an intermediate file for calculating \"c\". See fs cp help for details."<<endl;
  cout<<"\t<projectname>_<linked>.chunklengths.out: Created by stage 2 combine: an intermediate file for calculating \"c\". See fs cp help for details."<<endl;
  cout<<"\t<projectname>/stage3/*_linked_mcmc_run<r>_x<x>_y<y>_z<z>.xml: Created by stage 3 when MCMC fails convergence tests. This is a backup of where each MCMC run reached, and is used as a starting point for the next run."<<endl;
  cout<<"\t<projectname>/stage<X>/*.log: Log files created by each stage, 1,2,2a (combining stage2 output across chromosomes),3 and 4."<<endl;
}

void FsProject::stagesHelp(){
  cout<<stageshelpheader;
  cout<<"==== pre-stage0 ===="<<endl;
  cout<<"#### stage0 ####"<<endl;
  cout<<"Data conversion. Currently not implemented!"<<endl;
  cout<<"==== post-stage0 ===="<<endl;
  cout<<cmdInfo("countdata")<<endl;
  cout<<"==== pre-stage1 ===="<<endl;
  cout<<"Important note: stage1 is skipped when running in unlinked mode (no recombination file provided)"<<endl;
  cout<<cmdInfo("makes1")<<endl;
  cout<<"#### stage1 #### Chromopainter parameter inference"<<endl;
  cout<<cmdInfo("dos1")<<endl;
  cout<<cmdInfo("writes1")<<endl;
  cout<<"==== post-stage1 ===="<<endl;
  cout<<cmdInfo("combines1")<<endl;
  cout<<"==== pre-stage2 ===="<<endl;
  cout<<cmdInfo("makes2")<<endl;
  cout<<"#### stage2 #### Chromopainter painting"<<endl;
  cout<<cmdInfo("dos2")<<endl;
  cout<<cmdInfo("writes2")<<endl;
  cout<<"==== post-stage2 ===="<<endl;
  cout<<cmdInfo("combines2")<<endl;
  cout<<"==== pre-stage3 ===="<<endl;
  cout<<cmdInfo("makes3")<<endl;
  cout<<"#### stage3 #### FineSTRUCTURE MCMC inference"<<endl;
  cout<<cmdInfo("dos3")<<endl;
  cout<<cmdInfo("writes3")<<endl;
  cout<<"==== post-stage3 ===="<<endl;
  cout<<cmdInfo("combines3")<<endl;
  cout<<"==== pre-stage4 ===="<<endl;
  cout<<cmdInfo("makes4")<<endl;
  cout<<"#### stage4 #### FineSTRUCTURE tree inference"<<endl;
  cout<<cmdInfo("dos4")<<endl;
  cout<<cmdInfo("writes4")<<endl;
  cout<<"==== post-stage4 ===="<<endl;
  cout<<cmdInfo("combines4")<<endl;
  cout<<"Not a command, but if -go gets here, we will provide the GUI command line for visualising and exploring the results."<<endl;
}

void FsProject::getHelp(std::vector<std::string> args){
  unsigned int offset=1;
  while(args.size()>offset){
    /// Handle special cases:
    if(args[offset].compare("all")==0){ // all : do actions and parameters
      args.erase(args.begin()+offset);
      args.push_back("actions");
      args.push_back("parameters");
    }
    if(args[offset].compare("actions")==0 ||args[offset].compare("commands")==0){ // actions: get all actions
      args.erase(args.begin()+offset);
      for(unsigned int c1=0;c1<cmds.size();c1++){
	args.push_back(cmds[c1].getName());
      }
    }
    if(args[offset].compare("parameters")==0){ // parameters: get all parameters
      args.erase(args.begin()+offset);
      for(unsigned int c1=0;c1<pars.size();c1++){
	args.push_back(pars[c1].getName());
      }
    }
    // Separate helps for non-commands/parameters
    if(args[offset].compare("input")==0){ // input format help
      args.erase(args.begin()+offset);
      cout<<inputhelp0<<inputhelpidfile<<inputidfileexample<<endl;
      cout<<inputhelpphase<<inputphaseexample<<endl;
      cout<<inputhelprec<<inputrecexample<<endl<<inputhelp1;
    }else if(args[offset].compare("stages")==0){ // help on what happens in each stage
      args.erase(args.begin()+offset);
      stagesHelp();
    }else if(args[offset].compare("info")==0){
      args.erase(args.begin()+offset);
      cout<<fsprojecthelp<<endl;
    }else if(args[offset].compare("tools")==0){
      args.erase(args.begin()+offset);
      cout<<fstoolshelp<<endl;
    }else if(args[offset].compare("example")==0){ // Run through the example
      args.erase(args.begin()+offset);
      makeExample();
    }else if(args[offset].compare("output")==0){ // help on the files created
      args.erase(args.begin()+offset);
      outputHelp();
    }else{
      // Now we just give help on everything remaining
      cout<<"Help for "<<cmdInfo(args[offset])<<endl;
      args.erase(args.begin()+offset);
    }
  }
}

void FsProject::ensureGtRoot(){
  if(gtroot.compare("")==0){
    stringstream ss0;
    ss0<<dirname<<"/gt/"<<"gt";
    gtroot=ss0.str();
  }
}

void FsProject::ensureDonors(){
  if(gtdonors.size()==0 || gtrecips.size()==0){
    FsDonor donormodel(popidfile,donoridfile,dirname,0,(int)phasefiles.size(),
		       fsmode,verbose);
    if(donoridfile.compare("")==0){
      donoridfile=donormodel.getDonorFile();
    }
    donormodel.readDonorFile();

    gtdonors=donormodel.donorLabels();
    gtrecips=donormodel.recipLabels();
  }
  if(verbose) {
    cout<<"GT working with the following donors: ";
    for(unsigned int c1=0;c1<gtdonors.size();c1++){
      if(c1!=0) cout<<", ";
      cout<<gtdonors[c1];
    }
    cout<<endl<<" and the following recipients: ";
    for(unsigned int c1=0;c1<gtrecips.size();c1++){
      if(c1!=0) cout<<", ";
      cout<<gtrecips[c1];
    }
    cout<<endl;
  }
}

void FsProject::zcatSamples(vector<string> srcfiles, string destfile){

  ///////////////////////////////////
  char buffer[32768];
  if(verbose) cout<<"Creating compressed samples file "<<destfile<<" in directory "<<dirName(destfile)<<endl;
 
  ensureDirectory(dirName(destfile));
  gzFile destgz = gzopen(destfile.c_str(), "w");
  ///////////////////////////////////
  for(unsigned int c1=0;c1<srcfiles.size();++c1){
    srcfiles[c1].append(".samples.out.gz");
    if(verbose) cout<<"Reading compressed samples file "<<srcfiles[c1]<<endl;
    gzFile infile = gzopen(srcfiles[c1].c_str(), "rb");
    int num_read = 0;
    int found_newline=0;
    while ((num_read = gzread(infile, buffer, sizeof(buffer))) > 0) {
      if(found_newline==0 && c1!=0) { // remove the initial line (except on the first file)
	for(int c2=0;c2<num_read;++c2){// Hack-y way to do this; we only need this functionality here
	  if(buffer[c2]=='\n'){
	    for(int c3=c2+1;c3<num_read;++c3){
	      buffer[c3-c2-1] = buffer[c3];
	    }
	    num_read-=(c2+1);
	    found_newline=1;
	    break;
	  }
	}
      }
      gzwrite(destgz, buffer, num_read);
    }
    gzclose(infile);
  }
  gzclose(destgz);
}

void FsProject::ensureCombinedSamples(){
  ensureGtRoot();
  gtsamplesfiles.clear();
  for(unsigned int c1=0;c1<recombfiles.size();++c1){
    stringstream ss0;
    ss0<<gtroot<<"_samples_file"<<c1+1<<".samples.out.gz";
    gtsamplesfiles.push_back(ss0.str());
  }
  if(verbose){
    cout<<"Checking "<<recombfiles.size()<<" samples files with root "<<gtroot<<endl;
  }
  
  unsigned int noutputfiles=s7outputrootvec.size();
  unsigned int nrecfiles=recombfiles.size();
  //  if((((double)noutputfiles) % ((double)nrecfiles)) !=0) {
  if((noutputfiles % nrecfiles) !=0) {
    throw(runtime_error("ensureCombinedSamples: incompatability between samples and s7 output"));
  }
  // Assumption: Each batch of N files is associated with the same recombination file.
  int nlast=0,nthis=noutputfiles/nrecfiles;
  for(unsigned int rf=0;rf<recombfiles.size();++rf){
    // If the output file is empty or doesn't exist
    if( access( gtsamplesfiles[rf].c_str(), F_OK ) == -1 ) { // file does not exist
      vector<string> troots=vector<string>(s7outputrootvec.begin() + nlast,s7outputrootvec.begin() + nthis);
      zcatSamples(troots,gtsamplesfiles[rf]);
      nlast=nthis;
      nthis+=noutputfiles/nrecfiles;
    }else{
      cerr<<"WARNING: Combined samples file "<<gtsamplesfiles[rf]<<" already exists. Not overwriting it. Delete it manually if you want it to be recreated from the stage7 output."<<endl;
    }
  }
}

void FsProject::cleanupSamples(){
  // Removes all the run-specific samples.out.gz files
  for(unsigned int c1=0;c1<recombfiles.size();++c1){
    stringstream ss0;
    ss0<<gtroot<<"_samples_file"<<c1+1<<".samples.out.gz";
    unlink(ss0.str().c_str());
  }
} 

void FsProject::makeDataListFiles(){
  // Where do we make these lists?
  if(gtsamplesfiles.size()!=recombfiles.size()){
    throw(logic_error("GT samples recombfiles mismatch"));
  }
  if(gtsampleslist.compare("")==0){
    stringstream ss0;
    ss0<<gtroot<<"_sampleslist.txt";
    gtsampleslist=ss0.str();
  }
  if(gtrecomblist.compare("")==0){
    stringstream ss0;
    ss0<<gtroot<<"_recomblist.txt";
    gtrecomblist=ss0.str();
  }
  backupFileAndRemove(gtsampleslist);
  backupFileAndRemove(gtrecomblist);

  if(verbose) {
    cout<<"Creating samples list file "<<gtsampleslist<<" and recomb list file "<<gtrecomblist<<endl;
  }
  writeStringVectorToFile(gtsamplesfiles,gtsampleslist,
			  std::vector<std::string>(),false);
  writeStringVectorToFile(recombfiles,gtrecomblist,
			  std::vector<std::string>(),false);
 
}

void FsProject::makeGt(){
  // First, figure out the donor labels
  ensureDonors();
  // Second, create samples files for each chromosome
  ensureCombinedSamples();
  // Third, make the list of recombination files and samples files
  makeDataListFiles();
  // Then, make the K parameter files: one for each recipient population
  s9commands.clear();
  gtparamfiles.clear();
  s9logfiles.clear();
    for(unsigned int c1=0;c1<gtrecips.size();++c1){
    writeGtParams(gtrecips[c1]);
    stringstream ss;
    ss<<"R < GLOBETROTTER.R "<<gtparamfiles[c1]<<" "<<gtsampleslist<<" "<<gtrecomblist<<" --no-save";
    s9commands.push_back(ss.str());
  }
  //  R < GLOBETROTTER.R [parameter infile] [painting samples filelist infile] [recom rate filelist infile] --no-save
}

void FsProject::writeGtParams(string recip){
  string copyvectorsfile=popcproot;
  copyvectorsfile.append(".chunklengths.out");
  stringstream ss0;
  ss0<<gtroot<<"_"<<recip;
  string gtmain=ss0.str(),gtboot=ss0.str(),gtparam=ss0.str(),gtlog=ss0.str();
  gtmain.append(".main");
  gtboot.append(".boot");
  gtparam.append(".param");
  gtlog.append(".log");
  if(verbose) cout<<"Creating GT parameter file "<<gtparam<<endl;
  gtparamfiles.push_back(gtparam);
  s9logfiles.push_back(gtlog);
  
  filebuf fb;
  try{
    fb.open (gtparam.c_str(),ios::out);
  }catch(exception &x){
    cerr<<"Error opening file!"<<endl<<x.what()<<endl; 
    throw(runtime_error("-makes9: cannot write to file!"));
  }
  ostream os (&fb);
  os<<"prop.ind: 1"<<endl;
  os<<"bootstrap.date.ind: 0"<<endl;
  os<<"null.ind: 0"<<endl;
  os<<"input.file.ids: "<<popidfile<<endl;
  os<<"input.file.copyvectors: "<<popcproot<<".chunklengths.out"<<endl;//example/BrahuiYorubaSimulation.copyvectors.txt
  os<<"save.file.main: "<<gtmain<<endl;// example/BrahuiYorubaSimulation.globetrotter.main
  os<<"save.file.bootstraps: "<<gtboot<<endl; //example/BrahuiYorubaSimulation.globetrotter.boot
  os<<"copyvector.popnames:";
  for(unsigned int c1=0;c1<gtdonors.size();++c1) os<<" "<<gtdonors[c1];
  os<<endl; //Adygei Armenian ...
  os<<"surrogate.popnames:";
  for(unsigned int c1=0;c1<gtrecips.size();++c1){
    if(gtrecips[c1].compare(recip)!=0) os<<" "<<gtrecips[c1];
  }
  os<<endl; //Adygei Armenian ...
  os<<"target.popname: "<<recip<<endl;//BrahuiYorubaSimulation
  os<<"num.mixing.iterations: "<<gtnummixingiterations<<endl;
  os<<"props.cutoff: "<<gtpropscutoff<<endl;
  os<<"bootstrap.num: "<<gtbootstrapnum<<endl;
  os<<"num.admixdates.bootstrap: 1"<<endl;
  os<<"num.surrogatepops.perplot: 3"<<endl;
  os<<"curve.range: 1 50"<<endl;
  os<<"bin.width: "<<gtbinwidth<<endl;
  os<<"xlim.plot: 0 50"<<endl;
  os<<"prop.continue.ind: 0"<<endl;
  os<<"haploid.ind: "<<(ploidy==1)<<endl;

  fb.close();
  
}

void FsProject::makeExample(){
  try{
    writeStringToFile(inputidfileexample,"exampledata.idfile");
    writeStringToFile(inputphaseexample,"exampledata.phase");
    writeStringToFile(inputrecexample,"exampledata.recombfile");
    cout<<"Created three example input files:\n\
exampledata.idfile : the id file\n\
exampledata.phase : the phase file\n\
exampledata.recombfile : the recombination file.\n\
To run this example, try:\n";
    cout<<exampletext;
    cout<<"Note: because this example dataset has so few snps, we have to specify the \"chunks per region\" (via -s2chunksperregion) for it to work. Also try omitting this, you should be prompted with how to fix it!"<<endl;
    cout<<"Also note that if you don\'t change the file name, e.g. fs example2.cp, fs will try to continue your previous run, which may have completed. Restart it anew with -n."<<endl;
  }catch(runtime_error& e){
    cerr<<"Error creating example files!"<<endl<<e.what()<<endl; 
    throw(runtime_error("fsproject: cannot create examples"));
  }
}

////////////////////////////////
//
int FsProject::applyFiles(string readfile,string actionifdirexists){
    int dval=directoryExists(dirname);
    if(dval<0){// is a file, we are stuck
	cerr<<"ERROR: project directory with name "<<dirname<<" exists but is not a directory. Delete or move this manually and rerun!"<<endl;
	return(-1);
    }
    if(dval==1){ // is a directory, probably from a previous run
      if(actionifdirexists.compare("stop")==0){
	cerr<<"ERROR: project directory with name "<<dirname<<" already exists. Delete or move this manually and rerun!"<<endl;
	return(-1);
      }else if(actionifdirexists.compare("delete")==0){
	cerr<<"IMPORTANT: You have specified to remove the previously existing directory "<<dirname<<". This cannot be undone, and all backup information will be lost. You have until 0 to press Ctrl-C and cancel this command."<<endl;
	for(int i=3;i>=0;--i){
	  cout<<i<<"..."<<endl;
	  sleep(1);
	}
	cout<<"Deleting!"<<endl;
	deleteFolderTree(dirname);
	dval=directoryExists(dirname);
	deleteFolderTree(filename);
      }else if (readfile.compare("new")==0){
	cerr<<"WARNING: Removing old command files."<<endl;
	stringstream ss;
	ss<<dirname<<"/commandfiles";
	deleteFolderTree(ss.str());
      }
    }
    
    int fval=access( filename.c_str(), F_OK );
    if(fval != -1) { // file exists
      if((readfile.compare("detect")==0) | (readfile.compare("read")==0)) { // read it
	readFromFile();
	if(verbose) cout<<"Reading project file "<<filename<<"..."<<endl;
      }else{ // exists but told not to read it!
	cerr<<"WARNING: file "<<filename<<" exists but is being overwritten! You may wish to rerun without the \"-n\" option. If this command completes, you will also need to restore from the backup file. Ignore this warning if you meant to overwrite the project."<<endl;
      }
      if(dval!=1) { // directory existx
	if(verbose) cout<<"Project file exists but directory does not. Creating project directory "<<dirname<<"..."<<endl;
	ensureDirectory(dirname);
      }
    }else{ // no file 
      if(readfile.compare("read")==0){  // doesn't exist but told to read it!
	cout << "Told to read a file (-n option) but this is not possible."<<endl;
	string fexists="does not exist", dexists="does not exist";
	if(fval != -1) fexists="exists";
	if (dval==1) dexists="exists";
	cout<<" filename " << filename<<" "<<fexists<<endl;
	cout<<" directory " << dirname<<" "<<dexists<<endl;
	return(1);
      }
      if(verbose) cout<<"Creating new project directory "<<dirname<<"..."<<endl;
      ensureDirectory(dirname);
      if(verbose) cout<<"Creating new project file "<<filename<<"..."<<endl;
      try{
	filebuf fb;
	fb.open (filename.c_str(),ios::out);
	if (! fb.is_open() )
	  throw(runtime_error("fsproject: cannot write to file!"));
      }catch(exception &x){
	cerr<<"Error opening file!"<<endl<<x.what()<<endl; 
	throw(runtime_error("fsproject: cannot write to file!"));
      }
    }
    return(0);
}

int fsproject(int fsmode,int argc, char *argv[]) {
    unsigned int argon=0, argtmp=0;
    std::string filename=string("fsproject.cp");
    std::string dirname=string("./fsproject");
    bool verbose=false;
    int endstatus=0;
    std::string readfile=string("detect");
    std::string diraction=string("merge");
    std::vector<std::string> args(argv, argv+argc) ;
    args.erase(args.begin());

    FsProject *proj;
    try{ proj = new FsProject(filename,dirname,verbose);
    }catch (exception& e)  {
      cout << "Standard exception: " << e.what() << endl;
      return(1);
    }
    
    // Check for general options
    // Check if there are no options
    if(args.size()==0){
      proj->optionsHelp();return 0;
    }
    
    // check for, and read the filename 
    if(args[0].substr(0,1).compare("-")!=0){ // no leading \"-\"
      filename=string(args[0]);
      dirname=projectroot(filename);
      filename=projectfull(dirname);
      if(args[0].compare(filename)!=0){
	if(fsmode>0) cout<<"WARNING: you have not included the \".cp\" ending in the project name. This can cause problems if you try to call it the same name as a fs tool.\n";
      }
      if(verbose) cout<<"Using filename "<<filename<<" with directory "<<dirname<<endl;
      //interpret as a filename
      args.erase(args.begin());
    }else if((args[0].compare("-help")!=0 && args[0].compare("-h")!=0)){
      cout<<"Must provide the filename <projectname>.cp before any actions or parameter settting. See \"fs -h\" for help on this mode, or \"fs\" for general help."<<endl;
      return 0;
    }

    // Process the remaining parameters
    argon=0;
    while(argon<args.size()) {
      argtmp=argon;
      if(args[argon].compare("-h")==0 || args[0].compare("-help")==0){ // help!
	std::vector<std::string> args1(args.begin()+argon, args.begin()+args.size()) ; // the arguments of the command
	if(args1.size()==1) proj->optionsHelp();
	else proj->getHelp(args1); 
	return 0;
      }else if(args[argon].compare("-v")==0) {	// verbose
	verbose=true;
	cout<<"Verbose mode"<<endl;
	args.erase(args.begin() + argon);
      }else if(args[argon].compare("-n")==0) {	// new file
	readfile=string("new");
	diraction=string("merge");
	args.erase(args.begin() + argon);
      }else if(args[argon].compare("-N")==0) {	// new file
	readfile=string("new");
	diraction=string("delete");
	args.erase(args.begin() + argon);
      }else{
	argon++;
      }
    }


    ///////////////////////////
    // set up project structure
    if(verbose) cout<<"Using directory "<<dirname<<" and file "<<filename<<endl;
    proj->setFileName(filename);
    proj->setDirectoryName(dirname);
    proj->setVerbose(verbose);
    proj->setFsmode(fsmode);

    ///////////////////////////
    // read file if appropriate
    int fval;
    try{
      fval=proj->applyFiles(readfile,diraction);
    }catch(exception &x){
      cerr<<"fsproject: cannot write to file: "<<filename<<". Check permissions and file size."<<endl;
      exit(-1);
    }
    if(fval==1) return(fval);
    proj->addHistory(args);
    ////////////////////////////

    if(verbose) cout<<"Processing commands..."<<endl;
    argon=0;
    argtmp=0;
    int cmdon=0;
    int lastcmdon=cmdon;
    while(argon<args.size()){
      if(verbose) cout<<"Processing command argument "<<cmdon<<":";
      if(verbose) cout<<" \""<<args[argon];
      argtmp=argon+1;
      while(argtmp<args.size() && args[argtmp].at(0)!='-') {
	if(verbose) cout<<"\" \""<<args[argtmp];
	argtmp++;
      }
      if(verbose) cout<<"\""<<endl;
      std::vector<std::string> args1(args.begin()+argon, args.begin()+argtmp) ; // the arguments of the command
      args1 = getcommands(args1);/// extract the arguments split by comma
      try{
	proj->docmd(args1); // do the command
	argon=argtmp;
	cmdon++;
      }catch (runtime_error& e)  {
	string swhat=e.what();
	char * swhat_c=(char*) swhat.data();
	char * tmissing, * texists, * tincomplete;
	tmissing = strstr (swhat_c,"missing");
	texists = strstr (swhat_c,"exists");
	tincomplete = strstr (swhat_c,"incomplete");
	endstatus=-1;
	if ((proj->getHpc()==1)&&(tmissing!=NULL || texists!=NULL|| tincomplete!=NULL)) {
	  int stage=0;
	  string tcmdfile;
	  if(strstr (swhat_c,"stage1")!=NULL){
	    stage=1;
	  }else if(strstr (swhat_c,"stage2")!=NULL){
	    stage=2;
	  }else if(strstr (swhat_c,"stage3")!=NULL){
	    stage=3;
 	  }else if(strstr (swhat_c,"stage4")!=NULL){
	    stage=4;
 	  }else if(strstr (swhat_c,"stage6")!=NULL){
	    stage=6;
 	  }else if(strstr (swhat_c,"stage7")!=NULL){
	    stage=7;
 	  }else if(strstr (swhat_c,"stage9")!=NULL){
	    stage=9;
	  }else{
	    cerr << "ERROR: chromopainter failed. Check the input files and logs ... saving progress so far." <<endl;
	    argon=args.size();
	    break;
	  }
	  // If we get here, it is a valid hpc escape
	  cmdon++;
	  argon=argtmp;
	  endstatus=1;
	  tcmdfile=proj->getCommandFile(stage);
	  if(tmissing!=NULL) {
	    cout<<"HPC mode: "<<proj->getCommandFileLength()<<" commands for stage"<<stage<<" are in file "<<tcmdfile<<". Rerun when those commands have been completed and the results copied back to this directory."<<endl;
	  }else if(texists!=NULL) {
	    cout<<"HPC mode: "<<tcmdfile<<" already exists. Not overwriting it. Did it get processed?"<<endl;
	  }

	  if(tincomplete!=NULL) {
	    cout<<"HPC mode: The commands were not completed. You may wish to run -remakes"<<stage<<" to recreate only those commands that need to be rerun."<<endl;
	  }else{
	    cout<<"SUGGESTIONS:"<<endl;
	    cout<<"> cat "<<tcmdfile<<" | parallel --bar # for local execution in parallel"<<endl;
	    cout<<"> sbatch_array.sh -f "<<tcmdfile<<" -n "<<proj->recommendN()<<" -m "<<proj->recommendM()<<" # for SLURM-based HPC systems -n <n> is the number of cores requested on each HPC node and -m <m> is the number of commands sent to each node (not core). These values are suggested to keep the total number of jobs low ~(10-100)."<<endl;
	  }
	  if(texists!=NULL) {
	    cout<<"IF IT WAS RUN ALREADY:"<<endl;
	    cout<<"1. Is it still being executed? If so, wait for it to finish and rerun with -go."<<endl;
	    cout<<"2. Is it completed with some failures? Recreate a list of the remaining tasks with -remakes"<<stage<<"."<<endl;
	  }
      }else if (swhat.compare("chromopainter")==0){
	  cerr << "ERROR: chromopainter failed. Check the input files and logs ... saving progress so far." <<endl;
	}else if (swhat.compare("chromocombine")==0){
	  cerr << "ERROR: chromocombine failed. Check the input files and logs ... saving progress so far." <<endl;
	}else{
	  cerr << "ERROR: "<<e.what()<<" ... saving progress so far." <<endl;
	}
      }catch (std::exception& e)  {
	cerr << "ERROR: Failed to perform command:";
	for(unsigned int i=0;i<args1.size();i++){
	  cerr<<" "<<args1[i];
	}
	cerr<<endl;
	endstatus=-1;
	/// For the moment I'm a little confused as to how we end up here from a runtime_error.
	cerr << "ERROR DETAILS: General Exception: " << e.what() << endl <<"... saving progress so far."<<endl;
	break;
      }catch(...) {
	cerr<<"ERROR! Unknown error."<<endl;
      }
      if(lastcmdon==cmdon){
	cerr << "ERROR: Failed to perform command:";
	for(unsigned int i=0;i<args1.size();i++){
	  cerr<<" "<<args1[i];
	}
	cerr<<endl;
	argon=args.size();
      }else lastcmdon=cmdon;
    }

    //////////////////////////////
    /*
    if(verbose) cout<<"Configuring project file..."<<endl;
    try{ 
      proj->safeCreateFile();
      if(verbose) cout<<"Writing project file..."<<endl;
      proj->writeToFile();
    }catch (exception& e)  {
      cout << "Standard exception: " << e.what() << endl;
    }
    */
    if(endstatus<0 && fsmode>0 ) {
      cout<<"IMPORTANT: The run ended on an error. You should read the error above and correct it. Progress was saved BEFORE the last command was attempted."<<endl;
    }else if(endstatus>0){
      string tcmdfile=proj->getCommandFile();
      cout<<"IMPORTANT: The run ended with a requirement to run commands externally from file \""<<tcmdfile<<"\". Once you have done this, resume the analysis with \""<<proj->getExec()<<" "<<filename<<" -go\""<<endl;      
    }
    delete(proj);
    return(0);
}
