#ifndef FSCONSTANTS_H
#define FSCONSTANTS_H

    
static const char * fsprojecthelp=
  "\
***** FineSTRUCTURE and ChromoPainter *****\n\
USAGE: \"fs <projectname>.cp <options> <actions>\" \n\
SUMMARY:\n\
    * Creates <projectname>.cp and a directory <projectname>. \n\
    * Performs inference using ChromoPainter and FineSTRUCTURE. \n\
    * Can create commands to be run on an external HPC machine. \n\
    * Supports parallel processing in a single machine, obtain finestructure results from a single command.\n\
    * This is helper software to easily run chromopainter and finestructure, which are built in to this program.\n\
OVERVIEW OF USE: \n\
    * Provide SNP data with -phasefiles, recombination data with -recombfiles and individual data with -idfile\n\
    * Configure any advanced parameters\n\
    * And -go: fs will figure out the rest!\n\
    * If using -hpc mode, run \"fs\" on the head node, which will not run the model. Instead four \"stages\" of processing are created for running on a external HPC system. Then rerun this program once they have completed to generate the next stage.\n\
USE \"fs -h\" for help on the automatic mode.\n\
";


static const char * fsoptionshelp="\
***** Help for fs - running the whole chromopainter/finestructure inference pipeline in 'automatic' mode *****\n\
USAGE: \"fs <projectname>.cp <options> <actions>\" \n\
GENERAL OPTIONS FOR \"project\" tool: \n\
    -h/-help:    Show this help.\n\
    -help info: Show 'overview' help explaining how this software works.\n\
    -help actions: Show help for all actions.\n\
    -help parameters: Show help for all parameters.\n\
    -help <list of commands or parameters>: Show help on any specific commands or parameters.\n\
    -help input: Show examples and give details of the input file formats.\n\
    -help output: Details of the files that may be created.\n\
    -help stages: Detailed description of what happens in, and between, each stage of the computation.\n\
    -help tools: Show help on how to access the chromopainter/chromocombine/finestructure tools directly.\n\
    -help example: Create and show help for a simple example.\n\
    <tool> -h: Show help on a particular tool: one of fs,cp,combine. IMPORTANT NOTE: These have simplified interfaces with different names when running in automatic mode. The help for automatic mode parameters explains which parameters it changes.\n\
    -v      :    Verbose mode\n\
    -n      :    New settings file, overwriting any previous file\n\
    -N      :    Remove the project directory entirely\n\
    -<parameter>:<value> : Sets the internal parameter, exactly as if they were read in from -import. \n\
The colon is optional, unless <value> starts with a '-' symbol. Escape spaces and don't use quotes; \n\
e.g. '-s1args:-in\\ -iM'.\n\
    \n\
";


static const char * stageshelpheader=
    "\
***** Help on computational stages *****\n\
The computation for finestructure is split into 4 main stages.\n\
These are breakpoints at which we can export computation to a HPC machine. Before and after each, automatic mode will do the work necessary to construct the next stage. This includes the construction of the command lines to be executed; the command lines themselves are all that is run externally.\n\
pre-stage<x>: performed by fs. A -reset <x> command will result in this being redone.\n\
post-stage: performed by fs A -reset <x> command will use the output of the post-stage<x-1> processing. This means, for example, that we can avoid needing to duplicate the chromopainter (stage2) runs in order create a duplicated finestructure (stage3) run.\n\
stage: Either '-dop<x>' meaning that the previously generated commands are run internally (in parallel) or '-writep'<x>' meaning that they are written to file to be performed externally in HPC mode.\n\
DETAILS:\n\n\
";

static const char * fstoolshelp=
    "\
***** Help for tool mode *****\n\
USAGE: fs [tool] [OPTIONS] \n\
Using this interface you can access any of the advanced functionality of chromopainter and finestructure.\n\
    [tool] can be any of: \n\
    <projectname>.cp: automatic mode - creates and runs the commands below for you, organising what to do in a \n\
'project'. This should be the default unless you know what you are doing.\n\
    cp: chromopainter mode (commands exactly as chromopainter.) This can be used to perform more sophisticated analyses, such as admixture modelling via GLOBETROTTER.\n\
    combine: chromocombine mode (commands exactly as chromocombine) \n\
    fs: finestucture mode (commands exactly as finestructure) \n\
    \n\
    Run \"fs [tool] -h\" to obtain more detailed help on cp, combine, or fs tools.\n\
    Run \"fs -h\" to get the automatic mode help.\n\
    ";

static const char * inputhelp0=
"\
##########################\n\
HELP ON INPUT FORMATS.  \n\n\
This help, combined with looking at the example and the use of the provided scripts to convert your \n\
data, should be enough for most users to get started. \n\
NOTE: You can specify multiple phase and recombination files, one for each chromosome (at least, they are assumed unlinked.) Specify via:\n\
-phasefiles <list>.phase <of>.phase <files>.phase\n\
with corresponding:\n\
-recombfiles <list>.rec <of>.rec <files>.rec\n\n";

static const char * inputhelpidfile="\
##########################\n\
IDFILE FORMAT:\n\
This specifies the names of the individuals in the data, as well as (optionally) which population they are from and whether they are included.\n\
Format: N lines, one per individual, containing the following columns:\n\
<NAME> <POPULATION> <INCLUSION> <ignored extra info>\n\
    Where <NAME> and <POPULATION> are strings and <INCLUSION> is 1 to include an individual and 0 to exclude them. The second and third columns can be omitted (but the second must be present if the third is). Currently <POPULATION> is not used by this version of fs.\n\
EXAMPLE IDFILE:\n";
static const char * inputidfileexample="\
Ind1 Pop1 1\n\
Ind2 Pop1 1\n\
Ind3 Pop2 0\n\
Ind4 Pop2 1\n\
Ind5 Pop2 1\n\
";
static const char * inputhelpphase="\
##########################\n\
CHROMOPAINTER'S v2 'PHASE' FORMAT:\n\
This is heavily based on 'FastPhase' output. \n\
* The first line contains the number of *haplotypes* (i.e. for diploids, 2* the number of individuals).\n\
* The second line contains the number of SNPs.\n\
* The third line contains the letter P, followed by the basepair location of each SNP (space separated). These must match the recombination file. Within each chromosome, basepairs must be in order.\n\
* Each additional line contains a haplotype, in the order specified in the IDFILE. Diploids have two contiguous rows. Each character (allowing no spaces!) represents a *biallelic* SNP. Accepted characters are 0,1,A,C,G,T, with NO missing values!\n\
EXAMPLE PHASEFILE:\n";
static const char * inputphaseexample="\
10\n\
6\n\
P 100 200 300 400 500 600\n\
010101\n\
011101\n\
111101\n\
001101\n\
011000\n\
001100\n\
001001\n\
001011\n\
001001\n\
001111\n\
";
static const char * inputhelprec="\
##########################\n\
CHROMOPAINTERS RECOMBINATION FILE FORMAT:\n\
Required only if running in unlinked mode.\n\
This specifies the distance between SNPs in 'recombination rate' units. There should be a header line followed by one line for each SNP in haplotype infile. Each line should contain two columns, with the first column denoting the basepair position values given in haplotype infile, in the same order. The second column should give the genetic distance per basepair between the SNP at the position in the first column of the same row and the SNP at the position in the first column of the subsequent row. The last row should have a '0' in the second column (though this is not required – this value is simply ignored by the program). Genetic distance should be given in Morgans, or at least the relevant output files assume this value is in Morgans.\n\
If you are including genetic information from multiple chromosomes, put a '-9' (or any value < 0) next to the last basepair position of the preceeding chromosome.\n\
EXAMPLE RECOMBFILE:\n";
static const char * inputrecexample="\
start.pos recom.rate.perbp\n\
100 0.01\n\
200 0.02\n\
300 -9\n\
400 0.02\n\
500 0.05\n\
600 0\n\
";
static const char * inputhelp1="\
##########################\n\n\
See the dedicated ChromoPainter v1 manual for more details.\n\
";
static const char * exampletext="\
fs example.cp -idfile exampledata.idfile -phasefiles exampledata.phase -recombfiles exampledata.recombfile -s2chunksperregion 1 -go\n\
";


const std::string constnamepre0 = string("comment (information about the project, what has been done etc)\n");
const std::string constname0 = string("############## stage0 (data processing to create the required input files):\n");
const std::string constname1 = string("############## stage1 (chromopainter parameter estimation via EM)\n");
const std::string constname2 = string("############## stage2 (chromopainter full painting)\n");
const std::string constname3 = string("############## stage3 (finestructure mcmc)\n");
const std::string constname4 = string("############## stage4 (finestructure tree)\n");
const std::string constname5 = string("############## stage5 (reference population construction)\n");
const std::string constname6 = string("############## stage6 (chromopainter admixture model parameter estimation via EM)\n");
const std::string constname7 = string("############## stage7 (chromopainter admixture model full painting\n");
const std::string constname8 = string("############## stage8 (admixture model)\n");

const std::string constcommentspre = string("# This file contains the information about the project. It has the following sections:\n");
const std::string constcommentspost = string("#\n\
# You may modify this file by hand, but it is possible to break it if you are careless.  It is advisable instead to set parameters via the \"fs\" interface; either from the command line or by placing the set of parameters to set in a settings file and using -import.\n\
# Each section is separated by \"section:<sectionname>\".  Within each section the parameters are defined in the format:\n\
# paramname:value\n\
# A '#' symbol indicates the end of the readable line, everything afterwards is treated as a comment.\n\
# Everything from the colon to the end of the line (or first #) is treated as content. Parameters that are empty are assigned their default values as described in \"fs -h\".\n\n\
############## Command history:\n");

const std::string constccerror=string("ERROR: 'c' could not be inferred. This is usually because there are not enough chunks to perform cross validation.\n SOLUTIONS:\n\
a) use whole genomes as regions: 'fs <fsfile.cp> -s2combineargs:-C -go'. Try this if your regions have too few chunks (e.g. less than 20).\n\
b) Rerun stage2 with fewer chunks per region: 'fs <fsfile.cp> -reset 2 -s2chunksperregion <value> -go'. This sets chromopainter's chunks per region (cp -k) to <value>. Choose this to be less than the rowsums of each chromosome of each individual; ideally not less than 20.");
const std::string constccerror7=string("ERROR: 'c' could not be inferred. This is usually because there are not enough chunks to perform cross validation.\n SOLUTIONS:\n\
a) use whole genomes as regions: 'fs <fsfile.cp> -s7combineargs:-C -go'. Try this if your regions have too few chunks (e.g. less than 20).\n\
b) Rerun stage2 with fewer chunks per region: 'fs <fsfile.cp> -reset 7 -s7chunksperregion <value> -g'. This sets chromopainter's chunks per region (cp -k) to <value>. Choose this to be less than the rowsums of each chromosome of each individual; ideally not less than 20.");
const  std::string  constdupreset = string("If this was intended then try duplicating the project (-duplicate) or resetting it (-reset).");

#endif

