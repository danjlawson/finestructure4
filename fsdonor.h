#ifndef FSDONOR_H
#define FSDONOR_H
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <fstream>
#include <algorithm>
#include <string.h>
#include <map>

#include "finestructure/safegetline.h" // for cross-platform safely getting lines
#include "fsutils.h"

using namespace std;

namespace fines
{

/**
    @brief Manipulation, input and extraction of donor files
*/
class FsDonor
{
 public:
  FsDonor(std::string idfile,std::string donorfile,std::string cpdir, int runid,int nsamplefiles,int fsmode,bool verbose); /// constructor; calls readIdFileContents. runid needed to allow rerunning with different random subsets of donors. Set donorfile to "" for all pops vs all pops
  ~FsDonor();

  void readIdFileContents();///*read the file contents and create all the useful information for creating process-specific donor files
  void indexPopulations();///* Creates the indices we need to rapidly create populations
  void createDonorFile();///< Create a donor file for the processing
  void createSampleDonorFile();///< Create a donor file usable for external processing
  void readDonorFile();///< Read a donor file for the processing
  void createIdFileForInd(int i);///* Creates the id files and donor files needed for the analysis
  void createRequiredFiles();///* Creates the id files and donor files needed for the analysis
  std::vector<std::string> getCommandContent();///< Gets the flags needed to pass into a run for each individual
  void createDirectory();///< Create the directory for all the files
  std::vector<int> sampleInclusion(int i);///< samples the inclusion to make sure we have n_i-1 inds from every population
  int countPosition(vector<int> myinclusion,string myname);///< Count the Realised index (in terms of included inds) for a named individual of the position they appear in, for their own personal file
  void countRecipients(); ///< Count the number of recipients that are in the data
  inline int numRecipients() {
    return(Nrecip);
  }///< return the number of recipient individuals
  inline string getDonorFile(){
    return(donorfile);
  }
  inline string getReferenceDonorFile(){
    return(sampledonorfile);
  }
  inline vector<string> getReferenceIdFiles(){
    return(refidfiles);
  };///< return the list of the id files created for use by samples
  vector<string> donorLabels();///< return the labels of the donors
  vector<string> recipLabels();///< return the labels of the recipients
 protected:
  // variables passed into the object
  std::string idfile;///< File from which we read ids
  int runid;///< Id number for this donor object, so that we can create multiple random samples of donors
  int Nsamplefiles;///< Number of sample files that we will create, each with a different random individual removed from each population
  bool verbose;///< Verbosity
  int fsmode;///< 1 if running in classic fs mode, 0 if limited
  
  // content of the donor files
  std::map <string, int> dfiledonors; // donor file donor populations
  std::map <string, int> dfilerecips;///< Map from popnames to popmembers index
  std::vector <string> donornames; // donor file donor populations
  std::vector <string> recipnames;///< recip populations
  //  std::vector<std::string> dfilerecips; // donor file recipient populations
  
  // content of the ID files for each individual
  std::vector<std::string> allids; ///* Ids, including those that will be excluded from analysis
  std::vector<std::string> usedids; ///* Ids that we will use in general
  std::vector<int> allinclusion;///* Vector of inclusions

  // content of the ID files populations
  std::vector<std::string> allpops;///* Population labels for all individuals
  std::vector<std::string> usedpops;///* Population labels for only included individuals
  
  // derirved quantities
  int Nraw;///* number of individuals in the raw data file
  int Nused;///* number of individuals used
  int Nrecip;///* Number of recipient individuals in the id file referenced by the donor
  
  std::vector<std::string> popnames;///< Names of each population
  std::map <string, int> popmap;///< Map from popnames to popmembers index
  std::vector<std::vector<int> > popmembers;///< A list of the populations. Each contains a list of the raw indices, for each of the included individuals
  std::vector<int> realisedindex;///< Realised index for each used individual of the position they appear in, for their own personal file
  
  string cpdir;///< Directory where the main data are
  string dirname;///< Directory where we keep the idfiles

  std::vector<std::string> idfiles;///< idfiles that we create
  std::vector<std::string> refidfiles;///< reference idfiles that we create
  std::string donorfile;///< donorfiles that we create
  std::string sampledonorfile;///< donorfile to use for samples
  //  std::vector<int> usedpopindices;///* Indexes of the population ID for all used individuals
  //  std::vector<std::vector<int> > popIndices;///< List of indices of the populations
};

} // end namespace fines
#endif
