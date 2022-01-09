#ifndef CHROMOCOMBINEFUNS_H_
#define CHROMOCOMBINEFUNS_H_

#include <vector>
#include <string>
#include <iomanip>
#include <stdio.h>     /* for printf */
#include <stdlib.h>    /* for exit */
#include <unistd.h>    /* for getopt */

using namespace std;

class ChromoCombineWorker 
{
  public:
  ChromoCombineWorker();///<Constructor
  ~ChromoCombineWorker();///<Destructor
  vector<string> inputfilerootraw;///< What the user said they wanted to run on
  vector<string> inputfileroot;///< The files that this implies
  vector<long> validfile;///< Flag for whether the required set of files exist
  vector<string> possiblemiddles;///< the list of all types of file (e.g. chunkcount)
  vector<string> possibleendings;///< the list of all file endings we know about
  vector<string> requiredendings;///< the list of all file endings to be processed
  

  string outputfileroot;///<Root for the input files
  string outputfileending;///<Ending for the input files
  string outputfileoutending;///<Ending for the output files
  vector<string> rownames;///< row names IN THE COUNTS FILE
  vector<string> colnames;///< col names IN THE COUNTS FILE
  vector<string> maincolnames;///< col names as desired from the id file, or counts file, whichever is best
  string topleftentry;///< Entry that appears in the top left of each file

  vector<vector <double> > counts;///< matrix of counts
  vector<vector <double> > regioncounts;///< matrix of region counts
  vector<vector <double> > regioncountssq;///< matrix of region sums squared
  vector<vector <double> > lengths; ///< matrix of lengths
  vector<vector <double> > muts; ///< matrix of mutations
  vector<long> indexcounts;///< Number of times that we've seen each index

  vector<double> regionsums;///< the regionsums (just the 2nd column on regioncounts & regioncountssq)
  vector<vector <double> > empiricalvar;///< matrix of empirical variances
  vector<vector <double> > theoreticalvar;///< matrix of theoretical variances

  string forcefile;///< The name of the force file
  string forcefileoutput;///< The name of the force file output
  vector<string>forceconts;///< List of names for continents
  vector<vector<string> >forceinds;///< List of individuals within each continent (names)
  vector<vector<long> >forceindsid;///< List of individuals within each continent (indexes)
  vector<long> contin;///<index for the continent each individual is in, -ve (and with abs value equal to their order) if not in
  long numeffinds;///<effective number of individuals
  bool useforce;///< do we use a force file?

  vector<vector<double> > forcedregioncounts;///<regioncounts after application of the force file
  vector<vector<double> > forcedregioncountssq;///<regioncountssq after application of the force file
  vector<vector<double> > forcedregioncountssumsq;///<regioncounts^2 summed after application of the force file
  vector<double> forcedregionsums;///<regionsums after application of the force file

  bool wantcounts;///< do we want counts processed?
  bool wantlengths;///< do we want lengths processed?
  bool wantmuts;///< do we want mutations processed?
  long testrun;///< are we doing a test run?

  int verbose;///<Whether to chat about how things are going
  double cval;///< The value of c that we've calculated
  bool completegenomes;///< true if we should ignore regioncounts and regioncountssquared and instead use each input file as a region
  double completegenomesnormfactor;///<normalising factor in the complete genomes case
  double completegenomesnormfactorsq;///<normalising factor in the complete genomes case
  double completegenomesnormfactorcount;

  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  // FUNCTIONS


  ////////////////////////////////////////////
  // DEFINITIONS AND SETUP
  void makeEndings();///< Makes the endings, given what has been requested
  void clearAllData();///<Clears all data, but not the inputs
  ////////////////////////////////////////////
  // ORGANISING INPUT FILES
  void getFilesFromDirectory(string dir,vector<string> *dirlisting);///< dumpts all files from dir into dirlisting
  bool searchForEndings(string str);///< searches str for all interesting endings
  void getInterestingFiles(vector<string> allfilenames);///< gets all files from the vector that we are looking for
  void getFilesFromDirectories();///<(re)Fills the raw input files assuming they currently contain directories, and removes files that don't have the correct structure
  vector<string> getFileList(string root,vector<string> endings);///< returns the files that we will want, given a root
  void testFilesExist();///< Tests for the existence of the required files

  string getRoot(string in);///<Turns a file into a file root
  bool rationaliseInputFiles();///< forms a file list of those we have complete details of

  ////////////////////////////////////////////
  // FORCE FILES, READING AND APPLYING
  void addSuper(string name,vector<string> ind);///< just adds the superindividual to the list
  void makeSuper(std::string superstring); ///< makes some super individuals as desribed by superstring (which is all continents)
  void makeSuperFromFile(std::string filename); ///< makes some super individuals as desribed by superstring from a file name containing the string
  void getForceIds();///< Converts the names in the force file to ids
  void applyForceFile();///< populates the force* matrices and vectors from the current force ids
  bool writeForceFile(double forceC);///< Writes the force file

  ////////////////////////////////////////////
  // READING DATA
  std::string removequotes(std::string str);///< Gets rid of quotes around a string
  long getRowIndex(string dval);///<Looks through the column names for the current name and returns the index (updating the index counts), allowing for missing names
  long getColIndex(string dval);///<Looks through the column names for the current name and returns the index, insisting on finding the name.
  std::istream& safeGetline(std::istream& is, std::string& t);///< gets a line from a file safely

  bool readIds(string filename);//< Read the column headers for the data as a separate table. This allows the use of a donor file
  vector<long> setColumnIndices(long filetype);//< Sets thecolumn reordering index vector based on the predefined maincolnames. This is different to the table colnames if using an id file
  bool readTable(string filename,vector< vector <double> > *dest,vector<long> *indexes,long filetype);///< Reads a table from a file
  bool addToChunkCounts(vector< vector <double> > *tmpptr,vector<long> *indexes);
  bool addToRegionChunkCounts(vector< vector <double> > *tmpptr,vector<long> *indexes);
  bool addToRegionSquaredChunkCounts(vector< vector <double> > *tmpptr,vector<long> *indexes);
  bool addToChunkLengths(vector< vector <double> > *tmpptr,vector<long> *indexes);
  bool addToMutationProbs(vector< vector <double> > *tmpptr,vector<long> *indexes);
  
  void makeRegionChunkCounts(vector< vector <double> > *tmpptr,vector< vector <double> > *chunkcounts);///< makes tmpptr into a fake regionchunk inputfile assuming the whole of each data is a region
  void makeRegionSquaredChunkCounts(vector< vector <double> > *tmpptr,vector< vector <double> > *chunkcounts);///< makes tmpptr into a fake regionsquarechunk inputfile assuming the whole of each data is a region
  long readFilesOneRoot(string root);///< tries to read the root set of files and add them to the current matrices (creating them if necessary)
  long readFiles();///< Assumes that the files have been rationalised, reads into returns 1 on success
  long validateMatrices();///< Checks that the matrices all match up.  Returns -(index) of the failed matrix
  ////////////////////////////////////////////
  // PROCESSING DATA
  void calcTheoreticalVar();///<Calculates the theoretical variance
  void calcEmpiricalVar();///<Calculates the empirical variance
  void calcTheoreticalVarForce();///<Calculates the theoretical variance (applying ther force condition)
  void calcEmpiricalVarForce();///<Calculates the empirical variance (applying ther force condition)
  double calcC(bool useforce=false);///< returns the value of C based on the current data
  void finaliseData();///< Processes the data into output form, i.e. do the calculations for mutations and lengths

  vector<long> getRowOrder();///<Gets the row order based on the column name order
  void reorderMatrix(vector<vector<double> > *mat,vector<long> *torder);///<Reorders the matrix *rows* according to torder
  void reorderVector(vector<long > *mat,vector<long> *torder);///< reorder mat based on torder
  void reorderVector(vector<string > *mat,vector<long> *torder);///< reorder mat based on torder
  void sortRows();///<Sorts all the rows based on the column name order
  ////////////////////////////////////////////
  // WRITING OUTPUT 
  bool writeRowName(std::ostream * os,long rownum,long filetype);///< Writes the row name
  bool writeHeaderRow(std::ostream * os,long filetype);///< Writes the column headers 
  void writematrix(std::ostream * os,vector<vector<double> > *mat,long filetype);///<Writes a whole matrix (including header)
  bool writeChunkCountFile();///< Writes the current chunk file
  bool writeRegionChunkCountFile();///< Writes the current chunk file
  bool writeRegionChunkCountSqFile();///< Writes the current chunk file
  bool writeChunkLengthsFile();///< Writes the current chunk file
  bool writeMutationProbFile();///< Writes the current chunk file
  bool writeFiles();///< Tries to write the outputfiles as summed
};

#endif
