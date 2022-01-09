#ifndef FSUTILS_H
#define FSUTILS_H


#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <fstream>

using namespace std;

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems); ///< Split s into the vector elems
std::vector<std::string> split(const std::string &s, char delim); // Split s into a vector by delim
std::vector<std::string> split(const std::string &s, std::string delims); // Split s into a vector by delim
std::vector<std::string> getcommands(std::vector<std::string> svec); /// Converts commands split on spaces to commands split also on commas

std::string getcwd(); // Get the current working directory
int directoryExists(std::string pathname);// return 1 if directory exists, 0 if not (-1 if exists and isn't a directory)
void deleteFolderTree (string directory_name);// delete a folder tree and all files in them
void ensureDirectory(std::string pathname); // Create a directory, only if needed
  
std::string ssvec(std::vector<std::string> vec);///< Convert a string vector into a comma separated string
std::string ssvec(std::vector<int> vec);///< Convert an int vector into a comma separated string
std::string ssvec(std::vector<double> vec);/// < Convert a double vector into a comma separated string
std::istream& ignoreline(std::ifstream& in, std::ifstream::pos_type& pos);  // utility fn
std::string getLastLine(std::ifstream& in);// gets the last line in the file (as an ifstring)
std::string getLastLine(string in); // gets the last line in the file (for a named file)
std::string getLineContaining(std::string val, std::ifstream& in);// get the first line containing text val (return empty string if not found)
std::string getLineContaining(std::string val, string in);// get the first line containing text val from file named "in" (return empty string if not found)

std::string trim(const std::string& str,
                 const std::string& whitespace = " \t");///< trim whitespace from the front and back of a string
std::string reduce(const std::string& str,
                   const std::string& fill = " ",
                   const std::string& whitespace = " \t"); ///< trim ALL whitespace into the "fill" character

std::vector<std::string> unique(std::vector<std::string> in);///< Returns a unique version of the input

std::vector<std::string> getIdsFromFile(std::string  filename,bool keepall=true,int column=0);///< Extracts ids (or pops) from chromopainter idfile
std::vector<int> getChromoPainterHeaderInfo(std::string filename,int ploidy); ///< Extracts header information from the chromopainter header; nhaps,nsnps, cpv2?

std::string projectroot(std::string in); ///< Remove any ".cp" ending 
std::string projectfull(std::string inroot); ///< Add an ".cp" ending 
std::string baseName(std::string in);///< remove directories in string filename
std::string dirName(std::string in);///< Retain only the directory in a string filename

std::vector<int> sampleVec(std::vector<int> X,int m);///< samples m random elements of X
int countLines(std::string file);///< count ther number of lines in a file. Returns negative for file not existing

bool BothAreSpaces(char lhs, char rhs);///< A simple comparison for unique to remove duplicated spaces


class ParallelStream{
    std::ostringstream stdStream;
public:
    ParallelStream(){}
    template <class T>
    ParallelStream& operator<<(const T& inData){
        stdStream << inData;
        return *this;
    }
    std::string toString() const{
        return stdStream.str();
    }
};

bool isNumericallySorted(std::vector<string> in);///< test whether all integers from the string are sorted
int stringToInt(string in); // read an integer, throw an exception if the string is invalid
int forceStringToInt(string in);///< read integer from only the numeric values in the string
#endif
