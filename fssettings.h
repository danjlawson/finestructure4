#ifndef FSSETTINGS_H
#define FSSETTINGS_H


#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <fstream>

using namespace std;


class FsSettingsValue
{ ///< A class for holding settings read from file
 public:
  FsSettingsValue(std::string line=""); // set from a string in the form (name:arg1 arg2 arg3), or create a black setting object
  FsSettingsValue(std::vector<std::string> lineargs,bool throwerrors=false); // set from a vector in the form (name arg1 arg2...)
  ~FsSettingsValue();
  void setFromCmdArgs(std::vector<std::string> lineargs,bool throwerrors); 
  void set(std::string line,bool throwerrors=false); // set the setting value
  void set(std::string name, std::vector<std::string> lineargs,bool throwerrors=false);// set from a command line like string from a separated (name value) pair
  void set(std::vector<std::string> lineargs,bool throwerrors=false);// set from a command line like string (name:val)
  void setSection(std::string line);// assign it a section
  void clear();// remove all data from the setting
  std::string getName(); // return the name 
  std::string getSection(); // return the section
  std::string getLine(); // return the original string
  std::string getEndedLine(); // return the original string (with \n)

  std::string getVal(); // return the value (as a string)
  int getValAsInt(); // return the value
  double getValAsDouble(); // return the value
  std::vector<std::string> getValAsStringVec(char c=' '); // return the value
  std::vector<int> getValAsIntVec();///< return the value
  std::vector<double> getValAsDoubleVec();///< return the value
  bool success; // whether it contains a read value
 protected:
  std::string section;
  std::string name;
  std::string val;
  std::string line;
};

class FsSettingsReader
{
 public:
  FsSettingsReader(std::string f,bool verbose);///< Create from a file
  ~FsSettingsReader();///< Destructor
  void openFile();/// < open the file
  FsSettingsValue getNext();///< Gets the next settings value

 protected:
  
  bool verbose;///< Whether we output anything
  string filename;///< The file name
  std::ifstream fs;///< The stream from which we read the file
  string line; //< previously read line
  string section; //< the section we are currently in
};

#endif
