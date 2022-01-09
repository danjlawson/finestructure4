#ifndef FSCMDS_H
#define FSCMDS_H

#include <vector>
#include <string>
using namespace std;

/** 
    @brief Parameters that can be set manually on the command line. This probably should be a static class; these are not supposed to be modified
*/
class FsPar
{
 public:
  FsPar(std::string name, int stage,std::string help,int type=0);///< Create a parameter. 
  ~FsPar();
  inline std::string getName(){
    return(name);
  };
  inline std::string getHelp(){
    return(help);
  };
  inline int getStage(){
    return(stage);
  };
  inline int getType(){
    return(type);
  };
private:
  std::string name;///< the name of the parameter
  std::string help;///< the free text field
  int stage;///< stage is the maximum stage it can be set (-1 for any time)
  int type;///< type is 0 for normal parameters and 1 for derived parameters that shouldn't normally be set by command line (but can be)
};

class FsCmd
{
 public:
  FsCmd(std::string name, int stage,std::string shortargs,std::string help);///< Create a command. 
  ~FsCmd();
  inline std::string getName(){
    return(name);
  };
  inline std::string getShortargs(){
    return(shortargs);
  };
  inline std::string getHelp(){
    return(help);
  };
  inline int getStage(){
    return(stage);
  };
private:
  std::string name;///< the name of the command
  std::string shortargs;///< short version of the arguments
  std::string help;///< the free text field
  int stage;///< stage is the maximum stage it can be run (-1 for any time)
};

#endif
