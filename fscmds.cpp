#include "fscmds.h"

FsPar::FsPar(std::string name, int stage,std::string help,int type){
  //Create a parameter.  stage is the maximum stage it can be set (-1 for any time), help is the free text info for this parameter. type is 0 for normal parameters and 1 for derived parameters that shouldn't normally be set by command line (but can be)
  this->name=name;
  this->stage=stage;
  this->help=help;
  this->type=type;
}

FsPar::~FsPar(){
}

FsCmd::FsCmd(std::string name, int stage,std::string shortargs,std::string help){
  //Create a command.  stage is the maximum stage it can be run (-1 for any time), help is the free text info for this command.
  this->name=name;
  this->stage=stage;
  this->shortargs=shortargs;
  this->help=help;
}

FsCmd::~FsCmd(){
}
