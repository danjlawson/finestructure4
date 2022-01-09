

#include <iostream>
#include <limits>
#include <cstdio>
#include <string>
#include <stdexcept>
#include <exception>
#include <fstream>
#include <unistd.h>

/////////////////////////////
// Not cross platform?
#include <sys/types.h>
#include <sys/stat.h>

/////////////////////////////

#include "fsutils.h"
#include "fssettings.h"

using namespace std;

FsSettingsReader::FsSettingsReader(std::string f,bool verbose)
{
  filename=f;
  this->verbose=verbose;
  openFile();
}

FsSettingsReader::~FsSettingsReader() 
{
  fs.close();
}

void FsSettingsReader::openFile()
{
  fs.open(filename.c_str());
  if (!fs.is_open()) { // something went wrong!
    cerr<<"Could not open settings file "<<filename<<" for reading. Does it exist?"<<endl;
    throw(runtime_error("fssettings read error"));
  }
  if(verbose){
    cout<<"Opened settings file "<<filename<<endl;
  }
}

FsSettingsValue FsSettingsReader::getNext()
{
  FsSettingsValue ret;
  while ( getline (fs,line) ) {
    //    cout << "DEBUG: || " << line << " ||"<<endl;
    line=reduce(line.substr(0,line.find("#")));
    //    cout << "DEBUG: || " << line << " ||"<<endl;
    ret.set(line);
    if(ret.success){
      if(ret.getName().compare("section")==0){
	section=ret.getVal();
	ret.clear();
      }else {
	ret.setSection(section);
	return(ret);
      }
    }
  }
  return(ret);
}

///////////////////
// Setting values
FsSettingsValue::FsSettingsValue(std::string line)
{
  success=false;
  
  set(line);
}

FsSettingsValue::FsSettingsValue(std::vector<std::string> lineargs,bool throwerrors)
{
    vector<string> nargs=lineargs;
    nargs.erase(nargs.begin());
    set(lineargs[0],nargs,throwerrors);
}

FsSettingsValue::~FsSettingsValue()
{
}

void FsSettingsValue::clear()
{
  success=false;
  section=name=val=line="";
}
 

void FsSettingsValue::set(std::string name,std::vector<std::string> lineargs,bool throwerrors)
{
  ostringstream ss;
  ss<<name<<":";
  for(unsigned int c1=0;c1<lineargs.size();c1++){
    ss<<lineargs[c1];
    if(c1<lineargs.size()-1) ss<<" ";
  }
  set(ss.str(),throwerrors);
}

void FsSettingsValue::set(std::vector<std::string> lineargs,bool throwerrors)
{
  ostringstream ss;
  for(unsigned int c1=0;c1<lineargs.size();c1++){
    ss<<lineargs[c1]<<" ";
  }
  set(ss.str(),throwerrors);
}

void FsSettingsValue::set(std::string line,bool throwerrors)
{
  if(line.compare("")==0) {
    if(throwerrors) throw(runtime_error(line));
    return;
  }
  std::vector<std::string> splitval=split(line,':');
  if(splitval.size()==1 || splitval.size()>2){
    if(throwerrors) throw(runtime_error(line));
    return;
  }
  this->line=line;
  name=splitval[0];
  val=splitval[1];
  while(name.substr(0,1).compare("-")==0){ // strip leading -
    name=name.substr(1,name.size());
  }
  success=true;
}

void FsSettingsValue::setSection(std::string val)
{
  section=val;
}

std::string FsSettingsValue::getSection()
{
  return(section);
}

std::string FsSettingsValue::getName()
{
  return(name);
}

std::string FsSettingsValue::getLine()
{
  return(line);
}

std::string FsSettingsValue::getEndedLine()
{
  string s=line;
  s.append("\n");
  return(s);
}

std::string FsSettingsValue::getVal()
{
  return(val);
}

int FsSettingsValue::getValAsInt()
{
  int ret;
  istringstream ( val ) >> ret;
  return(ret);
}

double FsSettingsValue::getValAsDouble()
{
  double ret;
  istringstream ( val ) >> ret;
  return(ret);
}

std::vector<std::string> FsSettingsValue::getValAsStringVec(char c)
{
  std::vector<std::string> ret=split(val,c);
  ret=getcommands(ret);
  for(unsigned int c1=0;c1<ret.size();c1++){
    ret[c1]=trim(ret[c1]);
    //    cout<<"DEBUG: Ret "<<c1<<" = "<<ret[c1]<<endl;
  }
  return(ret);
}

std::vector<int> FsSettingsValue::getValAsIntVec()
{
  std::vector<std::string> rets=split(val,' ');
  rets=getcommands(rets);
  std::vector<int> ret;
  int v;
  for(unsigned int c1=0;c1<rets.size();c1++){
    rets[c1]=trim(rets[c1]);
    istringstream ( rets[c1] ) >> v;
    ret.push_back(v);
  }
  return(ret);
}

std::vector<double> FsSettingsValue::getValAsDoubleVec()
{
  std::vector<std::string> rets=split(val,' ');
  rets=getcommands(rets);
  std::vector<double> ret;
  double v;
  for(unsigned int c1=0;c1<rets.size();c1++){
    rets[c1]=trim(rets[c1]);
    istringstream ( rets[c1] ) >> v;
    ret.push_back(v);
  }
  return(ret);
}
