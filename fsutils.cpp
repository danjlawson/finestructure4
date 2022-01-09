

#include <algorithm>
#include <limits>
#include <cstdio>
#include <string>
#include <exception>
#include <stdexcept>
#include <fstream>
#include <unistd.h>

/////////////////////////////
// Not cross platform?
#include <sys/types.h>
#include <sys/stat.h>

/////////////////////////////

#include "fsutils.h"

using namespace std;

/////////////////////////////////////////////////
// Split s into the vector elems
bool BothAreSpaces(char lhs, char rhs) {
  return (lhs == rhs) && (lhs == ' ');
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}
std::vector<std::string> split(const std::string &s, string delims) {
    std::string item;
    std::vector<std::string> elems;
    size_t lastpos=0;
    size_t pos=0;
    bool finished=0;
    while (!finished) {
      pos=s.find_first_of(delims,lastpos);
      item=s.substr(lastpos,pos-lastpos);
      lastpos=pos+1;
      elems.push_back(item);
      if(pos==string::npos) finished=true;
    }
    return elems;
}

/////////////////////////////////////////////////////////////
// Split s into a vector by delim
std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}
/////////////////////////////////////////////////////////////
std::vector<std::string> getcommands(std::vector<std::string> svec) 
{ /// Converts commands split on spaces to commands split also on commas
  unsigned int i=0;
  while(i<svec.size()){
    std::vector<std::string> tvec=  split(svec[i],',');
    if(tvec.size()>1) { // replace the original arguments with the split
      svec.erase(svec.begin()+i);
      svec.insert(svec.begin()+i,tvec.begin(),tvec.end());
      i+=tvec.size();
    }else i++;
  }
  return(svec);
}
/////////////////////////////////////////////////
std::string getcwd() {// Get the current working directory
    char buf[FILENAME_MAX];
    char* succ = getcwd(buf, FILENAME_MAX);
    if( succ ) return std::string(succ);
    return "";  // raise a flag, throw an exception, ...
}

int directoryExists(std::string pathname){ // Check if name is a directory. If yes, return 1. If it doesn't exist, return 0. If it exists but is not a directory, return -1.
  struct stat info;
  if( stat( pathname.c_str(), &info ) != 0 ) // doesn't exsist
    return(0);
  else if( info.st_mode & S_IFDIR )  // S_ISDIR() doesn't exist on my windows 
    return(1); // does exist
  else
    return(-1); // exists, is not a directory
}

void deleteFolderTree (string directory_name) {
  // We rely on a system call. Its not used much so the slowness doesn't matter.
  stringstream cmd;
  cmd<<"rm -r "<<directory_name;
  system(cmd.str().c_str());
}

void ensureDirectory(std::string pathname){
  // Check if a directory exists. If not, create it. Note: Not recursive
  int already_exists=directoryExists(pathname);
  if(already_exists<0) 	throw(runtime_error(string("Cannot create directory ").append(pathname).append(" as it already exists and is not a directory.")));
  if(!directoryExists(pathname))   mkdir(pathname.c_str(),  (mode_t)0777);
}

/////////////////////////////////////////////////
///////////////////////////////////////////////// 
// Tools for checking chromopainter output
std::istream& ignoreline(std::ifstream& in, std::ifstream::pos_type& pos)
{ // utility fn
    pos = in.tellg();
    return in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}
std::string getLastLine(std::ifstream& in)
{// gets the last line in the file
    std::ifstream::pos_type pos = in.tellg();

    std::ifstream::pos_type lastPos;
    while (in >> std::ws && ignoreline(in, lastPos))
        pos = lastPos;

    in.clear();
    in.seekg(pos);

    std::string line;
    std::getline(in, line);
    return line;
}
std::string getLastLine(string in)
{// gets the last line in the file
  std::ifstream file(in.c_str());
  return(getLastLine(file));
}
std::string getLineContaining(std::string val, std::ifstream& in)
{// get the next line containing a value
    string line;
    size_t pos;
    while(getline(in, line)){
      pos=line.find(val);
      if(pos!=string::npos) return(line);
    }
    return string("");
}
std::string getLineContaining(std::string val, string in)
{// gets the last line in the file
  std::ifstream file(in.c_str());
  return(getLineContaining(val,file));
}


/////////////////
std::string trim(const std::string& str,
                 const std::string& whitespace)
{
    const size_t strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const size_t strEnd = str.find_last_not_of(whitespace);
    const size_t strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

std::string reduce(const std::string& str,
                   const std::string& fill,
                   const std::string& whitespace)
{
    // trim first
    std::string result = trim(str, whitespace);

    // replace sub ranges
    size_t beginSpace = result.find_first_of(whitespace);
    while (beginSpace != std::string::npos)
    {
        const size_t endSpace = result.find_first_not_of(whitespace, beginSpace);
        const size_t range = endSpace - beginSpace;

        result.replace(beginSpace, range, fill);

        const size_t newStart = beginSpace + fill.length();
        beginSpace = result.find_first_of(whitespace, newStart);
    }

    return result;
}

std::vector<std::string> unique(std::vector<std::string> in) 
{

  std::vector<string> result;
  result.reserve(in.size());
  for (std::vector<string>::iterator itor = in.begin();
                                    itor != in.end(); ++itor)
  {
    if (std::find(result.begin(), result.end(), *itor) == result.end()){
            result.push_back(*itor);
    }
  }
  return result;
}

std::string ssvec(std::vector<std::string> vec)
{
  std::ostringstream ss;
  if(vec.size()>0) ss<<vec[0];
  for(unsigned int i=1;i<vec.size();i++) ss<<","<<vec[i];
  return(ss.str());
}

std::string ssvec(std::vector<int> vec)
{
  std::ostringstream ss;
  if(vec.size()>0) ss<<vec[0];
  for(unsigned int i=1;i<vec.size();i++) ss<<","<<vec[i];
  return(ss.str());
}

std::string ssvec(std::vector<double> vec)
{
  std::ostringstream ss;
  if(vec.size()>0) ss<<vec[0];
  for(unsigned int i=1;i<vec.size();i++) ss<<","<<vec[i];
  return(ss.str());
}

std::vector<std::string> getIdsFromFile(std::string  filename,bool keepall,int column){
  // Gets the ids from a file
  // either keeping all of them (keepall) or only those that will appear in a ChromoPainter analysis (keepall=FALSE)
  vector<string>ret;
  vector<string>tline;
  std::ifstream file(filename.c_str());
  string line;
  while(getline(file, line)) {
    tline=split(line," \t");
    if(!keepall && tline.size()>1){
      bool tuse=true;
      if(tline.size()>2) {
	istringstream ( tline[2] ) >> tuse;
      }
      if(!tuse) continue;
    }
    if(tline.size()>0) {
      if(column>=(int)tline.size()){
	throw(runtime_error("donors missing"));
      }
      if(column==2 && tline.size()<3){ // Always included if we don't have an inclusion column
	ret.push_back(string("1"));
      }else{
	ret.push_back(tline[column]);
      }
    }
  }
  file.close();
  if(ret.size()==0) {
    cerr<<"ERROR: idfile \""<<filename.c_str()<<"\" could not be read!"<<endl;
    throw(runtime_error("data error"));
  }
  return(ret);
}

std::vector<int> getChromoPainterHeaderInfo(std::string filename,int ploidy)
{
  vector<int>ret;
  int nsnps=0;
  int nhaps=0;
  bool usev2=true;
  
  std::ifstream file(filename.c_str());
  vector<string> fvals;
  bool done=false;

  while(!done){
    std::string line;
    if(!getline(file, line)) {
      done=true;
      cerr<<"ERROR: phase file "<<filename.c_str()<<" could not be correctly read!"<<endl;
      throw(runtime_error("phase problem"));
    }
    if(line.at(0)=='P') {
      istringstream ( fvals.at(fvals.size()-2) ) >> nhaps;
      istringstream ( fvals.at(fvals.size()-1) ) >> nsnps;
      if(fvals.size()>2) {
	usev2=false;
	nhaps*=ploidy; /// correction for old method
      }
      done=true;
    }else{
      fvals.push_back(line);
    }
  }
  
  ret.push_back(nhaps);
  ret.push_back(nsnps);
  ret.push_back(usev2);
  return(ret);
}

///////////////////////
std::string projectroot(std::string in){
  size_t found=in.find(".cp");
  if(found!=std::string::npos){
    in=in.substr(0,found);
  }
  return(in);
 }

std::string projectfull(std::string inroot){
  inroot.append(".cp");
  return(inroot);
}

std::string baseName(std::string in){
  size_t found=in.find_last_of('/');
  if(found!=std::string::npos){
    in=in.substr(found+1,std::string::npos);
  }
  return(in);
}

std::string dirName(std::string in){
  size_t found=in.find_last_of('/');
  if(found!=std::string::npos){
    in=in.substr(0,found+1);
    return(in);
  }
  return(string(""));
}

std::vector<int> sampleVec(std::vector<int> X,int m){
  int n=X.size();
  for (int i = 0; i < m; i++)
    swap(X[i],X[i+(rand()%(n-i))]);
  X.resize(m);
  return(X);
}

int countLines(std::string file){
  if( access( file.c_str(), F_OK ) == -1 ) { // file does not exist
    return(-1);
  }

  int numlines=0;
  std::ifstream in(file.c_str());

  std::string unused;
  while ( std::getline(in, unused) )    ++numlines;
  return numlines;
  
}
  
bool isNumericallySorted(std::vector<string> in){
  std::vector<int> out;
  for(unsigned int c1=0;c1<in.size();++c1){
    out.push_back(forceStringToInt(in[c1]));
  }
  for(unsigned int c1=1;c1<out.size();++c1){
    if(out[c1-1]>out[c1]) return(false);
  }
  return(true);
}
  
int stringToInt(string in){
  int ret;
  istringstream ( in ) >> ret;
  stringstream ss;
  ss<<ret;
  if(ss.str().compare(in)!=0){
    throw(runtime_error("invalid integer"));
  }
  return(ret);
}

int forceStringToInt(string in){
  int out;
  ostringstream digits;
  for (unsigned int c1=0; c1 < in.length(); c1++)
    if (isdigit(in[c1])) digits<<in[c1];
  istringstream(digits.str())>>out;
  return(out);
}
