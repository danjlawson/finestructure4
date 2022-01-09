#include <sys/stat.h>
#include <fstream>
#include <sstream>
#include <cstring>
#include <stdlib.h>
#include <dirent.h>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "ChromoCombineFuns.h"

////////////////////////////////////////////
// DEFINITIONS AND SETUP

ChromoCombineWorker::ChromoCombineWorker(){
}

ChromoCombineWorker::~ChromoCombineWorker(){
}

void ChromoCombineWorker::makeEndings(){
  requiredendings.clear();
  possiblemiddles.clear();
  possibleendings.clear();

  possiblemiddles.push_back(string(".chunkcounts"));
  possiblemiddles.push_back(string(".regionchunkcounts"));
  possiblemiddles.push_back(string(".regionsquaredchunkcounts"));
  possiblemiddles.push_back(string(".chunklengths"));
  possiblemiddles.push_back(string(".mutationprobs"));
  possiblemiddles.push_back(string(".copyprobsperlocus"));  // .gz goes at the end
  possiblemiddles.push_back(string(".copyprobsperlocus"));
  possiblemiddles.push_back(string(".EMprobs"));
  possiblemiddles.push_back(string(".samples"));
  possiblemiddles.push_back(string(".prop"));

  if(wantcounts){
    requiredendings.push_back(string(".chunkcounts").append(outputfileending));
    requiredendings.push_back(string(".regionchunkcounts").append(outputfileending));
    requiredendings.push_back(string(".regionsquaredchunkcounts").append(outputfileending));
  }
  if(wantlengths) requiredendings.push_back(string(".chunklengths").append(outputfileending));
  if(wantmuts) requiredendings.push_back(string(".mutationprobs").append(outputfileending));

  for(unsigned long c1=0;c1<possiblemiddles.size();c1++){
    possibleendings.push_back(possiblemiddles[c1]);
    possibleendings.back().append(outputfileending);
  }
  possibleendings[5].append(".gz");
  cval=-1;// initialise c to *something*; -ve so we can tell its not calculated properly
}

void ChromoCombineWorker::clearAllData()
{
  rownames.clear();
  colnames.clear();
  topleftentry="";

  counts.clear();
  regioncounts.clear();
  regioncountssq.clear();
  lengths.clear();
  muts.clear();

  indexcounts.clear();

  regionsums.clear();
  empiricalvar.clear();
  theoreticalvar.clear();

  forceconts.clear();
  forceinds.clear();
  forceindsid.clear();
  contin.clear();
  numeffinds=0;

  forcedregioncounts.clear();
  forcedregioncountssq.clear();
  forcedregioncountssumsq.clear();
  forcedregionsums.clear();

  cval=-1;
  completegenomesnormfactor=0;
  completegenomesnormfactorsq=0;
  completegenomesnormfactorcount=0;
}

////////////////////////////////////////////
// ORGANISING INPUT FILES

void ChromoCombineWorker::getFilesFromDirectory(string dir,vector<string> *dirlisting){
  struct dirent *entry;
  DIR *dp;

  dp = opendir(dir.c_str());
  if (dp == NULL) {
    cerr<<"Failed to open directory"<<endl;
    throw(runtime_error(string("Invalid directory: ").append(dir)));
  }
  string tmp(dir);
    if(tmp[tmp.length()-1]!='\\' && tmp[tmp.length()-1]!='/'){
#if defined(__WXMSW__)
	tmp.append("\\");
#else
	tmp.append("/");
#endif
    }
  while((entry = readdir(dp))) {
    dirlisting->push_back(string(tmp).append(entry->d_name));
  }
  closedir(dp);
}

bool ChromoCombineWorker::searchForEndings(string str){
  size_t res=0;
  for(unsigned long c1=0;c1<possibleendings.size();c1++){
    res=str.find(possibleendings[c1]);
    if(res!=string::npos){// found the phrase
      // now check it is at the end
      if(str.size()==res+possibleendings[c1].size()) {return(true);
      }else return(false);
    }
  }
  return(false);
}

void ChromoCombineWorker::getInterestingFiles(vector<string> allfilenames){
  inputfilerootraw.clear();
  for(unsigned long c1=0;c1<allfilenames.size();c1++){
    if(searchForEndings(allfilenames[c1])) inputfilerootraw.push_back(allfilenames[c1]);
  }
}

void ChromoCombineWorker::getFilesFromDirectories(){
  vector<string> allfilenames;
  for(unsigned long c1=0;c1<inputfilerootraw.size();c1++) {
    getFilesFromDirectory(inputfilerootraw[c1], &allfilenames);
  }
  getInterestingFiles(allfilenames);
}

vector<string> ChromoCombineWorker::getFileList(string root,vector<string> endings){
  vector<string> ret;
  for(unsigned long endon=0;endon<endings.size();endon++){
    ret.push_back(root);
    ret.back().append(endings[endon]);
  }
  return(ret);
}

void ChromoCombineWorker::testFilesExist(){
  if(validfile.size()!=inputfileroot.size()) {
		cerr<<"Error: validated "<<validfile.size()<<" files, but have "<<inputfileroot.size()<<" to check!"<<endl;
		throw(runtime_error("Validation problem"));
  }
  for(unsigned long fron=0;fron<inputfileroot.size();fron++){
    validfile[fron]=0;
    vector<string> tlist=getFileList(inputfileroot[fron],requiredendings);
    //    if(!wantcounts)tlist.erase(tlist.begin(),tlist.begin()+3);
    for(unsigned long endon=0;endon<tlist.size();endon++){
      cout<<"testing file "<<tlist[endon].c_str()<<endl;
		long intStat;
		struct stat stFileInfo;
		// Attempt to get the file attributes
		intStat = stat(tlist[endon].c_str(),&stFileInfo);
		if(intStat != 0) validfile[fron]++; //failed to read this file
    }
  }
}

string ChromoCombineWorker::getRoot(string in){
  string ret=in;
  size_t res=0;
  for(unsigned long c1=0;c1<possibleendings.size();c1++){
    res=0;
    while(res!=string::npos){
      res=ret.find(possibleendings[c1]);
      if(res!=string::npos){
	ret.replace(res,possibleendings[c1].length(),"");
      }
    }
  }
  return(ret);
}

bool ChromoCombineWorker::rationaliseInputFiles(){
  inputfileroot.clear();
  for(unsigned long c1=0;c1<inputfilerootraw.size();c1++) {
    bool found=false;
    string strtest=getRoot(inputfilerootraw[c1]);
    for(unsigned long c2=0;c2<inputfileroot.size();c2++) {
      if(strtest.compare(inputfileroot[c2])==0){
	found=true;
	c2=inputfileroot.size();
      }
    }
    if(!found){
      inputfileroot.push_back(strtest);
      validfile.push_back(-1);
    }
  }
  return(true);
}

////////////////////////////////////////////
// FORCE FILES, READING AND APPLYING
void ChromoCombineWorker::addSuper(string name,vector<string> ind){
  for(long i=0;i<(long)forceinds.size();i++){ // check for duplication
    if(name==forceconts[i]){
      cerr<<"Error population "<<forceconts[i]<<" already exists, but attempted to add it again. Check whether it is duplicated in your -f forcefile"<<endl;
	  throw(runtime_error("Fixed population construction error"));    }
    for(long j=0;j<(long)forceinds[i].size();j++){
      for(long k=0;k<(long)ind.size();k++){
	if(ind[k] == forceinds[i][j]){
	  cerr<<"Error in population "<<name<<": adding individual "<<ind[k]<<" who already exists in population "<<forceconts[i]<<endl;
	  throw(runtime_error("Fixed population construction error"));
	}
      }
    }
  }
  forceconts.push_back(name);
  forceinds.push_back(ind);
}

void ChromoCombineWorker::makeSuper(std::string superstring){
  superstring.erase (std::remove (superstring.begin(), superstring.end(), ' '), superstring.end());
  superstring.erase (std::remove (superstring.begin(), superstring.end(), '\t'), superstring.end());
  superstring.erase (std::remove (superstring.begin(), superstring.end(), '"'), superstring.end());
  superstring.erase (std::remove (superstring.begin(), superstring.end(), '\n'), superstring.end()); // shouldn't be required as \n is newline
  superstring.erase (std::remove (superstring.begin(), superstring.end(), '\r'), superstring.end()); // saves us when moving reading dos files in unix

  long gnum=0;
  vector<string> ind;
  string dname;
  const char* test = ",()";
  size_t found=superstring.find_first_of(test), pos=0;
  while(found!=string::npos){
    if(superstring.at(pos)==')'){// should be a name, or the end
      addSuper(dname,ind);
      ind.clear();
      dname=superstring.substr(pos+1,found-pos-1);
    }else if(pos==0){// should be a name
      dname=superstring.substr(pos,found-pos);
    }
    if(superstring.at(pos)=='(') {// new population
      if(dname.size()==0){
	stringstream ss;
	ss<<string("Group")<<gnum;
	dname=ss.str();
      }
      ind.clear();
      //	    cout<<"NAME:"<<dname<<endl;
      if(found==0) found=superstring.find_first_of(test,pos+1);
      if (strchr(test, superstring.at(pos+1)) != NULL) {
	cerr<<"Invalid fixed file structure: two control characters appear together (population: \""<<dname<<"\")."<<endl;
	throw(runtime_error("State fixed file error: name not found!"));
      }	  
    }
    if(superstring.at(pos)=='(' || superstring.at(pos)==',') {
      if (strchr(test, superstring.at(pos+1)) != NULL) {
	cerr<<"Invalid fixed file structure: two control characters appear together (population: \""<<dname<<"\")."<<endl;
	throw(runtime_error("State fixed file error: name not found!"));
      }
      ind.push_back(superstring.substr(pos+1,found-pos-1));
    }
    pos=found;
    if(pos!=string::npos) found=superstring.find_first_of(test,pos+1);
  }
  if(ind.size()>0) addSuper(dname,ind);
}

void ChromoCombineWorker::makeSuperFromFile(std::string filename){
      string line, superstring;
      ifstream file;
      file.open(filename.data());//Open file
      if(!file.good()) throw(runtime_error("Invalid super individual file!"));
      while(1){
	safeGetline(file,line);
	if(line[0]!='#')superstring.append(line); // Ignore comment lines
	if (file.eof())
	      break;//Stop if end of file
      }
      file.close();
      makeSuper(superstring);
}

void ChromoCombineWorker::getForceIds(){
  contin=vector<long>(rownames.size(),-1);
  for(unsigned long c1=0;c1<forceinds.size();c1++){
    forceindsid.push_back(vector<long>(forceinds[c1].size(),0));
    for(unsigned long c2=0;c2<forceinds[c1].size();c2++){
      bool found=false;
      for(unsigned long rownum=0;rownum<rownames.size();rownum++){
	if(forceinds[c1][c2].compare(rownames[rownum])==0) {
	  found=true;
	  forceindsid[c1][c2]=rownum;
	  contin[rownum]=c1;
	  rownum=rownames.size();
	  numeffinds--;
	}
      }
      if(!found){
	  cerr<<"Error in force file: individual name: "<<forceinds[c1][c2]<< " not found!"<<endl;
	  cerr<<"Options are:";
	  for(unsigned long c3=0;c3<rownames.size();c3++) cerr<<rownames[c3]<<",";
	  cerr<<endl;
	  throw(runtime_error(string("Error in force file: individual name: ").append(forceinds[c1][c2]).append(" not found in count file!")));
      }
    }
  }
  long indexon=-(long)forceinds.size();
  for(unsigned long c1=0;c1<contin.size();c1++){
    if(contin[c1]<0) contin[c1] = indexon--;
  }
  numeffinds+=forceinds.size();
}

void ChromoCombineWorker::applyForceFile(){
  forcedregioncounts=vector<vector<double> >(numeffinds,vector<double>(numeffinds,0));
  forcedregioncountssq=forcedregioncounts;
  forcedregioncountssumsq=forcedregioncounts;
  forcedregionsums=vector<double>(numeffinds,0);
  for(unsigned long c1=0;c1<counts.size();c1++) {
      long index1=abs(contin[c1]);
      if(index1<0 || index1>=(long)forcedregioncounts.size()){
	cerr<<"Error: invalid index1. Forced index: "<< index1<<" >= maximum"<<forcedregioncounts.size()<<endl;
	throw(runtime_error("Error in applyForceFile"));
      }
      forcedregionsums[index1]+=regionsums[c1];
    for(unsigned long c2=0;c2<counts[c1].size();c2++) {
      long index2=abs(contin[c2]);
      if(index2<0 || index2>=(long)forcedregioncounts[index1].size()){
	cerr<<"Error: invalid index2 "<<c2<<", size: "<<index2<<" >= maximum "<<(long)forcedregioncounts[index1].size()<<endl;
	throw(runtime_error("Error in applyForceFile"));
      }
      forcedregioncountssq[index1][index2]+=regioncountssq[c1][c2+1];
      forcedregioncounts[index1][index2]+=regioncounts[c1][c2+1];
      forcedregioncountssumsq[index1][index2]+=regioncounts[c1][c2+1]*regioncounts[c1][c2+1];
    }
  }
}

bool ChromoCombineWorker::writeForceFile(double forceC) {
      vector<string> oldforcefile;
      string line;
      string outfile=forcefileoutput;
      if(outfile.size()==0)outfile=forcefile;
      if(verbose>0) cout<<"Writing force file to "<<outfile<<endl;
      // read in the old file
      ifstream file;
      file.open(forcefile.data());//Open file
      if(!file.good()) return(false);
      while(1){
	safeGetline(file,line);
	oldforcefile.push_back(line); // Ignore comment lines
	if (file.eof())
	      break;//Stop if end of file
      }
      file.close();
// write the updated file
  filebuf fb;
  try{
    fb.open (outfile.c_str(),ios::out);
    }catch(exception &x){
    cerr<<"Error opening file!"<<endl<<x.what()<<endl; return false;}
  ostream os (&fb);
  os<<setprecision(9);
  os<<"#Cfactor "<<forceC<<endl;
  os<<"#ChomoCombine: Forcefile generated for dataset root: "<<outputfileroot<<" with ending: "<<outputfileoutending<<endl;
  os<<"#ChomoCombine: For the following datasets, with ending: "<<outputfileending<<endl;
  for(unsigned long c1=0;c1<inputfileroot.size();c1++){
    os<<"#ChomoCombine: fileroot["<<c1<<"]="<<inputfileroot[c1]<<endl;
  }
  for(unsigned long c1=0;c1<oldforcefile.size();c1++) {
    if(oldforcefile[c1].substr(0,8).compare(string("#Cfactor"))==0) continue;
    if(oldforcefile[c1].substr(0,14).compare(string("#ChomoCombine:"))==0) continue;
    os<<oldforcefile[c1]<<endl;
  }
  fb.close();
  return(true);
}

////////////////////////////////////////////
// READING DATA

std::string ChromoCombineWorker::removequotes(std::string str){
	size_t found=str.find('"');
	while(found!=std::string::npos) {
		str.erase(found,1);
		found=str.find('"');
	}
	return(str);
    }

long ChromoCombineWorker::getRowIndex(string dval){// get the column index for the rows
  //////////////////////////////
  // **** THIS IS A VERY SLOW WAY TO DO IT!!!
  // FIXME

  for(unsigned long c1=0;c1<rownames.size();c1++) if(rownames[c1].compare(dval)==0) {
    indexcounts[c1]++;
    return(c1);
  }
  indexcounts.push_back(1);
  rownames.push_back(dval);
  return(rownames.size()-1);
}

long ChromoCombineWorker::getColIndex(string dval){// get the column index for the rows
  //////////////////////////////
  // **** THIS IS A VERY SLOW WAY TO DO IT!!!
  // FIXME

  for(unsigned long c1=0;c1<maincolnames.size();c1++) if(maincolnames[c1].compare(dval)==0) {
    return(c1);
  }
  cerr<<"ERROR: Column name "<<dval<<" not found in IDs. You have either specified an incomplete list of populations with \"-i <idfile>\" or you have used chromopainter in donor mode, without specifying the complete set of population names via \"-i <idfile>\"."<<endl;
  throw(runtime_error("getColIndex name matching error"));
  return(-1);
}

std::istream& ChromoCombineWorker::safeGetline(std::istream& is, std::string& t)
{
    string myline;
    if ( getline( is, myline ) ) {
       if ( myline.size() && myline[myline.size()-1] == '\r' ) {
           t = myline.substr( 0, myline.size() - 1 );
       }
       else {
           t = myline;
       }
    }else{
      t=myline;
    }
    return is;
}

bool ChromoCombineWorker::readIds(string filename) {
  ifstream file;
  size_t found;
  string line;
  file.open(filename.data());//Open file
  if(verbose>0) cout<<"Opening ID file "<<filename<<endl;
  if(!file.good()) return(false);
  maincolnames.clear();
  while (1)
    {
      safeGetline(file,line);//Read next line from file
      found=line.find_first_of(",\t ");
      string dval=line.substr(0,found);
      if(dval.size()>0)maincolnames.push_back(dval);
      found=string::npos;
      if (file.eof())
	break;//Stop if end of file
    }
  if(verbose>0) cout<<"Read "<<maincolnames.size()<<" population IDs from file."<<endl;
  if(maincolnames.size()>0) return(true);
  return(false);
}

vector<long> ChromoCombineWorker::setColumnIndices(long filetype){
  // Sets thecolumn reordering index vector based on the predefined maincolnames. This is different to the table colnames if using an id file
  vector<long> tindexes;
  long toff=0; // offset 
  if(filetype==1 || filetype==2){ // Account for the num.regions column
    tindexes.push_back(toff++);
  }
  
  if(colnames.size() > maincolnames.size()+toff) {
    cerr<<"ERROR: Data has too many columns ("<<colnames.size()<<"), compared to IDs ("<<maincolnames.size()<<")"<<endl;
    throw(runtime_error("Too many columns compared to IDs!"));
  }

  for(unsigned int c1=0;c1<colnames.size();++c1){
    //    cout<<"Read "<<colnames[c1]<<flush;
    if(colnames.size()==maincolnames.size()){
      if(maincolnames[c1]==colnames[c1]) {
	tindexes.push_back(c1+toff);
	//	cout<<" quickly as "<<tindexes.back()<<endl;
	continue;
      }
    }
    tindexes.push_back(getColIndex(colnames[c1])+toff);
    //    cout<<" as "<<tindexes.back()<<endl;
  }
  return(tindexes);
}

//
bool ChromoCombineWorker::readTable(string filename,vector< vector <double> > *dest,vector<long> *indexes,long filetype){
    ifstream file;
    size_t found;
    vector<long> colindexes;/// possible reordering of the columns based on predefined column names
    file.open(filename.data());//Open file
    if(verbose>0) cout<<"Opening file "<<filename<<endl;
    if(!file.good()) return(false);
    dest->clear();
    
    long lineon=0;
    bool usingrownames=true;
    bool applynames=false;
    //    if(colnames.size()==0){
    if(filetype==0){
      applynames=true;
      colnames.clear();
    }
    long ignore=0;// allow ignoring lines later
    string line;
    while (1) {
      //cout<<"Reading line "<<lineon<<", ignoring "<<ignore<<endl;
	safeGetline(file,line);//Read next line from file
	//cout<<"Read:"<<line<<": length "<<line.size()<<endl;
	if (file.eof()) break;//Stop if end of file
	if(line.size()==0) continue; // ignore empty lines
        if(lineon<ignore) {lineon++;continue;//ignore the first ignore lines
	}else{
	   if(lineon-ignore==0) {// read the column names
	      long colon=0;
		  while (1)
		  {
		    found=line.find_first_of(",\t ");
// check for cfactor
		    string dval=line.substr(0,found);
		    if(dval.compare("#Cfactor")==0 && colon==0){// its a c definition line
		      ignore++;
		      line.clear();
		      found=string::npos;
		      break;
		    }// column headers
		    else if(applynames) {
		      string dval=line.substr(0,found);
		      dval=removequotes(dval);
		      if(colon==0) {topleftentry=dval;
		      }else colnames.push_back(dval.c_str());
		    }else {line.clear(); break;}
		    if(found==string::npos) line.clear();
		    else line=line.substr(found+1,line.length());
		    if(line.length()==0) break;
		    colon++;
		  }
		  if(maincolnames.size()==0) maincolnames=colnames; // set up the columns
	      lineon++;continue;
	   }
	  if (line.size()==0 || line[0]=='#' )
	      continue;//Ignore empty lines or comments lines
	  else if(lineon-ignore>0)
	  {//read the row of data
	    if(colindexes.size()==0) colindexes=setColumnIndices(filetype);
	      long colon=0;
	      dest->push_back(vector<double> (maincolnames.size() + (filetype==1 || filetype==2)));
	      //cout<<"Dest size="<<dest->back().size()<<endl;
	      
		  while (1)
		  {//read each element of the matrix
		    found=line.find_first_of(",\t ");
		    string dval=line.substr(0,found);//this contains the element
		    dval=removequotes(dval);// get rid of any quotes around names
		    //cout<<"Column "<<colon<<" = "<<dval<<endl;
		    //cout<<"colindexes["<<colon<<"]="<<colindexes[colon]<<endl;
		    if(found==string::npos) line.clear();
		    else line=line.substr(found+1,line.length());
		    if(colon==0) {
		      if(colnames.size()>0) {
			indexes->push_back(getRowIndex(dval));
		      }else {
			// no col or row names; matrices are required to be the same order
			cerr<<"WARNING: No row and column names!"<<endl;
			usingrownames=false;
			indexes->push_back(indexes->size());
			dest->back()[colindexes[colon]]=atof(dval.c_str());
		      }
		      colon++;
		    }else {// colon>0
		      //cout<<"colon="<<colon<<" urn="<<usingrownames<<endl;
		      dest->back()[colindexes[colon-usingrownames]]=atof(dval.c_str());
		      if(line.length()==0) break;
		      colon++;
		    }
		  }
	      lineon++;
	      continue;
	  }
	}
    }
    file.close();
    if(numeffinds<(long)rownames.size()) numeffinds=(long)rownames.size();
    if(dest->size()==0){
      cerr<<"File "<<filename<<" has no rows.  Did ChromoPainter complete successfully?"<<endl; 
      return(false);
    }
    return(true);
}

bool ChromoCombineWorker::addToChunkCounts(vector< vector <double> > *tmpptr,vector<long> *indexes)
{
  unsigned long lengthfirst=0;

  if(indexes->size()!=tmpptr->size()) throw(runtime_error("Internal error adding to chunks!"));
  for(unsigned long c1=0;c1<tmpptr->size();c1++){
    if((long)counts.size()<=indexes->at(c1)) counts.resize(indexes->at(c1)+1,vector<double>(0));

    if(c1==0) {lengthfirst=tmpptr->at(0).size();
    }else if(tmpptr->at(c1).size()!=lengthfirst){ cerr<<"Error: ChunkCount rows not of equal size! (Expected "<<lengthfirst<<", received "<<tmpptr->at(indexes->at(c1)).size()<<")"<<endl;return(false);}
    if(counts[indexes->at(c1)].size()<tmpptr->at(c1).size()) {
      counts[indexes->at(c1)].resize(tmpptr->at(c1).size(),0);
    }
    for(unsigned long c2=0;c2<tmpptr->at(c1).size();c2++){
      counts[indexes->at(c1)][c2]+=tmpptr->at(c1)[c2];
    }
  }
  return(true);
}

bool ChromoCombineWorker::addToRegionChunkCounts(vector< vector <double> > *tmpptr,vector<long> *indexes)
{
  unsigned long lengthfirst=0;
  for(unsigned long c1=0;c1<tmpptr->size();c1++){
    if((long)regioncounts.size()<=indexes->at(c1)) regioncounts.resize(indexes->at(c1)+1,vector<double>(0));
    if(c1==0) {lengthfirst=tmpptr->at(0).size();
    }else if(tmpptr->at(c1).size()!=lengthfirst){ cerr<<"Error: RegionChunkCount rows not of equal size!"<<endl;return(false);}

    if(regioncounts[indexes->at(c1)].size()<tmpptr->at(c1).size()) regioncounts[indexes->at(c1)].resize(tmpptr->at(c1).size(),0);
    for(unsigned long c2=0;c2<tmpptr->at(c1).size();c2++){
      regioncounts[indexes->at(c1)][c2]+=tmpptr->at(c1)[c2];
    }
  }
  return(true);
}

bool ChromoCombineWorker::addToRegionSquaredChunkCounts(vector< vector <double> > *tmpptr,vector<long> *indexes)
{
  unsigned long lengthfirst=0;
  for(unsigned long c1=0;c1<tmpptr->size();c1++){
    if((long)regioncountssq.size()<=indexes->at(c1)) regioncountssq.resize(indexes->at(c1)+1,vector<double>(0));
    if(c1==0) {lengthfirst=tmpptr->at(0).size();
    }else if(tmpptr->at(c1).size()!=lengthfirst){ cerr<<"Error: RegionChunkCountSq rows not of equal size!"<<endl;return(false);}

    if(regioncountssq[indexes->at(c1)].size()<tmpptr->at(c1).size()) regioncountssq[indexes->at(c1)].resize(tmpptr->at(c1).size(),0);
    for(unsigned long c2=0;c2<tmpptr->at(c1).size();c2++){
      regioncountssq[indexes->at(c1)][c2]+=tmpptr->at(c1)[c2];
    }
  }
  return(true);
}

bool ChromoCombineWorker::addToChunkLengths(vector< vector <double> > *tmpptr,vector<long> *indexes)
{
  unsigned long lengthfirst=0;
  for(unsigned long c1=0;c1<tmpptr->size();c1++){
    if((long)lengths.size()<=indexes->at(c1)){
      lengths.resize(indexes->at(c1)+1,vector<double>(0));
    }
    if(c1==0) {lengthfirst=tmpptr->at(0).size();
    }else if(tmpptr->at(c1).size()!=lengthfirst){ cerr<<"Error: ChunkLength rows not of equal size!"<<endl;return(false);}

    if(lengths[indexes->at(c1)].size()<tmpptr->at(c1).size()){
      lengths[indexes->at(c1)].resize(tmpptr->at(c1).size(),0);
    }
    for(unsigned long c2=0;c2<tmpptr->at(c1).size();c2++){
      lengths[indexes->at(c1)][c2]+=tmpptr->at(c1)[c2];
    }
  }
  return(true);
}

bool ChromoCombineWorker::addToMutationProbs(vector< vector <double> > *tmpptr,vector<long> *indexes)
{
  unsigned long lengthfirst=0;
  for(unsigned long c1=0;c1<tmpptr->size();c1++){
    if((long)muts.size()<=indexes->at(c1)) muts.resize(indexes->at(c1)+1,vector<double>(0));
    if(c1==0) {lengthfirst=tmpptr->at(c1).size();
    }else if(tmpptr->at(c1).size()!=lengthfirst){ cerr<<"Error: MutationCounts rows not of equal size!"<<endl;return(false);}

    if(muts[indexes->at(c1)].size()<tmpptr->at(c1).size()) muts[indexes->at(c1)].resize(tmpptr->at(c1).size(),0);
    for(unsigned long c2=0;c2<tmpptr->at(c1).size();c2++){
      muts[indexes->at(c1)][c2]+=tmpptr->at(c1)[c2];
    }
  }
  return(true);
}

void ChromoCombineWorker::makeRegionChunkCounts(vector< vector <double> > *tmpptr,vector< vector <double> > *chunkcounts)
{
  tmpptr->clear();
  vector<double> rowsums(chunkcounts->size(),0.0);
  double allsums=0;
  double allsumssq=0;
  double ncount=0;
  for(unsigned long c1=0;c1<chunkcounts->size();c1++) {
    for(unsigned long c2=0;c2<chunkcounts->at(c1).size();c2++) {
      rowsums[c1]+=chunkcounts->at(c1)[c2]/(double) chunkcounts->at(c1).size();
      allsums+=chunkcounts->at(c1)[c2];
      allsumssq+=chunkcounts->at(c1)[c2] * chunkcounts->at(c1)[c2];
      ncount+=1.0;
    }
  }
  allsums=allsums/ncount;
  allsumssq=allsumssq/ncount;
  completegenomesnormfactor+=allsums;
  completegenomesnormfactorsq+=allsumssq;
  completegenomesnormfactorcount+=1.0;
  for(unsigned long c1=0;c1<chunkcounts->size();c1++) {
    tmpptr->push_back(vector<double>(chunkcounts->at(c1).size()+1,0));
    tmpptr->at(c1)[0]=1;
    for(unsigned long c2=0;c2<chunkcounts->at(c1).size();c2++) {
	if(rowsums[c1]==0) {
		tmpptr->at(c1)[c2+1]=0;
      	}else {
 	     tmpptr->at(c1)[c2+1]=chunkcounts->at(c1)[c2] / rowsums[c1];
	}
    }
  }
}

void ChromoCombineWorker::makeRegionSquaredChunkCounts(vector< vector <double> > *tmpptr,vector< vector <double> > *chunkcounts)
{
  tmpptr->clear();
  vector<double> rowsums(chunkcounts->size(),0.0);
  for(unsigned long c1=0;c1<chunkcounts->size();c1++) {
    for(unsigned long c2=0;c2<chunkcounts->at(c1).size();c2++) {
      rowsums[c1]+=chunkcounts->at(c1)[c2]/(double) chunkcounts->at(c1).size();
    }
  }
  for(unsigned long c1=0;c1<chunkcounts->size();c1++) {
    tmpptr->push_back(vector<double>(chunkcounts->at(c1).size()+1,0));
    tmpptr->at(c1)[0]=1;
    for(unsigned long c2=0;c2<chunkcounts->at(c1).size();c2++) {
	if(rowsums[c1]==0) {
		tmpptr->at(c1)[c2+1]=0;
      	}else {
	      tmpptr->at(c1)[c2+1]=chunkcounts->at(c1)[c2]*chunkcounts->at(c1)[c2]/ rowsums[c1]/ rowsums[c1];
	}
    }
  }
}

long ChromoCombineWorker::readFilesOneRoot(string root){
  vector<string> tlist=requiredendings;
  vector<string> allf=getFileList(root,tlist);
  vector< vector <double> > tmp;
  vector< vector <double> > tmpchunkcounts;
  vector<long> rownums;
  vector< vector <double> > *tmpptr=&tmp;
  for(unsigned long fon=0;fon<allf.size();fon++){
    long curtest=-1;
    for(unsigned long fon2=0;fon2<possibleendings.size();fon2++){
      if(requiredendings[fon].compare(possibleendings[fon2])==0){ curtest=fon2; fon2=possibleendings.size();}
    }
    if(curtest<0){
      cerr<<"Error reading files: File endings are confused?"<<endl;return(false);
    }
    int tabletype=curtest;
    if(fon==0 && !wantcounts) tabletype=0;
    if(!readTable(allf[fon],tmpptr,&rownums,tabletype)){
      cerr<<"Error reading files: file "<<allf[fon]<<" not read correctly.  Does it exist, and have the correct number of rows?"<<endl;return(0);
    }
    if(completegenomes) {
      if(curtest==0) tmpchunkcounts=tmp;
      else if(curtest==1) makeRegionChunkCounts(tmpptr,&tmpchunkcounts);
      else if(curtest==2) makeRegionSquaredChunkCounts(tmpptr,&tmpchunkcounts);
    }
    switch(curtest){
      case 0:	if(!addToChunkCounts(tmpptr,&rownums)) return(-1);break;
      case 1:	if(!addToRegionChunkCounts(tmpptr,&rownums)) return(-2);break;
      case 2:	if(!addToRegionSquaredChunkCounts(tmpptr,&rownums)) return(-3);break;
      case 3:	if(!addToChunkLengths(tmpptr,&rownums))return(-4);break;
      case 4:	if(!addToMutationProbs(tmpptr,&rownums))return(-5);break;
      default: cerr<<"Error reading files: invalid file ending found?"<<endl;return(-6);
    }
  }
  return(1);
}

long ChromoCombineWorker::readFiles(){
  long ret=1;
  for(unsigned long c1=0;c1<inputfileroot.size();c1++) {
    if(verbose>0) cout<<"Reading fileroot "<<c1+1<<" of "<<inputfileroot.size()<<" called "<<inputfileroot[c1]<<endl;
    ret=readFilesOneRoot(inputfileroot[c1]);
    if(ret<=0) return(ret);
  }
  return(ret);
}

long ChromoCombineWorker::validateMatrices(){
  if(wantcounts) {
    if(counts.size()==0) {
      cerr<<"Counts matrix is empty!"<<endl;
      return(-1);
    }
    if(counts.size() != regioncounts.size()) {
      cerr<<"Region counts has "<<regioncounts.size()<<" rows and counts has "<<counts.size()<<" rows!"<<endl;
      return(-2);
    }
    if(counts.size() != regioncountssq.size()) {
      cerr<<"Region counts Squared has "<<regioncountssq.size()<<" rows and counts has "<<counts.size()<<" rows!"<<endl;
      return(-3);
    }
  }
  if(wantlengths && wantcounts) if(counts.size() != lengths.size()) {
    cerr<<"Lengths has "<<lengths.size()<<" rows and counts has "<<counts.size()<<" rows!"<<endl;
    return(-4);
  }
  if(wantmuts && wantcounts) if(counts.size() != muts.size()) {
    cerr<<"Mutations has "<<muts.size()<<" rows and counts has "<<counts.size()<<" rows!"<<endl;
    return(-5);
  }
  for(unsigned long c1=0;c1<counts.size();c1++) {
    if(wantcounts) {
      if(counts.size()==0) {
	cerr<<"Counts matrix is empty!"<<endl;
	return(-1);
      }
      if(counts.size() != regioncounts.size()) {
	cerr<<"Region counts has "<<regioncounts.size()<<" rows and counts has "<<counts.size()<<" rows!"<<endl;
	return(-2);
      }
      if(counts.size() != regioncountssq.size()) {
	cerr<<"Region counts Squared has "<<regioncountssq.size()<<" rows and counts has "<<counts.size()<<" rows!"<<endl;
	if(counts[c1].size()==0) {
	  cerr<<"Counts matrix row "<<c1<<" is empty!"<<endl;
	  return(-1);
	}
      }
      if(1+counts[c1].size() != regioncounts[c1].size()) {
	cerr<<"Region counts row "<<c1<<" has "<<((long)regioncounts[c1].size())-1<<" columns and counts has "<<counts[c1].size()<<" columns!"<<endl;
	return(-2);
      }
      if(1+counts[c1].size() != regioncountssq[c1].size()) {
	cerr<<"Region counts Squared row "<<c1<<" has "<<((long)regioncountssq[c1].size())-1<<" columns and counts has "<<counts[c1].size()<<" columns!"<<endl;
	return(-3);
      }
    }
    if(wantlengths && wantcounts) if(counts[c1].size() != lengths[c1].size()) {
      cerr<<"Lengths row "<<c1<<" has "<<lengths[c1].size()<<" columns and counts has "<<counts[c1].size()<<" columns!"<<endl;
      return(-4);
    }
    if(wantmuts && wantcounts) if(counts[c1].size() != muts[c1].size()) {
      cerr<<"Mutations row "<<c1<<" has "<<muts[c1].size()<<" columns and counts has "<<counts[c1].size()<<" columns!"<<endl;
      return(-5);
    }
  }
  long firstcount=0;
  for(unsigned long c1=0;c1<indexcounts.size();c1++){
    if(c1==0)firstcount=indexcounts[c1];
    else if(firstcount!=indexcounts[c1]){
      cerr<<"WARNING: Seen individual 0 called "<<rownames[0]<<" a total of "<<firstcount<<" times, but individual "<<c1<<" called "<<rownames[c1]<<" a total of "<<indexcounts[c1]<<" times!"<<endl;
      return(1);
    }
  }
  return(0);
}

////////////////////////////////////////////
// PROCESSING DATA

void  ChromoCombineWorker::calcTheoreticalVar(){
  theoreticalvar=counts;
  for(unsigned long c1=0;c1<theoreticalvar.size();c1++) {
    double rowsum=0;
    for(unsigned long c2=0;c2<theoreticalvar[c1].size();c2++) {
      if(regionsums[c1]>1) rowsum+=regioncounts[c1][c2+1];
    }
    for(unsigned long c2=0;c2<theoreticalvar[c1].size();c2++) {
      if(rowsum==0 || regionsums[c1]==0) theoreticalvar[c1][c2]=0;
      else {
	double pest=regioncounts[c1][c2+1]/rowsum;
	theoreticalvar[c1][c2] =rowsum * pest * (1.0 - pest)/regionsums[c1];
      }
    }
  }
}

void ChromoCombineWorker::calcEmpiricalVar(){
  empiricalvar=counts;
  for(unsigned long c1=0;c1<empiricalvar.size();c1++) for(unsigned long c2=0;c2<empiricalvar[c1].size();c2++) {
    if(regionsums[c1]<=1) {
      empiricalvar[c1][c2]=0;
    }else{
      empiricalvar[c1][c2]=regioncountssq[c1][c2+1]/(regionsums[c1]-1.0) - regioncounts[c1][c2+1]*regioncounts[c1][c2+1]/regionsums[c1]/(regionsums[c1]-1.0);
    }
  }
}

void ChromoCombineWorker::calcTheoreticalVarForce(){
  theoreticalvar=vector<vector<double> >(numeffinds,vector<double>(numeffinds,0));
  for(unsigned long c1=0;c1<theoreticalvar.size();c1++) {
    double rowsum=0;
    for(unsigned long c2=0;c2<theoreticalvar[c1].size();c2++) {
      if(forcedregionsums[c1]>1) rowsum+=forcedregioncounts[c1][c2];
    }
    for(unsigned long c2=0;c2<theoreticalvar[c1].size();c2++) {
      if(rowsum==0 || forcedregionsums[c2]==0) theoreticalvar[c1][c2]=0;
      else {
	double pest=forcedregioncounts[c1][c2]/rowsum;
	theoreticalvar[c1][c2] =rowsum * pest * (1.0 - pest)/forcedregionsums[c1];
      }
    }
  }
}

void ChromoCombineWorker::calcEmpiricalVarForce(){
  empiricalvar=vector<vector<double> >(numeffinds,vector<double>(numeffinds,0));
  for(unsigned long c1=0;c1<empiricalvar.size();c1++) for(unsigned long c2=0;c2<empiricalvar[c1].size();c2++) {
    if(forcedregionsums[c1]<=1) {
      empiricalvar[c1][c2]=0;
    }else{
      empiricalvar[c1][c2]=forcedregioncountssq[c1][c2]/(forcedregionsums[c1]-1.0) - forcedregioncountssumsq[c1][c2]/forcedregionsums[c1]/(forcedregionsums[c1]-1.0);
    }
  }
}

double ChromoCombineWorker::calcC(bool useforce){
  if(useforce){
    calcEmpiricalVarForce();
    calcTheoreticalVarForce();
  }else{
    calcEmpiricalVar();
    calcTheoreticalVar();
  }
  double tsum=0;
  double tdenom=0;
  for(unsigned long c1=0;c1<empiricalvar.size();c1++) for(unsigned long c2=0;c2<empiricalvar[c1].size();c2++) {
      double val=0;

    if(empiricalvar.size() ==empiricalvar[c1].size()){// use symmetrisation
	if(theoreticalvar[c1][c2]+theoreticalvar[c2][c1]>1e-5) {
	  val=2.0 * (empiricalvar[c1][c2]+empiricalvar[c2][c1])/(theoreticalvar[c1][c2]+theoreticalvar[c2][c1]);
	  tsum+=val;
	  tdenom+=1.0;
	}
    }else{// use raw
	if(theoreticalvar[c1][c2]>1e-5) {
	  val=2.0 * (empiricalvar[c1][c2])/(theoreticalvar[c1][c2]);
	  tsum+=val;
	  tdenom+=1.0;
	}
    }// end if ,, else
  }// end for
  if(tdenom==0 || tsum==0) {
    if(verbose>=0)cerr<<"WARNING: No regions found, insufficient data for calculating c.  Try running chromopainter either with a smaller \"-k\" option, or run chromocombine with the \"-C\" option. See http://www.paintmychromosomes.com (faq page) for a discussion of this issue."<<endl;
    return(0);
  }
  return(tsum/tdenom);
}

void ChromoCombineWorker::finaliseData(){
  sortRows();
  if(!wantcounts) return;
  regionsums=vector<double>(regioncounts.size(),0);
  for(unsigned long c1=0;c1<regionsums.size();c1++) regionsums[c1]=regioncounts[c1][0];
  if(completegenomes) {/// correct the total number of chunks
    for(unsigned long c1=0;c1<regioncounts.size();c1++){
      for(unsigned long c2=1;c2<regioncounts[c1].size();c2++){
	regioncounts[c1][c2]*=completegenomesnormfactor/completegenomesnormfactorcount;
	regioncountssq[c1][c2]*=completegenomesnormfactor*completegenomesnormfactor/completegenomesnormfactorcount/completegenomesnormfactorcount;
      }
    }
  }
  return;
}

vector<long> ChromoCombineWorker::getRowOrder(){
  vector<long> torder,empty;
  if(rownames.size()!=maincolnames.size()) return(empty);/// don't reorder if there isn't a perfect match
  torder=vector<long> (maincolnames.size(),-1);
  for(unsigned long c1=0;c1<rownames.size();c1++) {
    bool found=false;
    for(unsigned long c2=0;c2<maincolnames.size();c2++) {
      if(maincolnames[c2].compare(rownames[c1])==0){
	torder[c2]=c1;
	found=true;
	c2=maincolnames.size();
      }
    }
    if(!found) {
	cerr<<"When matching rows and columns, I couldn't find row name "<<rownames[c1]<<" in the column names"<<endl;
      return(empty);// names don't match so don't sort
    }
  }
  for(unsigned long c1=0;c1<torder.size();c1++) {
	if(torder[c1]<0) {
		cerr<<"When matching rows and columns, I couldn't find column name "<<maincolnames[c1]<<" in the row names"<<endl;
		return(empty);
	}
  }
  return(torder);
}

void ChromoCombineWorker::reorderMatrix(vector<vector<double> > *mat,vector<long> *torder){
  vector<vector<double> > copymat=(*mat);
  for(unsigned long c1=0;c1<mat->size();c1++){
    mat->at(c1) = copymat[torder->at(c1)];
  }
}

void ChromoCombineWorker::reorderVector(vector<long > *mat,vector<long> *torder){
  vector<long> copymat=(*mat);
  for(unsigned long c1=0;c1<mat->size();c1++){
    mat->at(c1) = copymat[torder->at(c1)];
  }
}

void ChromoCombineWorker::reorderVector(vector<string > *mat,vector<long> *torder){
  vector<string> copymat=(*mat);
  for(unsigned long c1=0;c1<mat->size();c1++){
    mat->at(c1) = copymat[torder->at(c1)];
  }
}


void ChromoCombineWorker::sortRows(){
  vector<long> torder=getRowOrder();
  if(torder.size()==0) {
    if(verbose>=0){
    cout<<"Not sorting due to name differences between rows and columns."<<endl;
    cout<<"This could mean that your data has been read incorrectly!"<<endl;
	cout<<"Some suggestions are:"<<endl;
    cout<<"1. Did you split the data into many files? Check that all the files you expected to be created by chromopainter exist, and have the correct number of lines. You may wish to try rerunning stating the input files explicitly, in case autodetection of files has failed."<<endl;
    cout<<"2. Are there file line-ending issues? (did you run chromocombine"<<endl<<" on a Windows/UNIX/Mac machine and chromopainter elsewhere?)"<<endl;
    cout<<"3. Alternatively, you may be using chromopainter in donor population mode, in which case this warning is expected."<<endl;
    cout<<"Continuing, in case you intended this. Expect problems!"<<endl;
    }
  }else{
	  reorderMatrix(&counts,&torder);
	  reorderMatrix(&regioncounts,&torder);
	  reorderMatrix(&regioncountssq,&torder);
	  reorderMatrix(&lengths,&torder);
	  reorderMatrix(&muts,&torder);
	  reorderVector(&indexcounts,&torder);
	  reorderVector(&rownames,&torder);
  }
}

////////////////////////////////////////////
// WRITING OUTPUT

bool ChromoCombineWorker::writeRowName(std::ostream * os,long rownum,long filetype){
  if(rownum<(long)rownames.size()) *os<<rownames[rownum]<<" ";
  else *os<<rownum<<" ";
  return(true);
}

bool ChromoCombineWorker::writeHeaderRow(std::ostream * os,long filetype){
  *os<<topleftentry<<" ";
  if(filetype==1 || filetype==2) *os<<"num.regions ";
  if(maincolnames.size()>0){
    for(unsigned long c1=0;c1<maincolnames.size()-1;c1++) *os<<maincolnames[c1]<<" ";
    *os<<maincolnames[maincolnames.size()-1];
  }
  *os<<endl;
  return(true);
}

void ChromoCombineWorker::writematrix(std::ostream * os,vector<vector<double> > *mat,long filetype){
  writeHeaderRow(os,filetype);
  for(long c1=0;c1<(long)mat->size();c1++){
    writeRowName(os,c1,filetype);
    for(long c2=0;c2<((long)mat->at(c1).size())-1;c2++) {
      *os<<mat->at(c1)[c2]<<" ";
    }
    if(mat->at(c1).size()>0) {
      *os<<mat->at(c1)[mat->at(c1).size()-1]<<endl;
    }else *os<<endl;
  }

}

bool ChromoCombineWorker::writeChunkCountFile(){
  string filename=outputfileroot;
  filename.append(possiblemiddles[0]).append(outputfileoutending);
  filebuf fb;
  try{
    fb.open (filename.c_str(),ios::out);
  }catch(exception &x){
    cerr<<"Error opening file!"<<endl<<x.what()<<endl; return false;}

  ostream os (&fb);
  os<<setprecision(9);

  os<<"#Cfactor "<<cval<<endl;
  writematrix(&os,&counts,0);

  fb.close();
  return(true);
}

bool ChromoCombineWorker::writeRegionChunkCountFile(){
  string filename=outputfileroot;
  filename.append(possiblemiddles[1]).append(outputfileoutending);
  filebuf fb;
  try{
    fb.open (filename.c_str(),ios::out);
    }catch(exception &x){
    cerr<<"Error opening file!"<<endl<<x.what()<<endl; return false;}
  ostream os (&fb);
  os<<setprecision(9);
  writematrix(&os,&regioncounts,1);
  fb.close();
  return(true);
}

bool ChromoCombineWorker::writeRegionChunkCountSqFile(){
  string filename=outputfileroot;
  filename.append(possiblemiddles[2]).append(outputfileoutending);
  filebuf fb;
  try{
    fb.open (filename.c_str(),ios::out);
    }catch(exception &x){
    cerr<<"Error opening file!"<<endl<<x.what()<<endl; return false;}
  ostream os (&fb);
  os<<setprecision(9);
  writematrix(&os,&regioncountssq,2);
  fb.close();
  return(true);
}

bool ChromoCombineWorker::writeChunkLengthsFile(){
  string filename=outputfileroot;
  filename.append(possiblemiddles[3]).append(outputfileoutending);
  filebuf fb;
  try{
    fb.open (filename.c_str(),ios::out);
    }catch(exception &x){
    cerr<<"Error opening file!"<<endl<<x.what()<<endl; return false;}
  ostream os (&fb);
  os<<setprecision(9);
  writematrix(&os,&lengths,3);
  fb.close();
  return(true);
}

bool ChromoCombineWorker::writeMutationProbFile(){
  string filename=outputfileroot;
  filename.append(possiblemiddles[4]).append(outputfileoutending);
  filebuf fb;
  try{
    fb.open (filename.c_str(),ios::out);
    }catch(exception &x){
    cerr<<"Error opening file!"<<endl<<x.what()<<endl; return false;}
  ostream os (&fb);
  os<<setprecision(9);
  writematrix(&os,&muts,4);
  fb.close();
  return(true);
}


bool ChromoCombineWorker::writeFiles(){
  if(verbose>0) cout<<"Writing files"<<endl;
  if(wantcounts) {
    if(verbose>0) cout<<"Writing chunkcounts file"<<endl;
    if(!writeChunkCountFile()){
      cerr<<"Error writing ChunkCount file!"<<endl;
      return(false);
    }
    if(verbose>0) cout<<"Writing regionchunkcounts file"<<endl;
    if(!writeRegionChunkCountFile()){
      cerr<<"Error writing RegionChunkCount file!"<<endl;
      return(false);
    }
    if(verbose>0) cout<<"Writing regionsquaredchunkcounts file"<<endl;
    if(!writeRegionChunkCountSqFile()){
      cerr<<"Error writing RegionChunkCountsq file!"<<endl;
      return(false);
    }
  }
  if(wantlengths)if(!writeChunkLengthsFile()){
	if(verbose>0) cout<<"Writing chunklengths file"<<endl;
    cerr<<"Error writing ChunkLengths file!"<<endl;
    return(false);
  }
  if(wantmuts) if(!writeMutationProbFile()){
    if(verbose>0) cout<<"Writing mutationprob file"<<endl;
    cerr<<"Error writing MutationProb file!"<<endl;
    return(false);
  }
  return(true);
}
