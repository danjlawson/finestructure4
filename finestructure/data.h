#ifndef DATA_H
#define DATA_H
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <fstream>
#include <algorithm>
#include <string.h>

#include "safegetline.h" // for cross-platform safely getting lines

using namespace std;

namespace fines
{

/**
    @brief snip data converted to number of blocks copied by individuals (i.e. a matrix)
*/
class Data
{
public:
  Data(std::string filename,int ignore=0,bool yhead=false,bool xhead=false,bool square=true);///<Reads in the data from a file
  Data(std::vector < std::vector<double> > * inmat,std::vector <std::string> namesin,std::vector <std::string> colnamesin=std::vector<std::string>());///<Reads in the data from a matrix
    Data(Data * datin);///<Copies a data
    ~Data();
    void print(std::ostream * out);///< Writes data to a stream
    void output(std::string filename);///< Writes the data to a file
    inline double get(int i,int j,bool cf=false)
    {
      if(cf) return data[i][j]/datacfactor;
      return data[i][j];
    };///<Get accessor to the data
    inline double rowsum(int i,bool cf=false)
    {
      if(cf) return(datarowsums[i]/datacfactor);
        return datarowsums[i];
    };///<Get accessor to the data
    inline double colsum(int i,bool cf=false)
    {
      if(cf) return(datacolsums[i]/datacfactor);
        return datacolsums[i];
    };///<Get accessor to the data
    inline double allsum(bool cf=false)
    {
      if(cf) return(sumXall/datacfactor);
      return(sumXall);
    }///<Get accessor to the data
    inline std::vector<double>* get(int i)
    {
        return(&data[i]);
    };///<Get accessor to the data
    inline void set(int i,int j,double c)
    {
	sumXall+=c-data[i][j];
	datarowsums[i]+=c-data[i][j];
	datacolsums[j]+=c-data[i][j];
        data[i][j]=c;
    };///<Set accessor to the data
    inline int getN()
    {
        return n;
    };///<Returns the number of individuals
    inline int getDim()
    {
        return data.size();
    };///<Returns the number of used individual rows (i.e. merging super individuals)
    inline int getNCols()
    {
      if(square) return(getDim());
      if(data.size()==0) return(0); // if no rows, we don't report any columns
      return data[0].size();
    };///<Returns the number of used individual columns, if not square
    inline std::string getnames(int i)
    {
	if((int)names.size()>i && i>=0) {
	  return(names.at(i));
	}else{
	  std::stringstream out;
	  out << i;
	  return(out.str());
	}
    }///<Returns the name of the individual at i
    inline std::string getcolnames(int i)
    {
	if(square && (int)names.size()>i && i>=0) {
	  return(names.at(i));
	}else if(!square && (int)colnames.size()>i && i>=0) {
	  return(colnames.at(i));
	}else{
	  std::stringstream out;
	  out << "IND" << i;
	  return(out.str());
	}
    }///<Returns the name of the individual at i
    inline void forceNames(){
      while(names.size()<data.size()){
	  std::stringstream out;
	  out << "IND"<< names.size();
	names.push_back(out.str());
      }
      if(!square){
	if(data.size()>0) while(colnames.size()<data[0].size()){
	  std::stringstream out;
	  out << "IND"<< colnames.size();
	  colnames.push_back(out.str());
	}
      }
    }///< Creates the names IND1..N. Not necessary since we create on the fly.
    
    inline int nindiv(int i){
      if((int)neff.size()>i && i>=0) return(neff[i]);
      throw(string("Invalid individual count requested!")); return(-1);
    }///< Number of individuals in a population (counting the numer inside super individuals)
    
    inline int getIndex(std::string name){
	for(unsigned int c1=0;c1<names.size();c1++) if(name.compare(names[c1])==0) return(c1);
	return(-1);
    }///<Returns the index of a named individual (or -1 if not found)
    inline int getColIndex(std::string name){
      if(square) return(getIndex(name));
      for(unsigned int c1=0;c1<colnames.size();c1++) if(name.compare(colnames[c1])==0) return(c1);
      return(-1);
    }///<Returns the index of a named individual (or -1 if not found)
    inline std::vector<std::vector<double> > * getMatrix(){
	return(&data);
    }///< returns a pointer to the raw data matrix
    inline bool testData(int rowsexpected=-1){
    	if(rowsexpected<0) rowsexpected=n;
	int colsexpected=rowsexpected;
	if(!square) colsexpected = k;
    	for (unsigned int i=0;i<data.size();i++) {
        	if (colsexpected!=(int)data[i].size())
        	{
            		cerr<<"Data "<<i<<" is inconsistent: "<<data[i].size()<<"!="<<rowsexpected<<endl;
	  		cout<<"i="<<i<<" data="<<data[i][0]<<","<<data[i][1]<<","<<data[i][2]<<endl;
			std::string tfilename=filename;
			tfilename.append(".testdata");
			cout<<"Writing data to "<<tfilename<<" check it manually!"<<endl;
			cout<<"This is the symptom if individual names start with a number."<<endl<<"Try re-running with -X -Y"<<endl;
			output(tfilename);
            		return(false);
       		}
    	}
	return(true);
    }
    void applySuper(std::string name,std::vector<int> inds);///<Applies a superindividual merge, assigning "name" to the grouping
    void makeSuper(std::string superstring); ///< makes some super individuals as desribed by superstring
    void makeSuperFromFile(std::string filename); ///< makes some super individuals as desribed by superstring from a file name containing the string
    inline bool getIgnore(int i){
      if(i<0|| i>=(int)ignoresuper.size()) return(false);
      return(ignoresuper[i]);
    }
    inline int numIgnore(){
      return(numignore);
/*      int ret=0;
      for(unsigned int c1=0;c1<ignoresuper.size();c1++) ret+=ignoresuper[c1];
      return(ret);*/
    }
    inline int numSuper(){
      int nums=0;
      for(long i=0;i<(long) neff.size();i++) if(neff[i]>1) nums++;
      return(nums);
    }
    inline void setFileName(string newname){
      filename=newname;
    }
    inline string getFileName(){
      return(filename);
    }
    void addDataCfactor(string line);///< Adds the c factor from the datafile
    inline double getCfactor(){
      return(datacfactor);
    }
    void setCfactor(double c){ ///< set the c factor
      datacfactor=c;
    }
    void ensureCfactor(){ // Make sure that there is a valid cfactor
      if(datacfactor<=0) datacfactor=1.0;
    }
protected:
    inline std::string removequotes(std::string str){
	size_t found=str.find('"');
	while(found!=std::string::npos) {
		str.erase(found,1);
		found=str.find('"');
	}
	return(str);
    }
    unsigned int n;///<Number of individuals
    unsigned int k;///<Number of columns (only used if !square)
    bool square;///< Whether the data is a square matrix
    std::vector<std::vector<double> > data;///< data of copies
    std::vector<std::string> names;///<Names of the data
    std::vector<std::string> colnames;///<Names of the columns (used only if !square)
    std::vector<double> datarowsums;///<Data summed over rows
    std::vector<double> datacolsums;///<Data summed over columns
    std::vector<int> neff;///<Number of actual individuals used in this row
    std::vector<bool> ignoresuper;///< whether we should ignore for tree building
    double sumXall;///<Total number of copies
    double datacfactor;///<The c factor included in the data or force file
    int numignore;///<Total number of ignored super indivs
    string filename;///<The name of the data file that was read
};

} // end namespace fines
#endif
