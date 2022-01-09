#include "fsdonor.h"

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdexcept>

namespace fines
{
  /* fsutils.cpp has facilities to read id files.
*/

  
  FsDonor::FsDonor(std::string idfile,std::string donorfile,std::string cpdir,int runid,int nsamplefiles, int fsmode,bool verbose)
  {
    this->idfile=idfile;
    this->donorfile=donorfile;
    this->verbose=verbose;
    this->fsmode=fsmode;
    this->runid=runid;
    this->cpdir=cpdir;
    Nsamplefiles=nsamplefiles;
    Nrecip=0;
    Nraw=0;
    Nraw=0;
    realisedindex.clear();
    if(donorfile.compare("")!=0) {
      readDonorFile();
    }
    readIdFileContents();
    // read in the donor file
  }
  
  FsDonor::~FsDonor()
  {}

  void FsDonor::readIdFileContents(){
    
    allids=getIdsFromFile(idfile,true,0);
    usedids=getIdsFromFile(idfile,false,0);
    std::vector<std::string> allinclusionsasstr=getIdsFromFile(idfile,true,2);
    allinclusion.clear();
    for(unsigned int c1=0;c1<allinclusionsasstr.size();++c1){
      int tuse;
      istringstream ( allinclusionsasstr[c1] ) >> tuse;
      allinclusion.push_back(tuse);
    }
    try{
      allpops=getIdsFromFile(idfile,true,1);
      usedpops=getIdsFromFile(idfile,false,1);
    }catch(runtime_error& e){
      cerr<<"Provided donor file "<<idfile<<" does not contain population labels."<<endl;
      throw(runtime_error("fsproject: cannot create examples"));
    }
  
    Nraw = (int)allids.size();
    Nused = (int) usedids.size();
    indexPopulations();
  }

  void FsDonor::indexPopulations()
  {
    // first make a list of all unique population names
    if(recipnames.size()==0){
      popnames=unique(usedpops);
    }else{
      popnames=donornames;
      for(unsigned int c1=0;c1<usedpops.size();c1++){
	popnames.push_back(usedpops[c1]);
      }
      popnames=unique(popnames);
    }
    for(unsigned int c1=0;c1<popnames.size();c1++){
      popmap[popnames[c1]]=c1;
      popmembers.push_back(vector<int>());
    }
    
    // Loop over all individuals, and assign them to the appropriate population
    for(unsigned int c1=0;c1<allpops.size();++c1){
      if(allinclusion[c1]==1) {
	popmembers[popmap[allpops[c1]]].push_back(c1);
      }
    }
    
  }

  void FsDonor::createDirectory()
  {
    
    stringstream ss;
    ss<<cpdir<<"/stage5";
    ensureDirectory(ss.str());
    ss<<"/buildrefsample"<<runid;
    dirname=ss.str();
    ensureDirectory(dirname);
    if(runid==0){
      ss.str("");
      ss<<cpdir<<"/stage5/referenceidfiles";
      if(fsmode>0) ensureDirectory(ss.str());
      
    }
  }

  void FsDonor::createSampleDonorFile()
  {
    stringstream ssd;
    ssd <<cpdir<<"/stage5/externalsamples.donor";
    sampledonorfile=ssd.str();
    if(verbose) cout<<"Created reference donor file "<<sampledonorfile<<endl;
    filebuf fb;
    try{
      fb.open (sampledonorfile.c_str(),ios::out);
    }catch(exception &x){
      cerr<<"Error opening file!"<<endl<<x.what()<<endl; 
      throw(runtime_error("fsproject: cannot write to file!"));
    }
    ostream os (&fb);

    for(unsigned int c1=0;c1<popnames.size();++c1){
      os<<popnames[c1]<<" D"<<endl;
    }
    // 
    fb.close();
  }
  
  void FsDonor::createDonorFile()
  {
    stringstream ssd;
    ssd <<cpdir<<"_admixtureanalysis.donor";
    donorfile=ssd.str();
    if(verbose) cout<<"Created donor file "<<donorfile<<endl;
    filebuf fb;
    try{
      fb.open (donorfile.c_str(),ios::out);
    }catch(exception &x){
      cerr<<"Error opening file!"<<endl<<x.what()<<endl; 
      throw(runtime_error("fsproject: cannot write to file!"));
    }
    ostream os (&fb);

    dfiledonors.clear();
    dfilerecips.clear();
    donornames.clear();
    recipnames.clear();
    for(unsigned int c1=0;c1<popnames.size();++c1){
      os<<popnames[c1]<<" D"<<endl;
      dfiledonors[popnames[c1]]=1;
      donornames.push_back(popnames[c1]);
    }
    for(unsigned int c1=0;c1<popnames.size();++c1){
      os<<popnames[c1]<<" R"<<endl;
      dfilerecips[popnames[c1]]=1;
      recipnames.push_back(popnames[c1]);
    }
    // 
    fb.close();
  }

  void FsDonor::readDonorFile()
  {
    //      dfiledonors.push_back(popnames[c1]);
    //      dfiledonors.push_back(popnames[c1]);
    vector<string>tline;
    std::ifstream file(donorfile.c_str());
    string line;
    while(getline(file, line)) {
      tline=split(line," \t");
      if(tline.size()<=1){
	continue;
      }
      if(line.substr(0,1).compare("#")==0){
	continue;
      }
      if(tline[1].compare("D")==0){
	dfiledonors[tline[0]]=1;
	donornames.push_back(tline[0]);
	if(verbose) cout<<"From donorfile, read donor "<<tline[0]<<endl;
      }else if(tline[1].compare("R")==0){
	dfilerecips[tline[0]]=1;
	recipnames.push_back(tline[0]);
	if(verbose) cout<<"From donorfile, read recipient "<<tline[0]<<endl;
      }else{
	cerr<<"Donor line not reckognised: "<<line<<endl;
      }
      /*
      if ( popmap.find(tline[0]) == popmap.end() ) {
	// not found
	cerr<<"Error: A population called "<<tline[0]<<" was found in the donorfile "<<donorfile<<" but there are no individuals with this population label!"<<endl;
	throw(runtime_error("fsdonor: invalid donors supplied!"));
      }
      */
    }
    file.close();
    if(dfiledonors.size()==0) {
      cerr<<"ERROR: donorfile \""<<donorfile.c_str()<<"\" doesn't contain any donors!"<<endl;
      throw(runtime_error("data error"));
    }
    if(dfilerecips.size()==0) {
      cerr<<"ERROR: donorfile \""<<donorfile.c_str()<<"\" doesn't contain any recipients!"<<endl;
      throw(runtime_error("data error"));
    }
  }

    
  std::vector<int> FsDonor::sampleInclusion(int i)
  {
    int mypop=-1;
    if(i>=0) mypop=popmap[usedpops[i]];
    
    std::vector<int> ret=allinclusion;
    
    for(unsigned int c1=0;c1<recipnames.size();++c1){
      int refpop=popmap[popnames[c1]];
      if(mypop==refpop){ // remove self from donor pool, which is done for us
      }else{
	// sample one individual from the population to exclude
	int poprem = rand() % popmembers[c1].size();
	int indrem=popmembers[c1][poprem];
	ret[indrem]=0;
	if(verbose && i>=0) cout<<"Ind "<<i<<" called "<<usedids[i]<<" omitting ind "<<indrem<<" called "<<allids[indrem]<<" from population "<<c1<<" called "<<popnames[c1]<<endl;
      }
    }
    return(ret);
  }

  int FsDonor::countPosition(vector<int> myinclusion,string myname)
  {
    /// Count the position of this individual, as it will appear in the painting.
    // This has individuals ordered by population
    // as well as removing individuals from the model
    int ret=0;
    // This uses the ordering criteria for chromopainter which is a bit convoluted
    // We order the individuals by their population first, then their order in the file
    for(unsigned int tpop=0;tpop<recipnames.size();++tpop){ // for all populations
      string thispop=recipnames[tpop];
      for(unsigned int c1=0;c1<allpops.size();++c1){ // for all indivuals
	if((myinclusion[c1]>0) && allpops[c1].compare(thispop)==0) { // if including, and is this population
	  ++ret;
	  
	  if(allids[c1].compare(myname)==0) {// stop if we get to our own name
	    return(ret);
	  }
	}// end if
      }// end for over allpops
    }// end for over rpops
    return(ret);
  }
  
  void FsDonor::createIdFileForInd(int i)
  {

    stringstream ssid;

    if(i>=0) {
      int mypop=popmap[usedpops[i]];
      ssid<<dirname<<"/"<<"admixtureanalysis_"<<popnames[mypop]<<"_"<<usedids[i]<<"_idfile.ids";
      if(verbose) cout<<"Making sample idfile: "<<ssid.str()<<endl;
    }
    else {
      ssid<<cpdir<<"/stage5/referenceidfiles/admixtureanalysis_sample"<<-i<<"_idfile.ids";
      if(verbose) cout<<"Making reference idfile: "<<ssid.str()<<endl;
    }
    string idfile=ssid.str();

    filebuf fb;
    try{
      fb.open (idfile.c_str(),ios::out);
    }catch(exception &x){
      cerr<<"Error opening file!"<<endl<<x.what()<<endl; 
      throw(runtime_error("fsproject: cannot write to file!"));
    }
    ostream os (&fb);

    // Sample K-1 individuals from each population, except the one 
    std::vector<int> myinclusion=sampleInclusion(i);
    if(i>=0) {
      realisedindex.push_back(countPosition(myinclusion,usedids[i]));
      idfiles.push_back(idfile);
    }else{
      refidfiles.push_back(idfile);
    }
    
    for(int c1=0;c1<Nraw;++c1){
      os<<allids[c1]<<" "<<allpops[c1]<<" "<<myinclusion[c1]<<endl;
    }
    // 
    fb.close();
    //
  }
  
  void FsDonor::createRequiredFiles()
  {
    // Need to make the directory for all the files
    createDirectory();

    if(donorfile.compare("")==0) {
      createDonorFile();
    }
    createSampleDonorFile();
    countRecipients();
    for(int c1=0;c1<Nused;++c1){
      // Create the donor and ID file for each individual
      if(dfilerecips.count(usedpops[c1])) createIdFileForInd(c1);
    }
    if(fsmode>0 && runid==0){
      // Create donor files that can be used for samples
      for(int c1=1;c1<=Nsamplefiles;++c1){
	createIdFileForInd(-c1);
      }
    }
    
  }
  
  std::vector<std::string> FsDonor::getCommandContent()
  {
    //
    std::vector<std::string> cmds;
    // create the commands for each individual
    for(unsigned int c1=0;c1<(unsigned int) numRecipients();++c1){
      std::stringstream cmd;
      cmd<<"-t "<<idfiles[c1]<<" -f "<<donorfile<<" "<<realisedindex[c1]<<" "<<realisedindex[c1];
      cmds.push_back(cmd.str());
    }
    return(cmds);
  }

  void FsDonor::countRecipients()
  {
    int nrecips=0;
    for (std::map<string,int>::iterator it=dfilerecips.begin();
	 it!=dfilerecips.end();
	 ++it){
      int popindex=popmap.at(it->first); // we checked for existence in popindex alread
      nrecips+=popmembers[popindex].size();
    }
    if(verbose) cout<<"fsdonor: counted "<<nrecips<<" recipient individuals referenced in the donor file."<<endl;
    Nrecip=nrecips;
  }

  vector<string> FsDonor::recipLabels(){
    return(recipnames);
    /*    vector<string> v;
    for(map<string,int>::iterator it = dfilerecips.begin(); it != dfilerecips.end(); ++it) {
      v.push_back(it->first);
    }
    return(v);*/
  }

  vector<string> FsDonor::donorLabels(){
    return(donornames);
    /*    vector<string> v;
    for(map<string,int>::iterator it = dfiledonors.begin(); it != dfiledonors.end(); ++it) {
      v.push_back(it->first);
    }
    return(v);*/
  }
} // end namespace 
