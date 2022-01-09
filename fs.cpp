#include <string>
#include <cstring>
#include <fstream>
#include <sstream>

#include "cp/ChromoPainterMutEM.h"
#include "finestructure/fines.h"
#include "chromocombine/ChromoCombine.h"
#include "fsproject.h"
#include "fsconstants.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


using namespace std;

namespace fines
{

/////////////////////////////////////
string getFsVersion(){
	string ret;
	#ifdef PACKAGE_STRING // Omitted due to version matching issues
	ret.append(PACKAGE_STRING);
	ret.append("\n");
	ret.append("Report all bugs to ");
	ret.append(PACKAGE_BUGREPORT);
	ret.append("\n");
	
	#endif
	ret.append("fs");
	ret.append(" build date "); 
	ret.append(__DATE__);
	ret.append(" at ");
	ret.append(__TIME__);
	return(ret);
}

} // end namespace fines


using namespace fines;

int main(int argc, char *argv[])
{
    string comment="Command line: ";
    for(int c1=0;c1< argc;c1++) {comment.append(argv[c1]);comment.append(" ");}
    comment.append("\nVersion: ");
    comment.append(getFsVersion());

    if (argc<2) {return(fsproject(1,argc,&argv[0]));} // this is effectively help
        
// process the main arguments
/*    if(strcmp("project",argv[1]) ==0 ) {
      cout<<"Running in project mode"<<endl;
      return(fsproject(argc-1,&argv[1]));
      }*/
    if((strcmp("fs",argv[1]) ==0 )||(strcmp("finestructure",argv[1]) ==0 )) {
      cout<<"Running in fs mode"<<endl;
      return(finestructure(argc-1,&argv[1]));
    }else if((strcmp("cp",argv[1]) ==0 )||(strcmp("chromopainter",argv[1]) ==0 )) {
      cout<<"Running in cp mode"<<endl;
      return(chromopainter(argc-1,&argv[1]));
    }else if((strcmp("combine",argv[1]) ==0 )||(strcmp("chromocombine",argv[1]) ==0 )) {
      cout<<"Running in combine mode"<<endl;
      return(chromocombine(argc-1,&argv[1]));
    }else if((strcmp("version",argv[1]) ==0 )||(strcmp("-V",argv[1]) ==0 )||(strcmp("--V",argv[1]) ==0 )||(strcmp("--version",argv[1]) ==0)||(strcmp("-version",argv[1]) ==0 )) {
      cout<<getFsVersion()<<endl;return 0;
    /*}else if((strcmp("help",argv[1]) ==0 )||(strcmp("-h",argv[1]) ==0 )||(strcmp("--h",argv[1]) ==0 )||(strcmp("--help",argv[1]) ==0)||(strcmp("-help",argv[1]) ==0 )) {
       cout<<fshelp<<endl;return 0;*/
    }else if((strcmp("tools",argv[1]) ==0 )) {
       cout<<fstoolshelp<<endl;return 0;
    }else{
      //      cout<<"Running in automatic mode"<<endl;
      return(fsproject(1,argc,&argv[0]));
    }
    return(0);
}
