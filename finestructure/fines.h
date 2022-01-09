/*
 * fines.h
 *
 *  Created on: Jun 22, 2009
 *      Author: dan.lawson@bristol.ac.uk
 */

#ifndef FINES_H_
#define FINES_H_

#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <cstring>
#include <sstream>
#include <fstream>
#include <stdio.h>     /* for printf */
#include <stdlib.h>    /* for exit */
#include <unistd.h>    /* for getopt */

static const char * const fssuccesstext="FineStructure completed successfully.";

#define BETAMOD_EQUI 1
#define BETAMOD_CONST 2
#define BETAMOD_F 3
#define BETAMOD_F2 4
#define BETAMOD_COPYMAT 5
#define BETAMOD_F2_COPYMAT 6

#define INFDATA_COUNTS -1 
//<0 means ignore the lengths
#define INFDATA_LENGTHS 1
#define INFDATA_ALL 2
#define INFDATA_TOTALLENGTHS 3
#define INFDATA_ALLNOTLENGTHS 4

#define MODELTYPE_FINESTRUCTURE 1
#define MODELTYPE_NORMALISED 2
#define MODELTYPE_NORMALISEDMERGEONLY 3
#define MODELTYPE_INDIVIDUAL 4

#define TREEMOD_NOFLATTEN 0
#define TREEMOD_FLATTEN 1
#define TREEMOD_FLATTENHILLCLIMB 2

namespace fines
{

class ProgramOptions
{
public:
	ProgramOptions() :
		burnin(1000),
		additional(1000),
		thinin(1),
		verbose(0),
		silent(0),
		test_max(1500),
		userefonly(0),
		pcaprob(0.0),
		fixK(false),
		treescale(0),
		method("OMCMC"),
		usedata("Counts"),
		outfile(""){};
	int burnin;
	int additional;
	int thinin;
	int verbose;
	int silent;
	int test_max;
	int userefonly;
	double pcaprob;
	bool fixK;
	int treescale;
	std::string method;
	std::string usedata;
	std::string outfile;
};

 int finestructure(int argc, char *argv[]); //< The main finestructure utility

 std::vector<double> &splitasdouble(const std::string &s, char delim, std::vector<double> &elems);// A utility function
 

} // end namespace fines

#endif /* FINES_H_ */
