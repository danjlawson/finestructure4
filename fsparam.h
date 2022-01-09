#ifndef FSPARAM_H
#define FSPARAM_H
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <fstream>
#include <algorithm>
#include <string.h>

#include "finestructure/safegetline.h" // for cross-platform safely getting lines

using namespace std;

namespace fines
{

/**
    @brief Parameter structure for the fs pipeline
*/
class FsParam
{
public:
  FsParam();
  ~FsParam();
};

} // end namespace fines
#endif
