#ifndef GA_OPTIMIZE_H
#define GA_OPTIMIZE_H

#include <fstream>
#include "Cell.h"

namespace genetic
{
  Cell gaordering(ifstream& confFile, ifstream& dumpFile, std::vector<int> fiedlerorder);
  Cell gaoptimize(const int& seed, std::vector<int> fiedlerorder); 
};

#endif
