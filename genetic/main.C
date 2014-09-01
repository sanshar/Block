#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "GAInput.h"
#include "GAOptimize.h"
#include "fiedler.h"
#include "pario.h"
using namespace std;

#ifndef SERIAL
#include <boost/mpi.hpp>
namespace mpi = boost::mpi;
#endif

int main(int argc, char* argv[])
{
#ifndef SERIAL
  mpi::environment env(argc, argv);
  mpi::communicator world;
  if(world.rank() == 0) pout << "Parallel GA simulation" << endl;
#endif

  string confFileName;
  string dumpFileName;

  for(int i = 1; i < argc; ++i)
  {
    if(strcmp(argv[i], "-config")   == 0) confFileName = argv[++i];
    if(strcmp(argv[i], "-integral") == 0) dumpFileName = argv[++i];
  }

  ifstream confFile(confFileName.c_str());
  ifstream dumpFile(dumpFileName.c_str());

  std::vector<int> fiedlerv = get_fiedler(dumpFileName, dumpFile);

  genetic::Cell final = genetic::gaordering(confFile, dumpFile, fiedlerv);

#ifndef SERIAL
  if(world.rank() == 0)
#endif
  {
    pout << "##################### MINIMUM GENE REP. #####################" << endl;
    pout << "Gene with MinValue = " << final << endl;
    pout << "Effective Distance = " << sqrt(final.Fitness()) << endl;

    pout << "#################### DMRG REORDER FORMAT ####################" << endl;
    int n = genetic::Gene::Length() - 1;
    vector<int> gaorder(final.Gen().Sequence());

    for(int i = 0; i < n; ++i) pout << gaorder[i] + 1 << ",";
    pout << gaorder[n] + 1 << endl;
  }

  return 0;
}
