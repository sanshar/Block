#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#ifndef SERIAL
#include <boost/mpi.hpp>
#endif
#include "GAInput.h"
#include "GAOptimize.h"
#include "ReadIntegral.h"
using namespace std;

#ifndef SERIAL
namespace mpi = boost::mpi;
#endif

namespace genetic
{
  GAInput gainput;
  int Gene::m_length = 0;
  function<double(const Gene&)> Cell::Evaluate;
};

void gaordering(ifstream& confFile, ifstream& dumpFile, vector<int>& reorder, vector<int>& oldtonew)
{
  srand(genetic::gainput.random_seed);
  cout << "Random Seed = " << genetic::gainput.random_seed << endl;

  ifstream fdump(dumpFile.c_str());
  Matrix K; genetic::ReadIntegral(fdump, K);

  genetic::Cell final = genetic::Optimize(K);
  cout << "##################### MINIMUM GENE REP. #####################" << endl;
  cout << "Gene with MinValue = " << final << endl;
  cout << "Effective Distance = " << sqrt(final.Fitness()) << endl;

  cout << "#################### DMRG REORDER FORMAT ####################" << endl;
  vector<int> ordering(final.Gen().Sequence());
  for(int i = 0; i < genetic::Gene::Length()-1; ++i) cout << ordering[i] + 1 << ",";
  cout << ordering[genetic::Gene::Length()-1] + 1 << endl;

  return 0;
}
