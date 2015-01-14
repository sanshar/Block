#include <iostream>
#include <vector>
#include <algorithm>
#include "Evaluate.h"
#include "GAOptimize.h"
#include "Generation.h"
#include "GAInput.h"
#include "ReadIntegral.h"
using namespace std;

#include <boost/function.hpp>
#include <boost/bind.hpp>

#ifndef SERIAL
#include <boost/mpi.hpp>
namespace mpi = boost::mpi;
#endif

#include <newmat.h>
#include <newmatutils.h>

namespace genetic
{
  GAInput gainput;
  int Gene::m_length = 0;
  boost::function<double(const Gene&)> Cell::Evaluate;
};

genetic::Cell genetic::gaordering(ifstream& confFile, ifstream& dumpFile, std::vector<int> fiedlerorder)
{
#ifndef SERIAL
  mpi::communicator world;
#endif

  Matrix K;
  double ksum = 0.0;

#ifndef SERIAL
  if(world.rank() == 0)
  {
#endif
    if(confFile.is_open()) gainput.Configure(confFile);
    ReadIntegral(dumpFile, K);
    Gene::Length() = K.Nrows();
    if(gainput.max_cells == 0) gainput.max_cells = 2 * Gene::Length();
    for(int i = 0; i < K.Nrows(); ++i)
      for(int j = i + 1; j < K.Ncols(); ++j) ksum += K.element(i, j);
#ifndef SERIAL
  }
  mpi::broadcast(world, fiedlerorder, 0);
  mpi::broadcast(world, gainput, 0);
  mpi::broadcast(world, Gene::Length(), 0);
  mpi::broadcast(world, K, 0);
  mpi::broadcast(world, ksum, 0);
#endif

  Cell::Evaluate = boost::bind(genetic::Evaluate, 1.0/ksum, gainput.exponent, _1, K);
  Cell best;
#ifndef SERIAL
  int ntask = 1 + gainput.max_community / world.size();

  Cell comm_best = gaoptimize(genetic::gainput.random_seed+world.rank(), fiedlerorder);
  cout << "Order #" << world.rank() << ": " << comm_best << endl;
  for(int i = 1; i < ntask; ++i)
  {
    Cell comm_cell = gaoptimize(genetic::gainput.random_seed + i * world.size() + world.rank(), fiedlerorder);
    cout << "Order #" << i * world.size() + world.rank() << ": " << comm_cell << endl;
    if(comm_cell < comm_best) comm_best = comm_cell;
  }

  if(world.rank() == 0)
    mpi::reduce(world, comm_best, best, mpi::minimum<Cell>(), 0);
  else
    mpi::reduce(world, comm_best,       mpi::minimum<Cell>(), 0);

#else
  int ntask = gainput.max_community;
  best = gaoptimize(genetic::gainput.random_seed, fiedlerorder);
  cout << "Order #" << 0 << ": " << best << endl;
  for(int i = 1; i < ntask; ++i)
  {
    Cell comm_cell = gaoptimize(genetic::gainput.random_seed+i, fiedlerorder);
    cout << "Order #" << i << ": " << comm_cell << endl;
    if(comm_cell < best) best = comm_cell;
  }

#endif

  return best;
}

genetic::Cell genetic::gaoptimize(const int& seed, std::vector<int> fiedlerorder)
{
  srand(seed);
  Generation ancestor;
  if (gainput.fiedler==1)
     ancestor.AddFiedler(fiedlerorder);

  for(int g = 0; g < gainput.max_generation; ++g)
  {
    Generation nextgen;
    nextgen.Generate(ancestor);
    ancestor = nextgen;
  }

  return ancestor.Min();
}
