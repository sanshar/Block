#include "sweep.h"
#include "global.h"
#include "solver.h"
#include "initblocks.h"
#include "initblocks.h"
#include "MatrixBLAS.h"
#ifndef SERIAL
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif

using namespace boost;
using namespace std;


void SpinAdapted::Sweep::fullci(double sweep_tol)
{
  SweepParams sweepParams;
  sweepParams.set_sweep_parameters();
  SpinBlock system;
  InitBlocks::InitStartingBlock(system, true, sweepParams.get_forward_starting_size(),  sweepParams.get_backward_starting_size(), 0, false, true);

  int numsites = dmrginp.last_site();
  int forwardsites = numsites/2+numsites%2;
  int backwardsites = numsites - forwardsites;
  SpinQuantum hq(0,0,IrrepSpace(0));

  for (int i=0; i<forwardsites-1; i++) {
    SpinBlock sysdot(i+1, i+1);
    cout <<sysdot<<endl;
    SpinBlock newSystem;
    system.addAdditionalCompOps();
    newSystem.default_op_components(false, system, sysdot, false, true);
    newSystem.setstoragetype(DISTRIBUTED_STORAGE);
    newSystem.BuildSumBlock (NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, system, sysdot);
    system = newSystem;
    pout << system<<endl;
  }

  SpinBlock environment;
  InitBlocks::InitStartingBlock(environment, false, sweepParams.get_forward_starting_size(),  sweepParams.get_backward_starting_size(), 0, false, true);
  cout << environment<<endl;
  for (int i=0;i <backwardsites-1; i++) {
    SpinBlock envdot(numsites-2-i, numsites-2-i);
    cout << envdot<<endl;
    SpinBlock newEnvironment;
    environment.addAdditionalCompOps();
    newEnvironment.default_op_components(false, environment, envdot, true, true);
    newEnvironment.setstoragetype(DISTRIBUTED_STORAGE);
    newEnvironment.BuildSumBlock (NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, environment, envdot);
    environment = newEnvironment;
    pout << environment<<endl;
  }

  SpinBlock big;
  InitBlocks::InitBigBlock(system, environment, big); 
  int nroots = dmrginp.nroots(0);
  std::vector<Wavefunction> solution(nroots);
  std::vector<double> energies(nroots);
  double tol = sweepParams.get_davidson_tol();
  DiagonalMatrix e;
  Solver::solve_wavefunction(solution, energies, big, tol, BASIC, false, true, false, e, sweepParams.get_additional_noise());
  for (int i=0; i<nroots; i++) {
    pout << "fullci energy "<< energies[i]<<endl;
  }
}
