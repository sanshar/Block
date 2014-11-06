/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "sweep.h"
#include "global.h"
#include "solver.h"
#include "initblocks.h"
#include "initblocks.h"
#include "MatrixBLAS.h"
#include <boost/format.hpp>
#include "op_components.h"
#ifndef SERIAL
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif

using namespace boost;
using namespace std;


void SpinAdapted::Sweep::fullci(double sweep_tol)
{
  int integralIndex = 0;
  SweepParams sweepParams;
  sweepParams.set_sweep_parameters();

  SpinBlock system;
  InitBlocks::InitStartingBlock(system, true, 0, 0, sweepParams.get_forward_starting_size(),  sweepParams.get_backward_starting_size(), 0, false, true, integralIndex);
  int numsites = dmrginp.spinAdapted() ? dmrginp.last_site() : dmrginp.last_site()/2;
  int forwardsites = numsites/2+numsites%2;
  int backwardsites = numsites - forwardsites;
  SpinQuantum hq(0,SpinSpace(0),IrrepSpace(0));


  for (int i=0; i<forwardsites-1; i++) {
    SpinBlock sysdot(i+1, i+1, integralIndex, true);
    SpinBlock newSystem;
    system.addAdditionalCompOps();
    newSystem.set_integralIndex() = integralIndex;
    newSystem.default_op_components(false, system, sysdot, false, true, true);
    newSystem.setstoragetype(DISTRIBUTED_STORAGE);
    newSystem.BuildSumBlock (NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, system, sysdot);

    system = newSystem;

  }

  SpinBlock environment;
  InitBlocks::InitStartingBlock(environment, false, 0, 0, sweepParams.get_forward_starting_size(),  sweepParams.get_backward_starting_size(), 0, false, true, integralIndex);
  pout << environment<<endl;
  for (int i=0;i <backwardsites-1; i++) {
    SpinBlock envdot(numsites-2-i, numsites-2-i, integralIndex, true);
    SpinBlock newEnvironment;
    environment.addAdditionalCompOps();
    newEnvironment.set_integralIndex() = integralIndex;
    newEnvironment.default_op_components(false, environment, envdot, true, true, true);
    newEnvironment.setstoragetype(DISTRIBUTED_STORAGE);
    newEnvironment.BuildSumBlock (NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, environment, envdot);

    environment = newEnvironment;
  }

  pout <<"\t\t\t System Block :: "<< system;
  pout <<"\t\t\t Environment Block :: "<< environment;
  SpinBlock big;
  InitBlocks::InitBigBlock(system, environment, big); 

  //SpinBlock systemp(0,1), envtmp(2,3);
  //InitBlocks::InitBigBlock(systemp, envtmp, big); 

  int nroots = dmrginp.nroots(0);
  std::vector<Wavefunction> solution(nroots);
  std::vector<double> energies(nroots);
  double tol = sweepParams.get_davidson_tol();

  pout << "\t\t\t Solving the Wavefunction "<<endl;
  int currentState = 0;
  std::vector<Wavefunction> lowerStates;
  Solver::solve_wavefunction(solution, energies, big, tol, BASIC, false, true, false, sweepParams.get_additional_noise(), currentState, lowerStates);
  for (int i=0; i<nroots; i++) {
    pout << "fullci energy "<< energies[i]<<endl;
  }
  if (!mpigetrank())
  {
#ifndef MOLPRO
    FILE* f = fopen("dmrg.e", "wb");
#else
    std::string efile;
    efile = str(boost::format("%s%s") % dmrginp.load_prefix() % "/dmrg.e" );
    FILE* f = fopen(efile.c_str(), "wb");
#endif
    
    for(int j=0;j<nroots;++j) {
      double e = energies[j]; 
      fwrite( &e, 1, sizeof(double), f);
    }
    fclose(f);
  }

}

void SpinAdapted::Sweep::tiny(double sweep_tol)
{
#ifndef SERIAL
  if(mpigetrank() == 0) {
#endif
    pout.precision(12);

  int nroots = dmrginp.nroots(0);
  SweepParams sweepParams;
  sweepParams.set_sweep_parameters();
  SpinBlock system(0,dmrginp.last_site()-1, 0, true);
  const StateInfo& sinfo = system.get_stateInfo();
  SpinQuantum hq(0,SpinSpace(0),IrrepSpace(0));
  for (int i=0; i<sinfo.totalStates; i++) {
    if (sinfo.quanta[i] == dmrginp.molecule_quantum()) {
      Matrix& h = system.get_op_rep(HAM, hq)->operator_element(i,i);
      DiagonalMatrix energies(h.Nrows()); energies = 0.0;
      Matrix vec(h.Nrows(), h.Ncols()); vec = 0.0;
      diagonalise(h, energies, vec);
      
      for (int x=0; x<nroots; x++) 
	pout << "fullci energy  "<< energies(x+1)<<endl;

      if (mpigetrank() == 0)
      {
#ifndef MOLPRO
	FILE* f = fopen("dmrg.e", "wb");
#else
	std::string efile;
	efile = str(boost::format("%s%s") % dmrginp.load_prefix() % "/dmrg.e" );
	FILE* f = fopen(efile.c_str(), "wb");
#endif
	
	for(int j=0;j<nroots;++j) {
	  double e = energies(j+1); 
	  fwrite( &e, 1, sizeof(double), f);
	}
	fclose(f);
      }

      return;
    }
  }

  pout << "The wavefunction symmetry is not possible with the orbitals supplied."<<endl;
  abort();
#ifndef SERIAL
  }
#endif
}

