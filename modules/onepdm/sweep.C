/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "sweeponepdm.h"

#ifdef MOLPRO
#include "global/CxOutputStream.h"
#define pout if (dmrginp.outputlevel() < 0) xout
#endif

namespace SpinAdapted{
void SweepOnepdm::BlockAndDecimate (SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, const bool &useSlater, const bool& dot_with_sys)
{
  //mcheck("at the start of block and decimate");
  // figure out if we are going forward or backwards
  dmrginp.guessgenT -> start();
  bool forward = (system.get_sites() [0] == 0);
  SpinBlock systemDot;
  SpinBlock envDot;
  int systemDotStart, systemDotEnd;
  int systemDotSize = sweepParams.get_sys_add() - 1;
  if (forward)
  {
    systemDotStart = *system.get_sites().rbegin () + 1;
    systemDotEnd = systemDotStart + systemDotSize;
  }
  else
  {
    systemDotStart = system.get_sites() [0] - 1;
    systemDotEnd = systemDotStart - systemDotSize;
  }
  vector<int> spindotsites(2); 
  spindotsites[0] = systemDotStart;
  spindotsites[1] = systemDotEnd;
  systemDot = SpinBlock(systemDotStart, systemDotEnd);

  SpinBlock environment, environmentDot, newEnvironment;

  int environmentDotStart, environmentDotEnd, environmentStart, environmentEnd;

  const int nexact = forward ? sweepParams.get_forward_starting_size() : sweepParams.get_backward_starting_size();

  system.addAdditionalCompOps();
  InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, sweepParams.get_sys_add(), dmrginp.direct(), DISTRIBUTED_STORAGE, true, true);
  
  InitBlocks::InitNewEnvironmentBlock(environment, systemDot, newEnvironment, system, systemDot,
				      sweepParams.get_sys_add(), sweepParams.get_env_add(), forward, dmrginp.direct(),
				      sweepParams.get_onedot(), nexact, useSlater, true, true, true);
  SpinBlock big;
  newSystem.set_loopblock(true);
  system.set_loopblock(false);
  newEnvironment.set_loopblock(false);
  InitBlocks::InitBigBlock(newSystem, newEnvironment, big); 

  const int nroots = dmrginp.nroots();
  std::vector<Wavefunction> solutions(nroots);

  
  for(int i=0;i<nroots;++i)
    {
      StateInfo newInfo;
      solutions[i].LoadWavefunctionInfo (newInfo, newSystem.get_sites(), i);
    }
  
#ifndef SERIAL
  mpi::communicator world;
  mpi::broadcast(world,solutions,0);
#endif

#ifdef SERIAL
  const int numprocs = 1;
#endif
#ifndef SERIAL
  const int numprocs = world.size();
#endif


  compute_onepdm(solutions, system, systemDot, newSystem, newEnvironment, big, numprocs);

}

double SweepOnepdm::do_one(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize)
{

  SpinBlock system;
  const int nroots = dmrginp.nroots();
  std::vector<double> finalEnergy(nroots,0.);
  std::vector<double> finalEnergy_spins(nroots,0.);
  double finalError = 0.;

  Matrix onepdm(2*dmrginp.last_site(), 2*dmrginp.last_site());onepdm=0.0;
  for (int i=0; i<nroots; i++)
    for (int j=0; j<=i; j++)      
      save_onepdm_binary(onepdm, i ,j);

  sweepParams.set_sweep_parameters();
  // a new renormalisation sweep routine
  pout << ((forward) ? "\t\t\t Starting renormalisation sweep in forwards direction" : "\t\t\t Starting renormalisation sweep in backwards direction") << endl;
  pout << "\t\t\t ============================================================================ " << endl;
  
  InitBlocks::InitStartingBlock (system,forward, sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), restartSize, restart, warmUp);

  sweepParams.set_block_iter() = 0;
 
  pout << "\t\t\t Starting block is :: " << endl << system << endl;

  SpinBlock::store (forward, system.get_sites(), system); // if restart, just restoring an existing block --
  sweepParams.savestate(forward, system.get_sites().size());
  bool dot_with_sys = true;

  sweepParams.set_guesstype() = TRANSPOSE;
	  
  SpinBlock newSystem;
  BlockAndDecimate (sweepParams, system, newSystem, warmUp, dot_with_sys);
  pout.precision(12);
  pout << "\t\t\t The lowest sweep energy : "<< sweepParams.get_lowest_energy()[0]+dmrginp.get_coreenergy()<<endl;
  pout << "\t\t\t ============================================================================ " << endl;

  for (int i=0; i<nroots; i++)
    for (int j=0; j<=i; j++) {
      load_onepdm_binary(onepdm, i ,j);
      accumulate_onepdm(onepdm);
      save_onepdm_spatial_text(onepdm, i ,j);
      save_onepdm_text(onepdm, i ,j);
      save_onepdm_spatial_binary(onepdm, i ,j);
    }
  return sweepParams.get_lowest_energy()[0];
}
}
