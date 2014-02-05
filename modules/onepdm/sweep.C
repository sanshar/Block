/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "density.h"
#include "sweeponepdm.h"
#include "pario.h"
#include "onepdm.h"

namespace SpinAdapted{
void SweepOnepdm::BlockAndDecimate (SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, const bool &useSlater, const bool& dot_with_sys, int state)
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
    systemDotStart = dmrginp.spinAdapted() ? *system.get_sites().rbegin () + 1 : (*system.get_sites().rbegin ())/2 + 1 ;
    systemDotEnd = systemDotStart + systemDotSize;
  }
  else
  {
    systemDotStart = dmrginp.spinAdapted() ? system.get_sites()[0] - 1 : (system.get_sites()[0])/2 - 1 ;
    systemDotEnd = systemDotStart - systemDotSize;
  }
  vector<int> spindotsites(2); 
  spindotsites[0] = systemDotStart;
  spindotsites[1] = systemDotEnd;
  systemDot = SpinBlock(systemDotStart, systemDotEnd);

  SpinBlock environment, environmentDot, newEnvironment;
  int environmentDotStart, environmentDotEnd, environmentStart, environmentEnd;

  const int nexact = forward ? sweepParams.get_forward_starting_size() : sweepParams.get_backward_starting_size();
  
  InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, sweepParams.get_sys_add(), dmrginp.direct(), DISTRIBUTED_STORAGE_FOR_ONEPDM, true, true);
  
  InitBlocks::InitNewEnvironmentBlock(environment, systemDot, newEnvironment, system, systemDot,
				      sweepParams.get_sys_add(), sweepParams.get_env_add(), forward, dmrginp.direct(),
				      sweepParams.get_onedot(), nexact, useSlater, true, true, true);
  SpinBlock big;
  newSystem.set_loopblock(true);
  system.set_loopblock(false);
  newEnvironment.set_loopblock(false);
  InitBlocks::InitBigBlock(newSystem, newEnvironment, big); 

  const int nroots = dmrginp.nroots();
  std::vector<Wavefunction> solution(1);

  DiagonalMatrix e;
  GuessWave::guess_wavefunctions(solution[0], e, big, sweepParams.get_guesstype(), true, state, true, 0.0); 

#ifndef SERIAL
  mpi::communicator world;
  mpi::broadcast(world, solution, 0);
#endif

  std::vector<Matrix> rotateMatrix;
  DensityMatrix tracedMatrix;
  tracedMatrix.allocate(newSystem.get_stateInfo());
  tracedMatrix.makedensitymatrix(solution, big, std::vector<double>(1,1.0), 0.0, 0.0, false);
  rotateMatrix.clear();
  if (!mpigetrank())
    double error = newSystem.makeRotateMatrix(tracedMatrix, rotateMatrix, sweepParams.get_keep_states(), sweepParams.get_keep_qstates());
  

#ifndef SERIAL
  mpi::broadcast(world,rotateMatrix,0);
#endif
#ifdef SERIAL
  const int numprocs = 1;
#endif
#ifndef SERIAL
  const int numprocs = world.size();
#endif

  Matrix onepdm;
  load_onepdm_binary(onepdm, state ,state);

  if (sweepParams.get_block_iter() == 0) {
    //this is inface a combination of  2_0_0, 1_1_0 and 0_2_0
    pout << "compute 2_0_0"<<endl;
    compute_one_pdm_2_0_0(solution[0], solution[0], big, onepdm);
    pout << "compute 1_1_0"<<endl;
    compute_one_pdm_1_1_0(solution[0], solution[0], big, onepdm);
  }

  pout << "compute 0_2_0"<<endl;
  compute_one_pdm_0_2_0(solution[0], solution[0], big, onepdm);
  pout << "compute 1_1"<<endl;
  compute_one_pdm_1_1(solution[0], solution[0], big, onepdm);

  if (sweepParams.get_block_iter()  == sweepParams.get_n_iters() - 1) {
    pout << "compute 0_2"<<endl;
    compute_one_pdm_0_2(solution[0], solution[0], big, onepdm);
  }

  accumulate_onepdm(onepdm);
  save_onepdm_binary(onepdm, state, state);

  SaveRotationMatrix (newSystem.get_sites(), rotateMatrix, state);

  solution[0].SaveWavefunctionInfo (big.get_stateInfo(), big.get_leftBlock()->get_sites(), state);

  newSystem.transform_operators(rotateMatrix);

}

double SweepOnepdm::do_one(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize, int state)
{

  SpinBlock system;
  const int nroots = dmrginp.nroots();
  std::vector<double> finalEnergy(nroots,0.);
  std::vector<double> finalEnergy_spins(nroots,0.);
  double finalError = 0.;

  int pdmsize = dmrginp.spinAdapted() ? 2*dmrginp.last_site() : dmrginp.last_site();
  Matrix onepdm(pdmsize, pdmsize);onepdm=0.0;

  save_onepdm_binary(onepdm, state ,state);

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
  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); )
    {
      pout << "\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
      pout << "\t\t\t ----------------------------" << endl;
      if (forward)
	pout << "\t\t\t Current direction is :: Forwards " << endl;
      else
	pout << "\t\t\t Current direction is :: Backwards " << endl;

      if (sweepParams.get_block_iter() == 0)
	sweepParams.set_guesstype() = TRANSPOSE;
      else
	sweepParams.set_guesstype() = TRANSFORM;

      pout << "\t\t\t Blocking and Decimating " << endl;

      SpinBlock newSystem;
      BlockAndDecimate (sweepParams, system, newSystem, warmUp, dot_with_sys, state);
      pout.precision(12);

      system = newSystem;

      pout << system<<endl;
      
      SpinBlock::store (forward, system.get_sites(), system);	 	

      pout << "\t\t\t saving state " << system.get_sites().size() << endl;
      ++sweepParams.set_block_iter();
      //sweepParams.savestate(forward, system.get_sites().size());
    }
  pout << "\t\t\t The lowest sweep energy :  ----"<<endl;
  pout << "\t\t\t ============================================================================ " << endl;


  load_onepdm_binary(onepdm, state ,state);
  accumulate_onepdm(onepdm);
  save_onepdm_spatial_text(onepdm, state, state);
  save_onepdm_text(onepdm, state, state);
  save_onepdm_spatial_binary(onepdm, state, state);

  return sweepParams.get_lowest_energy()[0];
}
}
