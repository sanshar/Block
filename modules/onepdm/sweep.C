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
  systemDot = SpinBlock(systemDotStart, systemDotEnd, system.get_integralIndex(), true);

  SpinBlock environment, environmentDot, newEnvironment;
  int environmentDotStart, environmentDotEnd, environmentStart, environmentEnd;

  const int nexact = forward ? sweepParams.get_forward_starting_size() : sweepParams.get_backward_starting_size();
  
  newSystem.set_integralIndex() = system.get_integralIndex();
  newSystem.default_op_components(dmrginp.direct(), system, systemDot, false, false, true);
  newSystem.erase(CRE_CRE_DESCOMP);
  newSystem.erase(CRE_CRE);
  newSystem.erase(HAM);
  newSystem.setstoragetype(DISTRIBUTED_STORAGE_FOR_ONEPDM);
  newSystem.BuildSumBlock (NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, system, systemDot);
  if (dmrginp.outputlevel() > 0) {
    pout << "\t\t\t NewSystem block " << endl << newSystem << endl;
    newSystem.printOperatorSummary();
  }

  
  InitBlocks::InitNewEnvironmentBlock(environment, systemDot, newEnvironment, system, systemDot, sweepParams.current_root(), sweepParams.current_root(),
				      sweepParams.get_sys_add(), sweepParams.get_env_add(), forward, dmrginp.direct(),
				      sweepParams.get_onedot(), nexact, useSlater, system.get_integralIndex(), false, false, true);
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
  DensityMatrix tracedMatrix(newSystem.get_stateInfo());
  tracedMatrix.allocate(newSystem.get_stateInfo());
  tracedMatrix.makedensitymatrix(solution, big, std::vector<double>(1,1.0), 0.0, 0.0, false);
  rotateMatrix.clear();
  if (!mpigetrank())
    double error = makeRotateMatrix(tracedMatrix, rotateMatrix, sweepParams.get_keep_states(), sweepParams.get_keep_qstates());
  

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
  Matrix pairmat;
  if (dmrginp.hamiltonian() == BCS)
    load_pairmat_binary(pairmat, state ,state);

  if (sweepParams.get_block_iter() == 0) {
    //this is inface a combination of  2_0_0, 1_1_0 and 0_2_0
    p2out << "\t\t\t compute 2_0_0"<<endl;
    compute_one_pdm_2_0_0(solution[0], solution[0], big, onepdm);
    if (dmrginp.hamiltonian() == BCS)
      compute_pair_2_0_0(solution[0], solution[0], big, pairmat);
    p2out << "\t\t\t compute 1_1_0"<<endl;
    compute_one_pdm_1_1_0(solution[0], solution[0], big, onepdm);
    if (dmrginp.hamiltonian() == BCS)    
      compute_pair_1_1_0(solution[0], solution[0], big, pairmat);
  }

  p2out << "\t\t\t compute 0_2_0"<<endl;
  compute_one_pdm_0_2_0(solution[0], solution[0], big, onepdm);
  if (dmrginp.hamiltonian() == BCS)  
    compute_pair_0_2_0(solution[0], solution[0], big, pairmat);  
  p2out << "\t\t\t compute 1_1"<<endl;
  compute_one_pdm_1_1(solution[0], solution[0], big, onepdm);
  if (dmrginp.hamiltonian() == BCS)  
    compute_pair_1_1(solution[0], solution[0], big, pairmat);

  if (sweepParams.get_block_iter()  == sweepParams.get_n_iters() - 1) {
    p2out << "\t\t\t compute 0_2"<<endl;
    compute_one_pdm_0_2(solution[0], solution[0], big, onepdm);
    if (dmrginp.hamiltonian() == BCS)    
      compute_pair_0_2(solution[0], solution[0], big, pairmat);    
  }

  accumulate_onepdm(onepdm);
  save_onepdm_binary(onepdm, state, state);

  if (dmrginp.hamiltonian() == BCS) {
    accumulate_onepdm(pairmat);
    save_pairmat_binary(pairmat, state, state);
  }

  SaveRotationMatrix (newSystem.get_sites(), rotateMatrix, state);

  solution[0].SaveWavefunctionInfo (big.get_stateInfo(), big.get_leftBlock()->get_sites(), state);

  newSystem.transform_operators(rotateMatrix);

}

double SweepOnepdm::do_one(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize, int state)
{
  Timer sweeptimer;
  int integralIndex = 0;
  SpinBlock system;
  const int nroots = dmrginp.nroots();
  std::vector<double> finalEnergy(nroots,0.);
  std::vector<double> finalEnergy_spins(nroots,0.);
  double finalError = 0.;

  int pdmsize = dmrginp.spinAdapted() ? 2*dmrginp.last_site() : dmrginp.last_site();
  Matrix onepdm(pdmsize, pdmsize);onepdm=0.0;
  Matrix pairmat;
  if (dmrginp.hamiltonian() == BCS) {
    pairmat.ReSize(pdmsize, pdmsize);
    pairmat = 0.0;
    save_pairmat_binary(pairmat, state, state);
  }

  save_onepdm_binary(onepdm, state ,state);

  sweepParams.set_sweep_parameters();
  // a new renormalisation sweep routine
  pout << ((forward) ? "\t\t\t Starting renormalisation sweep in forwards direction" : "\t\t\t Starting renormalisation sweep in backwards direction") << endl;
  pout << "\t\t\t ============================================================================ " << endl;
  
  InitBlocks::InitStartingBlock (system,forward, sweepParams.current_root(), sweepParams.current_root(), sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), restartSize, restart, warmUp, integralIndex);

  sweepParams.set_block_iter() = 0;
 
  pout << "\t\t\t Starting block is :: " << endl << system << endl;

  SpinBlock::store (forward, system.get_sites(), system, sweepParams.current_root(), sweepParams.current_root()); // if restart, just restoring an existing block --
  sweepParams.savestate(forward, system.get_sites().size());
  bool dot_with_sys = true;

  sweepParams.set_guesstype() = TRANSPOSE;
  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); )
    {
      pout << "\n\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
      pout << "\t\t\t ----------------------------" << endl;
      if (forward)
	{ p1out << "\t\t\t Current direction is :: Forwards " << endl; }
      else
	{ p1out << "\t\t\t Current direction is :: Backwards " << endl; }

      if (sweepParams.get_block_iter() == 0)
	sweepParams.set_guesstype() = TRANSPOSE;
      else
	sweepParams.set_guesstype() = TRANSFORM;

      p1out << "\t\t\t Blocking and Decimating " << endl;

      SpinBlock newSystem;
      BlockAndDecimate (sweepParams, system, newSystem, warmUp, dot_with_sys, state);
      pout.precision(12);

      system = newSystem;

      pout << system<<endl;
      
      SpinBlock::store (forward, system.get_sites(), system, sweepParams.current_root(), sweepParams.current_root());	 	

      p1out << "\t\t\t saving state " << system.get_sites().size() << endl;
      ++sweepParams.set_block_iter();
      //sweepParams.savestate(forward, system.get_sites().size());
    }
  pout << "\t\t\t The lowest sweep energy : "<< sweepParams.get_lowest_energy()[0] << endl;
  pout << "\t\t\t ============================================================================ " << endl;


  load_onepdm_binary(onepdm, state ,state);
  accumulate_onepdm(onepdm);
  save_onepdm_spatial_text(onepdm, state, state);
  save_onepdm_text(onepdm, state, state);
  save_onepdm_spatial_binary(onepdm, state, state);

  if (dmrginp.hamiltonian() == BCS) {
    load_pairmat_binary(pairmat, state, state);
    accumulate_onepdm(pairmat);
    // FIXME write out text version
    // only <D{ia}D{jb}> is in the matrix
    save_pairmat_text(pairmat , state, state);
  }

  ecpu = sweeptimer.elapsedcputime(); ewall = sweeptimer.elapsedwalltime();
  pout << "\t\t\t Elapsed Sweep CPU  Time (seconds): " << setprecision(3) << ecpu << endl;
  pout << "\t\t\t Elapsed Sweep Wall Time (seconds): " << setprecision(3) << ewall << endl;

  return sweepParams.get_lowest_energy()[0];
}
}
