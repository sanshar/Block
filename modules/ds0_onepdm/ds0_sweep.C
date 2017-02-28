/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "density.h"
#include "ds0_sweeponepdm.h"
#include "pario.h"
#include "ds0_onepdm.h"

namespace SpinAdapted{
namespace ds0_onepdm{


void BlockAndDecimate (SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, const bool &useSlater, const bool& dot_with_sys, int state, int stateB)
{

  if (dmrginp.outputlevel() > 0) {
    mcheck("at the start of block and decimate");
    pout << "\t\t\t dot with system "<<dot_with_sys<<endl;
  }

  // figure out if we are going forward or backwards
  dmrginp.guessgenT -> start();
  bool forward = (system.get_sites() [0] == 0);
  SpinBlock systemDot;
  SpinBlock  environment;
  SpinBlock  newEnvironment;
  SpinBlock big;
  int systemDotStart, systemDotEnd;

  int systemDotSize = sweepParams.get_sys_add() - 1;

  if (forward)
  {
    systemDotStart = dmrginp.spinAdapted() ? *system.get_sites().rbegin () + 1: (*system.get_sites().rbegin ())/2 + 1 ;
    systemDotEnd = systemDotStart + systemDotSize;
  }
  else
  {
    systemDotStart = dmrginp.spinAdapted() ? system.get_sites()[0] - 1: (system.get_sites()[0])/2 - 1 ;
    systemDotEnd = systemDotStart - systemDotSize;
  }

    systemDot = SpinBlock(systemDotStart, systemDotEnd, system.get_integralIndex(), true);

  const int nexact = forward ? sweepParams.get_forward_starting_size() : sweepParams.get_backward_starting_size();

  system.addAdditionalCompOps(); 

    InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, state, stateB,
                                   sweepParams.get_sys_add(), dmrginp.direct(), system.get_integralIndex(), DISTRIBUTED_STORAGE, true, true, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);

    InitBlocks::InitNewEnvironmentBlock(environment, systemDot, newEnvironment, system, systemDot,
                                      state, stateB,
                                      sweepParams.get_sys_add(), sweepParams.get_env_add(), forward, dmrginp.direct(),
                                      sweepParams.get_onedot(), nexact, useSlater, environment.get_integralIndex(), true, true, true, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
 
  newSystem.set_loopblock(true);
  system.set_loopblock(false);
  newEnvironment.set_loopblock(false);
  InitBlocks::InitBigBlock(newSystem, newEnvironment, big); 

  std::vector<Wavefunction> solution(2);

  DiagonalMatrix e;

  GuessWave::guess_wavefunctions(solution[0], e, big, sweepParams.get_guesstype(), true, stateB, true, 0.0, true);
  GuessWave::guess_wavefunctions(solution[1], e, big, sweepParams.get_guesstype(), true, state, true, 0.0, false);

#ifndef SERIAL
  mpi::communicator world;
  mpi::broadcast(world, solution, 0);
#endif
  std::vector<Matrix> rotateMatrix;
  std::vector<Matrix> rotateMatrixB;

    DensityMatrix tracedMatrix(newSystem.get_braStateInfo());
    tracedMatrix.allocate(newSystem.get_braStateInfo());
    tracedMatrix.makedensitymatrix(std::vector<Wavefunction>(1,solution[1]), big, std::vector<double>(1,1.0), 0.0, 0.0, false);
    rotateMatrix.clear();

    DensityMatrix tracedMatrixB(newSystem.get_ketStateInfo());
    tracedMatrixB.allocate(newSystem.get_ketStateInfo());
    tracedMatrixB.makedensitymatrix(std::vector<Wavefunction>(1,solution[0]), big, std::vector<double>(1,1.0), 0.0, 0.0, false);
    rotateMatrixB.clear();

    if (!mpigetrank()){
     double error = makeRotateMatrix(tracedMatrix, rotateMatrix, sweepParams.get_keep_states(), sweepParams.get_keep_qstates());
     double error2 = makeRotateMatrix(tracedMatrixB, rotateMatrixB, sweepParams.get_keep_states(), sweepParams.get_keep_qstates());
    }

#ifndef SERIAL
    mpi::broadcast(world,rotateMatrix,0);
    if(state!=stateB)
    mpi::broadcast(world,rotateMatrixB,0);
#endif

#ifdef SERIAL
  const int numprocs = 1;
#endif
#ifndef SERIAL
  const int numprocs = world.size();
#endif



  Matrix onepdm; 
  load_onepdm_binary(onepdm, state ,stateB); 

//  if (sweepParams.get_block_iter() == 0) {
    pout << "compute 2_0"<<endl;
    compute_one_pdm_2_0(solution[1], solution[0], big, onepdm);
//}

    pout << "compute 1_1"<<endl;
    compute_one_pdm_1_1(solution[1], solution[0], big, onepdm);

  if (sweepParams.get_block_iter()  == sweepParams.get_n_iters() - 1) {
    pout << "compute 0_2"<<endl;
    compute_one_pdm_0_2(solution[1], solution[0], big, onepdm);
}
  accumulate_onepdm(onepdm); 
  save_onepdm_binary(onepdm, state, stateB); 

  SaveRotationMatrix (newSystem.get_sites(), rotateMatrix, state);
  SaveRotationMatrix (newSystem.get_sites(), rotateMatrixB, stateB);

  solution[1].SaveWavefunctionInfo (big.get_braStateInfo(), big.get_leftBlock()->get_sites(), state);
  solution[0].SaveWavefunctionInfo (big.get_ketStateInfo(), big.get_leftBlock()->get_sites(), stateB);

  newSystem.transform_operators(rotateMatrix,rotateMatrixB);

}

double do_one(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize,int state, int stateB)
{
  pout.precision(12);

  SpinBlock system;

  Matrix onepdm(2*dmrginp.last_site(), 2*dmrginp.last_site());onepdm=0.0;
  save_onepdm_binary(onepdm, state ,stateB);

  sweepParams.set_sweep_parameters();
  // a new renormalisation sweep routine
  pout << ((forward) ? "\t\t\t Starting renormalisation sweep in forwards direction" : "\t\t\t Starting renormalisation sweep in backwards direction") << endl;
  pout << "\t\t\t ============================================================================ " << endl;
  
    InitBlocks::InitStartingBlock( system, forward, state, stateB, sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), restartSize, restart, warmUp,0); 
    pout << "\t\t\t Starting block is :: " << endl << system << endl;

    sweepParams.set_block_iter() = 0;   

//    SpinBlock::store (forward, system.get_sites(), system, state, stateB ); 
 
   sweepParams.savestate(forward, system.get_sites().size());
   bool dot_with_sys = true;

  sweepParams.set_guesstype() = TRANSPOSE;
  // Loop over all block sites
  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); ) {
    pout << "\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
    pout << "\t\t\t ----------------------------" << endl;
    if (forward) pout << "\t\t\t Current direction is :: Forwards " << endl;
    else pout << "\t\t\t Current direction is :: Backwards " << endl;

    if (!warmUp && sweepParams.get_block_iter() != 0)
            sweepParams.set_guesstype() = TRANSFORM;
    else
      sweepParams.set_guesstype() = TRANSPOSE;

    pout << "\t\t\t Blocking and Decimating " << endl;

  SpinBlock newSystem;

  BlockAndDecimate (sweepParams, system, newSystem, warmUp, dot_with_sys, state,stateB);

  system = newSystem;
  pout << system<<endl;

//  SpinBlock::store (forward, system.get_sites(), system, state, stateB);

  pout << "\t\t\t saving state " << system.get_sites().size() << endl;
  ++sweepParams.set_block_iter();
}

  load_onepdm_binary(onepdm, state ,stateB);
  accumulate_onepdm(onepdm);
  save_onepdm_spatial_text(onepdm, state, stateB);
  save_onepdm_text(onepdm, state, stateB);
  save_onepdm_spatial_binary(onepdm, state, stateB);
}
}
}










	






















