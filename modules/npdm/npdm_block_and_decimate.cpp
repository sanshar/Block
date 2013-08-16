/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


//FIXME includes not all needed
#include "global.h"
#include "solver.h"
#include "initblocks.h"
#include "rotationmat.h"
#include "davidson.h"
#include "linear.h"
#include "guess_wavefunction.h"
#include "density.h"
#include "davidson.h"
#include "pario.h"
#include "twopdm.h"
#include "threepdm.h"
#include "fourpdm.h"
#include "npdm_driver.h"
#include "npdm_block_and_decimate.h"

#ifndef SERIAL
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif
using namespace boost;
using namespace std;

namespace SpinAdapted{

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm::BlockAndDecimate (std::string npdm_mode, SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, 
                             const bool &useSlater, const bool& dot_with_sys, int state)
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
  //if (useSlater) {
    systemDot = SpinBlock(systemDotStart, systemDotEnd);
    //SpinBlock::store(true, systemDot.get_sites(), systemDot);
    //}
    //else
    //SpinBlock::restore(true, spindotsites, systemDot);
  SpinBlock environment, environmentDot, newEnvironment;

  int environmentDotStart, environmentDotEnd, environmentStart, environmentEnd;

  const int nexact = forward ? sweepParams.get_forward_starting_size() : sweepParams.get_backward_starting_size();

  system.addAdditionalCompOps();
//FIXME MAW change depending on forward or backward which operators are assigned to which mpi procs
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
//#ifdef SERIAL
//  const int numprocs = 1;
//#endif
//#ifndef SERIAL
//  const int numprocs = world.size();
//#endif
// >>>>>>>>> MAW
//  if (sweepParams.get_block_iter() == 0)
//    compute_twopdm_initial(solution, system, systemDot, newSystem, newEnvironment, big, numprocs, state);
//
//  compute_twopdm_sweep(solution, system, systemDot, newSystem, newEnvironment, big, numprocs, state);
//
//  if (sweepParams.get_block_iter()  == sweepParams.get_n_iters() - 1)
//    compute_twopdm_final(solution, system, systemDot, newSystem, newEnvironment, big, numprocs, state);
// <<<<<<<<<< MAW

//MAW
  int sweepPos = sweepParams.get_block_iter();
  int endPos = sweepParams.get_n_iters()-1;
  if ( npdm_mode == "twopdm" ) {
//    compute_twopdm_sweep(solution, big, state, sweepPos, endPos);
    Npdm_driver twopdm_driver(solution, big);
    twopdm_driver.compute_npdm_sweep(state, sweepPos, endPos);
  }
  else if ( npdm_mode == "threepdm" ) {
    compute_threepdm_sweep(solution, big, state, sweepPos, endPos);
  }
  else if ( npdm_mode == "fourpdm" ) {
    compute_fourpdm_sweep(solution, big, state, sweepPos, endPos);
  }
  else assert(false);
//MAW

  SaveRotationMatrix (newSystem.get_sites(), rotateMatrix, state);

  //for(int i=0;i<dmrginp.nroots();++i)
  solution[0].SaveWavefunctionInfo (big.get_stateInfo(), big.get_leftBlock()->get_sites(), state);

  newSystem.transform_operators(rotateMatrix);

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

}
