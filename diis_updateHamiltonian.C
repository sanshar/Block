#include "diis.h"
#include "global.h"
#include "solver.h"
#include "initblocks.h"
#include "MatrixBLAS.h"
#include "density.h"
#include "rotationmat.h"
#include "davidson.h"
#include "linear.h"
#include "guess_wavefunction.h"
#include <include/sortutils.h>

#ifndef SERIAL
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif

using namespace boost;
using namespace std;



void SpinAdapted::DIIS::updateHamiltonian(SweepParams &sweepParams, bool forward)
{
  bool restart = false;
  int restartSize = 0;
  bool warmUp = false;

  SpinBlock system;
  const int nroots = dmrginp.nroots(sweepParams.get_sweep_iter());
  std::vector<double> finalEnergy(nroots,0.);
  std::vector<double> finalEnergy_spins(nroots,0.);
  double finalError = 0.;


  sweepParams.set_sweep_parameters();

  // a new renormalisation sweep routine
  pout << ((forward) ? "\t\t\t Starting renormalisation sweep in forwards direction" : "\t\t\t Starting renormalisation sweep in backwards direction") << endl;
  pout << "\t\t\t ============================================================================ " << endl;
  mcheck("");
  InitBlocks::InitStartingBlock (system,forward, sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), restartSize, restart, warmUp);
  if(!restart)
    sweepParams.set_block_iter() = 0;
 
  pout << "\t\t\t Starting block is :: " << endl << system << endl;
  if (!restart) 
    SpinBlock::store (forward, system.get_sites(), system); // if restart, just restoring an existing block --
  sweepParams.savestate(forward, system.get_sites().size());
  bool dot_with_sys = true;
  vector<int> syssites;
  {
    syssites = system.get_sites();
  }
  mcheck("at the very start of sweep");

  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); )
    {
      
      pout << "\t\t\t Sweep Iteration :: " << sweepParams.get_block_iter() << endl;
      pout << "\t\t\t ----------------------------" << endl;
      if (forward)
	pout << "\t\t\t Current direction is :: Forwards " << endl;
      else
	pout << "\t\t\t Current direction is :: Backwards " << endl;

      //if (SHOW_MORE) pout << "system block" << endl << system << endl;
  
      if (dmrginp.no_transform())
	      sweepParams.set_guesstype() = BASIC;
      else if (!warmUp && sweepParams.get_block_iter() != 0) 
  	    sweepParams.set_guesstype() = TRANSFORM;
      else if (!warmUp && sweepParams.get_block_iter() == 0 && 
                ((dmrginp.algorithm_method() == TWODOT_TO_ONEDOT && dmrginp.twodot_to_onedot_iter() != sweepParams.get_sweep_iter()) ||
                  dmrginp.algorithm_method() != TWODOT_TO_ONEDOT))
        sweepParams.set_guesstype() = TRANSPOSE;
      else
        sweepParams.set_guesstype() = BASIC;
      
      pout << "\t\t\t Blocking and Decimating " << endl;
	  
      SpinBlock newSystem;

      BlockAndDecimate (sweepParams, system, newSystem, warmUp, dot_with_sys);

      
      system = newSystem;
      pout << system<<endl;

      SpinBlock::store (forward, system.get_sites(), system);	 	
      syssites = system.get_sites();
      pout << "\t\t\t saving state " << syssites.size() << endl;
      ++sweepParams.set_block_iter();
      sweepParams.savestate(forward, syssites.size());
      mcheck("at the end of sweep iteration");
    }



  return ;
}

void SpinAdapted::DIIS::BlockAndDecimate (SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, const bool &useSlater, const bool& dot_with_sys)
{
  mcheck("at the start of block and decimate");
  pout << "dot with system "<<dot_with_sys<<endl;
  // figure out if we are going forward or backwards
  dmrginp.guessgenT.start();
  bool forward = (system.get_sites() [0] == 0);
  SpinBlock systemDot;
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
  int environmentDotSize = sweepParams.get_env_add() - 1;
  if (forward)
  {
    environmentDotStart = systemDotEnd + 1;
    environmentDotEnd = environmentDotStart + environmentDotSize;
  }
  else
  {
    environmentDotStart = systemDotEnd - 1;
    environmentDotEnd = environmentDotStart - environmentDotSize;
  }
  vector<int> envdotsites(2); 
  envdotsites[0] = environmentDotStart;
  envdotsites[1] = environmentDotEnd;

  environmentDot = SpinBlock(environmentDotStart, environmentDotEnd);
  
  const int nexact = forward ? sweepParams.get_forward_starting_size() : sweepParams.get_backward_starting_size();

  InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, sweepParams.get_sys_add(), dmrginp.direct(), 
				 DISTRIBUTED_STORAGE, dot_with_sys, true);

  InitBlocks::InitNewEnvironmentBlock(environment, systemDot, newEnvironment, system, systemDot,
				      sweepParams.get_sys_add(), sweepParams.get_env_add(), forward, dmrginp.direct(),
				      sweepParams.get_onedot(), nexact, useSlater, !dot_with_sys, true, dot_with_sys);
  
  SpinBlock big;

  newSystem.set_loopblock(true);
  system.set_loopblock(false);
  newEnvironment.set_loopblock(false);
  if (!sweepParams.get_onedot())
    environment.set_loopblock(false);
  InitBlocks::InitBigBlock(newSystem, newEnvironment, big); 
  

  
  std::vector<Wavefunction> waves(1);

  
  DiagonalMatrix e; e.ReSize(big.get_stateInfo().totalStates); e=0;


  if(buildup && currentIndex == 0)
    GuessWave::guess_wavefunctions(waves, e, big, sweepParams.get_guesstype(), true, true);
  else {
    if(forward)
      LoadInterpolatedWavefunction (waves[0], sweepParams.get_block_iter());
    else
      LoadInterpolatedWavefunction (waves[0], sweepParams.get_n_iters() - sweepParams.get_block_iter()-1);
  }

  Normalise(waves[0]);

  StateInfo newenvState, envstate = newEnvironment.get_stateInfo(), sysdotstate = systemDot.get_stateInfo(), sysstate = system.get_stateInfo();
  TensorProduct( envstate, sysdotstate, newenvState, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT); 
  newenvState.CollectQuanta();
  StateInfo bigstate;
  TensorProduct(newenvState, sysstate, bigstate, PARTICLE_SPIN_NUMBER_CONSTRAINT); 
  Wavefunction guesswf;
  if(!forward && !(buildup&& currentIndex==0)) {
    guesswf = waves[0];
    //we have to transpose the wavefunction
    GuessWave::onedot_transpose_wavefunction(bigstate, big.get_stateInfo(),
					     guesswf, waves[0]);
  }
  cout << waves[0]<<endl;
  //this can be deleted later
  multiply_h h_multiply(big, true);
  std::vector<Wavefunction> sigmaV(1,waves[0]);
  sigmaV[0].Clear();
  h_multiply(waves[0], sigmaV[0]);
  double E = DotProduct(waves[0], sigmaV[0]);
  ScaleAdd(-E, waves[0], sigmaV[0]);
  double error = DotProduct(sigmaV[0], sigmaV[0]);
  cout << "Energy :: "<<E<<" Error :: "<<error<<endl;
  //*************

  waves[0].SaveWavefunctionInfo(big.get_stateInfo(), big.get_leftBlock()->get_sites(), 0);

  transformBlock(big, waves, sweepParams);

}

