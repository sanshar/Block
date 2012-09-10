/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is free software: you can redistribute it and/or modify         
it under the terms of the GNU General Public License as published by         
the Free Software Foundation, either version 3 of the License, or            
(at your option) any later version.                                          
                                                                             
This program is distributed in the hope that it will be useful,              
but WITHOUT ANY WARRANTY; without even the implied warranty of               
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
GNU General Public License for more details.                                 
                                                                             
You should have received a copy of the GNU General Public License            
along with this program.  If not, see <http://www.gnu.org/licenses/>.        

This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "sweepgenblock.h"
#include "guess_wavefunction.h"
#include "density.h"
#include "davidson.h"

namespace SpinAdapted{
void SweepGenblock::BlockAndDecimate (SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, const bool &useSlater, const bool& dot_with_sys, int state)
{
  if (dmrginp.outputlevel() != 0) 
    mcheck("at the start of block and decimate");
  // figure out if we are going forward or backwards
  pout << "\t\t\t Performing Blocking"<<endl;
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

  const int nexact = forward ? sweepParams.get_forward_starting_size() : sweepParams.get_backward_starting_size();

  system.addAdditionalCompOps();
  InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, sweepParams.get_sys_add(), dmrginp.direct(), DISTRIBUTED_STORAGE, dot_with_sys, true);


  pout << "\t\t\t System  Block"<<newSystem;
  if (dmrginp.outputlevel() != 0)
    newSystem.printOperatorSummary();

  std::vector<Matrix> rotateMatrix;


  if (!dmrginp.get_fullrestart()) {
    //this should be done when we actually have wavefunctions stored, otherwise not!!
    SpinBlock environment, environmentDot, newEnvironment;
    int environmentDotStart, environmentDotEnd, environmentStart, environmentEnd;
    InitBlocks::InitNewEnvironmentBlock(environment, systemDot, newEnvironment, system, systemDot,
					sweepParams.get_sys_add(), sweepParams.get_env_add(), forward, dmrginp.direct(),
					sweepParams.get_onedot(), nexact, useSlater, true, true, true);
    SpinBlock big;
    InitBlocks::InitBigBlock(newSystem, newEnvironment, big); 
    DiagonalMatrix e;
    std::vector<Wavefunction> solution(1);
    
    GuessWave::guess_wavefunctions(solution[0], e, big, sweepParams.get_guesstype(), true, state, true, 0.0); 
    solution[0].SaveWavefunctionInfo (big.get_stateInfo(), big.get_leftBlock()->get_sites(), state);


    DensityMatrix tracedMatrix;
    tracedMatrix.allocate(newSystem.get_stateInfo());
    tracedMatrix.makedensitymatrix(solution, big, std::vector<double>(1, 1.0), 0.0, 0.0, false);
    rotateMatrix.clear();
    if (!mpigetrank())
      double error = newSystem.makeRotateMatrix(tracedMatrix, rotateMatrix, sweepParams.get_keep_states(), sweepParams.get_keep_qstates());
    
  }
  else
    LoadRotationMatrix (newSystem.get_sites(), rotateMatrix, state);

#ifndef SERIAL
  mpi::communicator world;
  broadcast(world, rotateMatrix, 0);
#endif

  if (!dmrginp.get_fullrestart())
    SaveRotationMatrix (newSystem.get_sites(), rotateMatrix, state);

  pout <<"\t\t\t Performing Renormalization "<<endl<<endl;
  newSystem.transform_operators(rotateMatrix);
  if (dmrginp.outputlevel() != 0) 
    mcheck("after rotation and transformation of block");
  if (dmrginp.outputlevel() != 0) 
    pout <<newSystem<<endl;
  if (dmrginp.outputlevel() != 0)
    newSystem.printOperatorSummary();
  //mcheck("After renorm transform");
}

double SweepGenblock::do_one(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize, int state)
{

  SpinBlock system;
  const int nroots = dmrginp.nroots();
  std::vector<double> finalEnergy(nroots,0.);
  std::vector<double> finalEnergy_spins(nroots,0.);
  double finalError = 0.;

  sweepParams.set_sweep_parameters();
  // a new renormalisation sweep routine
  pout << ((forward) ? "\t\t\t Starting renormalisation sweep in forwards direction" : "\t\t\t Starting renormalisation sweep in backwards direction") << endl;
  pout << "\t\t\t ============================================================================ " << endl;
  
  InitBlocks::InitStartingBlock (system,forward, sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), restartSize, restart, warmUp);
  if(!restart)
    sweepParams.set_block_iter() = 0;

  if (dmrginp.outputlevel() != 0) 
    pout << "\t\t\t Starting block is :: " << endl << system << endl;
  //if (!restart) 
    SpinBlock::store (forward, system.get_sites(), system); // if restart, just restoring an existing block --
  sweepParams.savestate(forward, system.get_sites().size());
  bool dot_with_sys = true;
  if (restart)
  {
    if (forward && system.get_complementary_sites()[0] >= dmrginp.last_site()/2)
      dot_with_sys = false;
    if (!forward && system.get_sites()[0]-1 < dmrginp.last_site()/2)
      dot_with_sys = false;
  }

  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); )
    {
      pout << "\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
      pout << "\t\t\t ----------------------------" << endl;
      if (dmrginp.outputlevel() != 0) {
      if (forward)
	pout << "\t\t\t Current direction is :: Forwards " << endl;
      else
	pout << "\t\t\t Current direction is :: Backwards " << endl;
      }
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
      
      if (dmrginp.outputlevel() != 0) 
	pout << "\t\t\t Blocking and Decimating " << endl;
	  
      SpinBlock newSystem;

      BlockAndDecimate (sweepParams, system, newSystem, warmUp, dot_with_sys, state);

      
      system = newSystem;

      //system size is going to be less than environment size
      if (forward && system.get_complementary_sites()[0] >= dmrginp.last_site()/2)
	dot_with_sys = false;
      if (!forward && system.get_sites()[0]-1 < dmrginp.last_site()/2)
	dot_with_sys = false;

      SpinBlock::store (forward, system.get_sites(), system);	 	

      if (dmrginp.outputlevel() != 0) 
	pout << "\t\t\t saving state " << system.get_sites().size() << endl;
      ++sweepParams.set_block_iter();
      //if (sweepParams.get_onedot())
      //pout << "\t\t\tUsing one dot algorithm!!"<<endl; 
      sweepParams.savestate(forward, system.get_sites().size());
    }
  pout << "\t\t\t Finished Generate-Blocks Sweep. " << endl;
  pout << "\t\t\t ============================================================================ " << endl;

  // update the static number of iterations

  ++sweepParams.set_sweep_iter();

  return finalEnergy[0];
}
}
