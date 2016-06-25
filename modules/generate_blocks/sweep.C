/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "operatorfunctions.h"
#include "sweepgenblock.h"
#include "guess_wavefunction.h"
#include "density.h"
#include "davidson.h"
#include "pario.h"
#ifdef USE_BTAS
#include "overlaptensor.h"
#endif
#include "sweep.h"

namespace SpinAdapted{
void SweepGenblock::BlockAndDecimate (SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, const bool &useSlater, const bool& dot_with_sys, int stateA, int stateB)
{
  if (dmrginp.outputlevel() > 0) 
    mcheck("at the start of block and decimate");
  p1out << "\t\t\t Performing Blocking"<<endl;
  dmrginp.guessgenT -> start();
  // figure out if we are going forward or backwards  
  bool forward = (system.get_sites() [0] == 0);
  SpinBlock systemDot;
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
  dmrginp.sysdotmake->start();
  systemDot = SpinBlock(systemDotStart, systemDotEnd, system.get_integralIndex(), stateA==stateB);
  dmrginp.sysdotmake->stop();

  const int nexact = forward ? sweepParams.get_forward_starting_size() : sweepParams.get_backward_starting_size();

  dmrginp.guessgenT -> stop();
  dmrginp.datatransfer -> start();
  system.addAdditionalCompOps();
  dmrginp.datatransfer -> stop();
  dmrginp.initnewsystem->start();
  InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, stateA, stateB, sweepParams.get_sys_add(), dmrginp.direct(), system.get_integralIndex(), DISTRIBUTED_STORAGE, dot_with_sys, true);
  dmrginp.initnewsystem->stop();

  pout << "\t\t\t System  Block"<<newSystem;
  newSystem.printOperatorSummary();

  std::vector<Matrix> leftrotateMatrix, rightrotateMatrix;

  LoadRotationMatrix (newSystem.get_sites(), leftrotateMatrix, stateA);
  LoadRotationMatrix (newSystem.get_sites(), rightrotateMatrix, stateB);

#ifndef SERIAL
  mpi::communicator world;
  broadcast(world, leftrotateMatrix, 0);
  broadcast(world, rightrotateMatrix, 0);
#endif

  p1out <<"\t\t\t Performing Renormalization "<<endl<<endl;

  dmrginp.operrotT->start();
  if (stateB == stateA)
    newSystem.transform_operators(leftrotateMatrix);
  else
    newSystem.transform_operators(leftrotateMatrix, rightrotateMatrix);
  dmrginp.operrotT->stop();


  if (dmrginp.outputlevel() > 0) 
    //mcheck("after rotation and transformation of block");
  p2out <<newSystem<<endl;
  newSystem.printOperatorSummary();
  //mcheck("After renorm transform");

  p2out << *dmrginp.guessgenT<<" "<<*dmrginp.multiplierT<<" "<<*dmrginp.operrotT<< "  "<<globaltimer.totalwalltime()<<" timer "<<endl;
  p2out << *dmrginp.makeopsT<<"  "<<*dmrginp.initnewsystem<<"  "<<*dmrginp.sysdotmake<<"  "<<*dmrginp.buildcsfops<<" makeops "<<endl;
  p2out << *dmrginp.datatransfer<<" datatransfer "<<endl;
  p2out <<"oneindexopmult   twoindexopmult   Hc  couplingcoeff"<<endl;  
  p2out << *dmrginp.oneelecT<<" "<<*dmrginp.twoelecT<<" "<<*dmrginp.hmultiply<<" "<<*dmrginp.couplingcoeff<<" hmult"<<endl;
  p2out << *dmrginp.buildsumblock<<" "<<*dmrginp.buildblockops<<" build block"<<endl;
  p2out << *dmrginp.blockintegrals<<"  "<<*dmrginp.blocksites<<"  "<<*dmrginp.statetensorproduct<<"  "<<*dmrginp.statecollectquanta<<"  "<<*dmrginp.buildsumblock<<" "<<*dmrginp.builditeratorsT<<"  "<<*dmrginp.diskio<<" build sum block"<<endl;
  p2out << "addnoise  S_0_opxop  S_1_opxop   S_2_opxop"<<endl;
  p3out << *dmrginp.addnoise<<" "<<*dmrginp.s0time<<" "<<*dmrginp.s1time<<" "<<*dmrginp.s2time<<endl;

}

double SweepGenblock::do_one(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize, int stateA, int stateB)
{
  Timer sweeptimer;
  int integralIndex = 0;

  SpinBlock system;
  const int nroots = dmrginp.nroots();
  std::vector<double> finalEnergy(nroots,0.);
  std::vector<double> finalEnergy_spins(nroots,0.);
  double finalError = 0.;

  sweepParams.set_sweep_parameters();
  // a new renormalisation sweep routine
  pout << ((forward) ? "\t\t\t Starting renormalisation sweep in forwards direction" : "\t\t\t Starting renormalisation sweep in backwards direction") << endl;
  pout << "\t\t\t ============================================================================ " << endl;
  
  InitBlocks::InitStartingBlock (system,forward, stateA, stateB, sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), restartSize, restart, warmUp, integralIndex);
  if(!restart)
    sweepParams.set_block_iter() = 0;

  p2out << "\t\t\t Starting block is :: " << endl << system << endl;
  //if (!restart) 
  SpinBlock::store (forward, system.get_sites(), system, stateA, stateB); // if restart, just restoring an existing block --
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
      pout << "\n\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
      pout << "\t\t\t ----------------------------" << endl;
      if (forward)
	{ p1out << "\t\t\t Current direction is :: Forwards " << endl; }
      else
	{ p1out << "\t\t\t Current direction is :: Backwards " << endl; }
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
      
      p1out << "\t\t\t Blocking and Decimating " << endl;
	  
      SpinBlock newSystem;

      BlockAndDecimate (sweepParams, system, newSystem, warmUp, dot_with_sys, stateA, stateB);

      
      system = newSystem;

      //system size is going to be less than environment size
      if (forward && system.get_complementary_sites()[0] >= dmrginp.last_site()/2)
	dot_with_sys = false;
      if (!forward && system.get_sites()[0]-1 < dmrginp.last_site()/2)
	dot_with_sys = false;

      SpinBlock::store (forward, system.get_sites(), system, stateA, stateB);	 	

      p1out << "\t\t\t saving state " << system.get_sites().size() << endl;
      ++sweepParams.set_block_iter();
      //if (sweepParams.get_onedot())
      //pout << "\t\t\tUsing one dot algorithm!!"<<endl; 
      sweepParams.savestate(forward, system.get_sites().size());
    }
  pout << "\t\t\t Finished Generate-Blocks Sweep. " << endl;
  pout << "\t\t\t ============================================================================ " << endl;

  // update the static number of iterations

  ++sweepParams.set_sweep_iter();

  ecpu = sweeptimer.elapsedcputime(); ewall = sweeptimer.elapsedwalltime();
  pout << "\t\t\t Elapsed Sweep CPU  Time (seconds): " << setprecision(3) << ecpu << endl;
  pout << "\t\t\t Elapsed Sweep Wall Time (seconds): " << setprecision(3) << ewall << endl;

  return finalEnergy[0];
}


void SweepGenblock::do_one(SweepParams &sweepParams, const bool &forward, int stateA, int stateB)
{
  Timer sweeptimer;
  int integralIndex = 0;
  SpinBlock system;

  sweepParams.set_sweep_parameters();
  // a new renormalisation sweep routine
  pout << ((forward) ? "\t\t\t Starting renormalisation sweep in forwards direction" : "\t\t\t Starting renormalisation sweep in backwards direction") << endl;
  pout << "\t\t\t ============================================================================ " << endl;
  
  InitBlocks::InitStartingBlock (system,forward, stateA, stateB, sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), 0, false, false, integralIndex);

  sweepParams.set_block_iter() = 0;

  p2out << "\t\t\t Starting block is :: " << endl << system << endl;

  bool dot_with_sys = true;

  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); )
    {
      pout << "\n\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
      pout << "\t\t\t ----------------------------" << endl;
      if (forward)
	{ p1out << "\t\t\t Current direction is :: Forwards " << endl; }
      else
	{ p1out << "\t\t\t Current direction is :: Backwards " << endl; }

  
      if (dmrginp.no_transform())
	      sweepParams.set_guesstype() = BASIC;
      else if ( sweepParams.get_block_iter() != 0) 
  	    sweepParams.set_guesstype() = TRANSFORM;
      else if ( sweepParams.get_block_iter() == 0 )
        sweepParams.set_guesstype() = TRANSPOSE;
      else
        sweepParams.set_guesstype() = BASIC;
      
      p1out << "\t\t\t Blocking and Decimating " << endl;
	  
      SpinBlock newSystem;

      BlockAndDecimate (sweepParams, system, newSystem, false, dot_with_sys, stateA, stateB);

      system = newSystem;

      SpinBlock::store(forward, system.get_sites(), system, stateA, stateB);

      //system size is going to be less than environment size
      if (forward && system.get_complementary_sites()[0] >= dmrginp.last_site()/2)
	dot_with_sys = false;
      if (!forward && system.get_sites()[0]-1 < dmrginp.last_site()/2)
	dot_with_sys = false;

      ++sweepParams.set_block_iter();
    }
  pout << "\t\t\t Finished Generate-Blocks Sweep. " << endl;
  pout << "\t\t\t ============================================================================ " << endl;

  // update the static number of iterations

  ++sweepParams.set_sweep_iter();

  ecpu = sweeptimer.elapsedcputime(); ewall = sweeptimer.elapsedwalltime();
  pout << "\t\t\t Elapsed Sweep CPU  Time (seconds): " << setprecision(3) << ecpu << endl;
  pout << "\t\t\t Elapsed Sweep Wall Time (seconds): " << setprecision(3) << ewall << endl;

}


}
