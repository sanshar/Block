/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "guess_wavefunction.h"
#include "sweepResponse.h"
#include "global.h"
#include "solver.h"
#include "initblocks.h"
#include "MatrixBLAS.h"
#include <boost/format.hpp>
#ifndef SERIAL
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif
#include "rotationmat.h"
#include "density.h"
#include "pario.h"
#include "davidson.h"
#include "sweep.h"
#include "initblocks.h"

using namespace boost;
using namespace std;


void SpinAdapted::SweepResponse::BlockAndDecimate (SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, const bool &useSlater, const bool& dot_with_sys, int targetState, int correctionVector, int baseState)
{
  if (dmrginp.outputlevel() > 0) {
    mcheck("at the start of block and decimate");
    pout << "\t\t\t dot with system "<<dot_with_sys<<endl;
  }
  pout <<endl<< "\t\t\t Performing Blocking"<<endl;
  // figure out if we are going forward or backwards
  dmrginp.guessgenT -> start();

  bool forward = (system.get_sites() [0] == 0);
  SpinBlock systemDot, environmentDot;
  int systemDotStart, systemDotEnd, environmentDotStart, environmentDotEnd;
  int systemDotSize = sweepParams.get_sys_add() - 1;
  int environmentDotSize = sweepParams.get_env_add() - 1;
  if (forward)
  {
    systemDotStart = dmrginp.spinAdapted() ? 
      *system.get_sites().rbegin () + 1 : (*system.get_sites().rbegin ())/2 + 1 ;
    systemDotEnd = systemDotStart + systemDotSize;
    environmentDotStart = systemDotEnd + 1;
    environmentDotEnd = environmentDotStart + environmentDotSize;
  }
  else
  {
    systemDotStart = dmrginp.spinAdapted() ? 
      system.get_sites()[0] - 1 : (system.get_sites()[0])/2 - 1 ;
    systemDotEnd = systemDotStart - systemDotSize;
    environmentDotStart = systemDotEnd - 1;
    environmentDotEnd = environmentDotStart - environmentDotSize;
  }
  systemDot = SpinBlock(systemDotStart, systemDotEnd, system.get_integralIndex(), true);
  environmentDot = SpinBlock(environmentDotStart, environmentDotEnd, system.get_integralIndex(), true);
  SpinBlock environment, newEnvironment;

  SpinBlock big;  // new_sys = sys+sys_dot; new_env = env+env_dot; big = new_sys + new_env then renormalize to find new_sys(new)
  Sweep::makeSystemEnvironmentBigBlocks(system, systemDot, newSystem, environment, environmentDot, 
					newEnvironment, big, sweepParams, dot_with_sys, useSlater, system.get_integralIndex(),  
					targetState, targetState);


  //analyse_operator_distribution(big);
  dmrginp.guessgenT -> stop();
  dmrginp.multiplierT -> start();
  std::vector<Matrix> rotatematrix;

  if (dmrginp.outputlevel() > 0)
    mcheck(""); 
  if (dmrginp.outputlevel() == 0) {
    if (!dot_with_sys && sweepParams.get_onedot()) pout << "\t\t\t System  Block"<<system;    
    else pout << "\t\t\t System  Block"<<newSystem;
    pout << "\t\t\t Environment Block"<<newEnvironment<<endl;
    pout << "\t\t\t Solving wavefunction "<<endl;
  }

  std::vector<Wavefunction> lowerStates;


  //make the baseState
  int originalOutputlevel = dmrginp.outputlevel();
  dmrginp.setOutputlevel() = -1;
  lowerStates.resize(2); //only make the 
  DiagonalMatrix e;
  guessWaveTypes guesstype = sweepParams.get_block_iter() == 0 ? TRANSPOSE : TRANSFORM;

  std::vector<int> systemsites;
  if (dmrginp.spinAdapted())
    systemsites = system.get_sites();
  else {
    if (forward) {
      systemsites.push_back(0); systemsites.push_back(*system.get_sites().rbegin()/2);
    }
    else {
      systemsites.push_back(system.get_sites()[0]/2); systemsites.push_back(*system.get_sites().rbegin()/2);
    }
  }

  //<target|O|baseState>
  {  
    //now one needs to make |phi_0> = O|psi_0> so that the |phi_0> has the same dimensions as our target state
    SpinBlock overlapBig;
    SpinBlock overlapsystem, overlapenvironment, overlapnewsystem, overlapnewenvironment;
    Sweep::makeSystemEnvironmentBigOverlapBlocks(systemsites, systemDot, environmentDot,
						 overlapsystem, overlapenvironment, overlapnewsystem, 
						 overlapnewenvironment, overlapBig, sweepParams, 
						 dot_with_sys, useSlater, system.get_integralIndex(), targetState, baseState);
    GuessWave::guess_wavefunctions(lowerStates[0], e, overlapBig, guesstype, 
				   sweepParams.get_onedot(), baseState, dot_with_sys, 0.0);
    
    if (mpigetrank() == 0) {
      Wavefunction temp; 
      temp.set_onedot( sweepParams.get_onedot());
      temp.AllowQuantaFor(*overlapBig.get_braStateInfo().leftStateInfo,*overlapBig.get_braStateInfo().rightStateInfo,
			  dmrginp.effective_molecule_quantum_vec());
      temp.Clear();
      overlapBig.multiplyOverlap(lowerStates[0], &temp, MAX_THRD);
      lowerStates[0] = temp;
    }
    else 
      lowerStates[0].CleanUp();

    overlapsystem.clear(); overlapenvironment.clear(); 
    overlapnewsystem.clear(); overlapnewenvironment.clear();  
  }


  //<target|O|cV>
  {  
    //now one needs to make |phi_0> = O|psi_0> so that the |phi_0> has the same dimensions as our target state
    SpinBlock overlapBig;
    SpinBlock overlapsystem, overlapenvironment, overlapnewsystem, overlapnewenvironment;
    Sweep::makeSystemEnvironmentBigOverlapBlocks(systemsites, systemDot, environmentDot,
						 overlapsystem, overlapenvironment, overlapnewsystem, 
						 overlapnewenvironment, overlapBig, sweepParams, 
						 dot_with_sys, useSlater, system.get_integralIndex(), targetState, correctionVector);
    GuessWave::guess_wavefunctions(lowerStates[1], e, overlapBig, guesstype, 
				   sweepParams.get_onedot(), correctionVector, dot_with_sys, 0.0);


    if (mpigetrank() == 0) {
      Wavefunction temp; 
      temp.set_onedot( sweepParams.get_onedot());
      temp.AllowQuantaFor(*overlapBig.get_braStateInfo().leftStateInfo,*overlapBig.get_braStateInfo().rightStateInfo,
			  dmrginp.effective_molecule_quantum_vec());
      temp.Clear();
      overlapBig.multiplyOverlap(lowerStates[1], &temp, MAX_THRD);
      lowerStates[1] = temp;
    }
    else 
      lowerStates[1].CleanUp();

    overlapsystem.clear(); overlapenvironment.clear(); 
    overlapnewsystem.clear(); overlapnewenvironment.clear();  
  
  }
  dmrginp.setOutputlevel() = originalOutputlevel;

  newSystem.RenormaliseFrom (sweepParams.set_lowest_energy(), sweepParams.set_lowest_energy_spins(),
			     sweepParams.set_lowest_error(), rotatematrix, 
			     sweepParams.get_keep_states(), sweepParams.get_keep_qstates(), 
			     sweepParams.get_davidson_tol(), big, sweepParams.get_guesstype(), 
			     sweepParams.get_noise(), sweepParams.get_additional_noise(), 
			     sweepParams.get_onedot(), system, systemDot, environment, 
			     dot_with_sys, useSlater, sweepParams.get_sweep_iter(), targetState, 
			     lowerStates, correctionVector);



  if (dmrginp.outputlevel() > 0)
    mcheck("");
  environment.clear();
  newEnvironment.clear();

  pout <<"\t\t\t Performing Renormalization "<<endl;
  pout << "\t\t\t Total discarded weight "<<sweepParams.get_lowest_error()<<endl<<endl;


  dmrginp.multiplierT -> stop();
  dmrginp.operrotT -> start();
  newSystem.transform_operators(rotatematrix);
  dmrginp.operrotT -> stop();

  //<target|O|cV>
  //save the updated overlap spinblock
  dmrginp.setOutputlevel() = -1; 
  {
    int originalOutputlevel = dmrginp.outputlevel();
    //dmrginp.setOutputlevel() = -1;

    SpinBlock overlapBig;
    SpinBlock overlapsystem, overlapenvironment, overlapnewsystem, overlapnewenvironment;
    SpinBlock overlapsystemDot(systemDotStart, systemDotEnd, newSystem.get_integralIndex(), true);
    SpinBlock overlapenvironmentDot(environmentDotStart, environmentDotEnd, newSystem.get_integralIndex(), true);
    guessWaveTypes guesstype = sweepParams.get_block_iter() == 0 ? TRANSPOSE : TRANSFORM;
      
    DiagonalMatrix e;
    Sweep::makeSystemEnvironmentBigOverlapBlocks(systemsites, overlapsystemDot, overlapenvironmentDot,
						 overlapsystem, overlapnewsystem, overlapenvironment, overlapnewenvironment,
						 overlapBig, sweepParams, true, useSlater, newSystem.get_integralIndex(), 
						 targetState, correctionVector);

    Wavefunction iwave;
    GuessWave::guess_wavefunctions(iwave, e, overlapBig, guesstype, sweepParams.get_onedot(), correctionVector, true, 0.0);
    std::vector<Matrix> ketrotatematrix;
    DensityMatrix tracedMatrix;
    tracedMatrix.allocate(overlapnewsystem.get_ketStateInfo());
    operatorfunctions::MultiplyProduct(iwave, Transpose(const_cast<Wavefunction&> (iwave)), tracedMatrix, 1.0);
    int largeNumber = 1000000;
    if (!mpigetrank())
      double error = makeRotateMatrix(tracedMatrix, ketrotatematrix, largeNumber, sweepParams.get_keep_qstates());
    
#ifndef SERIAL
    mpi::communicator world;
    broadcast(world, ketrotatematrix, 0);
#endif
    
    iwave.SaveWavefunctionInfo (overlapBig.get_ketStateInfo(), overlapBig.get_leftBlock()->get_sites(), correctionVector);
    SaveRotationMatrix (overlapnewsystem.get_sites(), ketrotatematrix, correctionVector);
    
    overlapnewsystem.transform_operators(rotatematrix, ketrotatematrix);

    SpinBlock::store(forward, overlapnewsystem.get_sites(), overlapnewsystem, targetState, correctionVector);
    overlapsystem.clear(); overlapenvironment.clear(); overlapnewsystem.clear(); overlapnewenvironment.clear();

  }


  //<target|O|base>
  //save the updated overlap spinblock
  {
    int originalOutputlevel = dmrginp.outputlevel();
    //dmrginp.setOutputlevel() = -1;

    SpinBlock overlapBig;
    SpinBlock overlapsystem, overlapenvironment, overlapnewsystem, overlapnewenvironment;
    SpinBlock overlapsystemDot(systemDotStart, systemDotEnd, newSystem.get_integralIndex(), true);
    SpinBlock overlapenvironmentDot(environmentDotStart, environmentDotEnd, newSystem.get_integralIndex(), true);
    guessWaveTypes guesstype = sweepParams.get_block_iter() == 0 ? TRANSPOSE : TRANSFORM;
      
    DiagonalMatrix e;
    Sweep::makeSystemEnvironmentBigOverlapBlocks(systemsites, overlapsystemDot, overlapenvironmentDot,
						 overlapsystem, overlapnewsystem, overlapenvironment, overlapnewenvironment,
						 overlapBig, sweepParams, true, useSlater, newSystem.get_integralIndex(), 
						 targetState, baseState);

    Wavefunction iwave;
    GuessWave::guess_wavefunctions(iwave, e, overlapBig, guesstype, sweepParams.get_onedot(), baseState, true, 0.0);
    std::vector<Matrix> ketrotatematrix;
    DensityMatrix tracedMatrix;
    tracedMatrix.allocate(overlapnewsystem.get_ketStateInfo());
    operatorfunctions::MultiplyProduct(iwave, Transpose(const_cast<Wavefunction&> (iwave)), tracedMatrix, 1.0);
    int largeNumber = 1000000;
    if (!mpigetrank())
      double error = makeRotateMatrix(tracedMatrix, ketrotatematrix, largeNumber, sweepParams.get_keep_qstates());
    
#ifndef SERIAL
    mpi::communicator world;
    broadcast(world, ketrotatematrix, 0);
#endif
    
    iwave.SaveWavefunctionInfo (overlapBig.get_ketStateInfo(), overlapBig.get_leftBlock()->get_sites(), baseState);
    SaveRotationMatrix (overlapnewsystem.get_sites(), ketrotatematrix, baseState);
    
    overlapnewsystem.transform_operators(rotatematrix, ketrotatematrix);
    SpinBlock::store(forward, overlapnewsystem.get_sites(), overlapnewsystem, targetState, baseState);
    overlapsystem.clear(); overlapenvironment.clear(); overlapnewsystem.clear(); overlapnewenvironment.clear();

  }
  dmrginp.setOutputlevel() = originalOutputlevel;


  if (dmrginp.outputlevel() > 0)
    mcheck("after rotation and transformation of block");

  if (dmrginp.outputlevel() > 0){
    pout << *dmrginp.guessgenT<<" "<<*dmrginp.multiplierT<<" "<<*dmrginp.operrotT<< "  "<<globaltimer.totalwalltime()<<" timer "<<endl;
    pout << *dmrginp.makeopsT<<" makeops "<<endl;
    pout << *dmrginp.datatransfer<<" datatransfer "<<endl;
    pout <<"oneindexopmult   twoindexopmult   Hc  couplingcoeff"<<endl;  
    pout << *dmrginp.oneelecT<<" "<<*dmrginp.twoelecT<<" "<<*dmrginp.hmultiply<<" "<<*dmrginp.couplingcoeff<<" hmult"<<endl;
    pout << *dmrginp.buildsumblock<<" "<<*dmrginp.buildblockops<<" build block"<<endl;
    pout << "addnoise  S_0_opxop  S_1_opxop   S_2_opxop"<<endl;
    pout << *dmrginp.addnoise<<" "<<*dmrginp.s0time<<" "<<*dmrginp.s1time<<" "<<*dmrginp.s2time<<endl;
  }

}



double SpinAdapted::SweepResponse::do_one(SweepParams &sweepParams, const bool &warmUp, const bool &forward, 
					  const bool &restart, const int &restartSize, int targetState, 
					  int correctionVector, int baseState)
{
  int integralIndex = 0;
  SpinBlock system;

  std::vector<double> finalEnergy(1,1.0e10);
  double finalError = 0.;

  sweepParams.set_sweep_parameters();
  // a new renormalisation sweep routine
  pout << endl;
  if (forward)
    pout << "\t\t\t Starting sweep "<< sweepParams.set_sweep_iter()<<" in forwards direction"<<endl;
  else
    pout << "\t\t\t Starting sweep "<< sweepParams.set_sweep_iter()<<" in backwards direction" << endl;
  pout << "\t\t\t ============================================================================ " << endl;

  InitBlocks::InitStartingBlock (system,forward, targetState, targetState,
				 sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), 
				 restartSize, restart, warmUp, integralIndex);

  if(!restart)
    sweepParams.set_block_iter() = 0;

 
  if (dmrginp.outputlevel() > 0)
    pout << "\t\t\t Starting block is :: " << endl << system << endl;

 // if restart, just restoring an existing block --
  SpinBlock::store (forward, system.get_sites(), system, targetState, targetState);
  SpinBlock::store (forward, system.get_sites(), system, targetState, baseState);
  SpinBlock::store (forward, system.get_sites(), system, targetState, correctionVector);
  sweepParams.savestate(forward, system.get_sites().size());

  bool dot_with_sys = true;
  vector<int> syssites = system.get_sites();

 // get_n_iters() returns the number of blocking iterations needed in one sweep
  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); )
    {
      pout << "\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
      pout << "\t\t\t ----------------------------" << endl;
      if (dmrginp.outputlevel() > 0) {
	    if (forward)
	      pout << "\t\t\t Current direction is :: Forwards " << endl;
	    else
	      pout << "\t\t\t Current direction is :: Backwards " << endl;
      }


      if (sweepParams.get_block_iter() == 0 && sweepParams.get_sweep_iter() == 1)
	sweepParams.set_guesstype() = BASIC;
      else if (sweepParams.get_block_iter() != 0) 
	sweepParams.set_guesstype() = TRANSFORM;
      else
        sweepParams.set_guesstype() = TRANSPOSE;


      
      if (dmrginp.outputlevel() > 0)
         pout << "\t\t\t Blocking and Decimating " << endl;
	  
      SpinBlock newSystem; // new system after blocking and decimating

      //Need to substitute by:
      if (warmUp )
	StartUp(sweepParams, system, newSystem, dot_with_sys, 
		targetState, correctionVector, baseState);
      else {
	BlockAndDecimate(sweepParams, system, newSystem, warmUp, 
			 dot_with_sys, targetState, correctionVector, baseState);
      }
      
      //Need to substitute by?

      if (!warmUp ){

	//this criteria should work for state average or state specific because the lowest sweep energy is always the lowest of the average
	finalError = max(sweepParams.get_lowest_error(),finalError);
	finalEnergy[0] = min(sweepParams.get_lowest_energy()[0], finalEnergy[0]);
	pout << "final energy "<<finalEnergy[0]<<"  "<<sweepParams.get_lowest_energy()[0]<<endl;
      }
      
      system = newSystem;
      if (dmrginp.outputlevel() > 0){
	    pout << system<<endl;
	    system.printOperatorSummary();
      }
      
      //system size is going to be less than environment size
      if (forward && system.get_complementary_sites()[0] >= dmrginp.last_site()/2)
	    dot_with_sys = false;
      if (!forward && system.get_sites()[0]-1 < dmrginp.last_site()/2)
	    dot_with_sys = false;

      SpinBlock::store (forward, system.get_sites(), system, targetState, targetState);
      syssites = system.get_sites();
      if (dmrginp.outputlevel() > 0)
	pout << "\t\t\t saving state " << syssites.size() << endl;
      ++sweepParams.set_block_iter();
      
#ifndef SERIAL
      mpi::communicator world;
      world.barrier();
#endif
      sweepParams.savestate(forward, syssites.size());
      if (dmrginp.outputlevel() > 0)
         mcheck("at the end of sweep iteration");
    }

  //when we are doing twodot, we still need to do the last sweep to make sure that the
  //correctionVector and base wavefunction are propogated correctly across sweeps
  //especially when we switch from twodot to onedot algorithm
  if (!sweepParams.get_onedot() && !warmUp) {
      pout << "\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
      pout << "\t\t\t ----------------------------" << endl;
      if (dmrginp.outputlevel() > 0) {
	    if (forward)
	      pout << "\t\t\t Current direction is :: Forwards " << endl;
	    else
	      pout << "\t\t\t Current direction is :: Backwards " << endl;
      }
    sweepParams.set_onedot() = true;
    sweepParams.set_env_add() = 0;
    bool dot_with_sys = true;
    WavefunctionCanonicalize(sweepParams, system, warmUp, dot_with_sys, targetState, correctionVector, baseState);
    sweepParams.set_onedot() = false;
    sweepParams.set_env_add() = 1;
  }

  pout << "\t\t\t Largest Error for Sweep with " << sweepParams.get_keep_states() << " states is " << finalError << endl;
  pout << "\t\t\t Sweep Energy for Sweep with " << sweepParams.get_keep_states() << " states is " << finalEnergy[0] << endl;
  sweepParams.set_largest_dw() = finalError;
  

  pout << "\t\t\t ============================================================================ " << endl;

  // update the static number of iterations

  ++sweepParams.set_sweep_iter();

  return finalError;
}


void SpinAdapted::SweepResponse::StartUp (SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, const bool& dot_with_sys, int targetState, int correctionVector, int baseState)
{
  if (dmrginp.outputlevel() > 0) {
    mcheck("at the start of block and decimate");
    pout << "\t\t\t dot with system "<<dot_with_sys<<endl;
  }
  pout <<endl<< "\t\t\t Performing Blocking"<<endl;
  // figure out if we are going forward or backwards
  dmrginp.guessgenT -> start();

  bool forward = (system.get_sites() [0] == 0);
  SpinBlock systemDot, environmentDot;
  int systemDotStart, systemDotEnd, environmentDotStart, environmentDotEnd;
  int systemDotSize = sweepParams.get_sys_add() - 1;
  int environmentDotSize = sweepParams.get_env_add() - 1;
  if (forward)
  {
    systemDotStart = dmrginp.spinAdapted() ? 
      *system.get_sites().rbegin () + 1 : (*system.get_sites().rbegin ())/2 + 1 ;
    systemDotEnd = systemDotStart + systemDotSize;
    environmentDotStart = systemDotEnd + 1;
    environmentDotEnd = environmentDotStart + environmentDotSize;
  }
  else
  {
    systemDotStart = dmrginp.spinAdapted() ? 
      system.get_sites()[0] - 1 : (system.get_sites()[0])/2 - 1 ;
    systemDotEnd = systemDotStart - systemDotSize;
    environmentDotStart = systemDotEnd - 1;
    environmentDotEnd = environmentDotStart - environmentDotSize;
  }
  systemDot = SpinBlock(systemDotStart, systemDotEnd, system.get_integralIndex(), true);
  environmentDot = SpinBlock(environmentDotStart, environmentDotEnd, system.get_integralIndex(), true);
  SpinBlock environment, newEnvironment;

  SpinBlock big;  // new_sys = sys+sys_dot; new_env = env+env_dot; big = new_sys + new_env then renormalize to find new_sys(new)
  bool haveNormOps = dot_with_sys, haveCompOps = true;

  std::vector<Matrix> brarotateMatrix, ketrotateMatrix;
  //<target| H |target>
  {
    InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, correctionVector, correctionVector, 
				   sweepParams.get_sys_add(), dmrginp.direct(), system.get_integralIndex(), 
				   DISTRIBUTED_STORAGE, haveNormOps, haveCompOps);
    
    LoadRotationMatrix (newSystem.get_sites(), brarotateMatrix, correctionVector);
    Wavefunction targetWave;
    StateInfo s;
    targetWave.LoadWavefunctionInfo(s, newSystem.get_sites(), correctionVector);
    targetWave.SaveWavefunctionInfo(s, newSystem.get_sites(), targetState);
    SaveRotationMatrix (newSystem.get_sites(), brarotateMatrix, targetState);
    
#ifndef SERIAL
    mpi::communicator world;
    broadcast(world, brarotateMatrix, 0);
#endif
    
    dmrginp.operrotT -> start();
    newSystem.transform_operators(brarotateMatrix);
    dmrginp.operrotT -> stop();
  }

  std::vector<int> systemsites;
  if (dmrginp.spinAdapted())
    systemsites = system.get_sites();
  else {
    if (forward) {
      systemsites.push_back(0); systemsites.push_back(*system.get_sites().rbegin()/2);
    }
    else {
      systemsites.push_back(system.get_sites()[0]/2); systemsites.push_back(*system.get_sites().rbegin()/2);
    }
  }
  
  //<target| H |cV>
  {
    SpinBlock newOverlapSystemBlock, overlapSystemBlock;
    systemDot = SpinBlock(systemDotStart, systemDotEnd, system.get_integralIndex(), true);
    
    if ( (system.get_sites().size() == 1 && dmrginp.spinAdapted()) || (!dmrginp.spinAdapted()&&system.get_sites().size() == 2) ) {
      int restartSize = 0; bool restart=false, warmUp = false;
      InitBlocks::InitStartingBlock(overlapSystemBlock, forward, targetState, correctionVector, 
				    sweepParams.get_forward_starting_size(), 
				    sweepParams.get_backward_starting_size(), restartSize, 
				    restart, warmUp, system.get_integralIndex());
    }
    else {
      overlapSystemBlock.set_integralIndex() = system.get_integralIndex();
      SpinBlock::restore(forward, systemsites, overlapSystemBlock, targetState, correctionVector);
    }
    newOverlapSystemBlock.set_integralIndex() = system.get_integralIndex();
    newOverlapSystemBlock.initialise_op_array(OVERLAP, false);
    newOverlapSystemBlock.setstoragetype(DISTRIBUTED_STORAGE);
    newOverlapSystemBlock.BuildSumBlock (NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, overlapSystemBlock, 
					 systemDot);
    
    dmrginp.operrotT -> start();
    newOverlapSystemBlock.transform_operators(brarotateMatrix, brarotateMatrix);
    dmrginp.operrotT -> stop();
    SpinBlock::store(forward, newSystem.get_sites(), newOverlapSystemBlock, targetState, correctionVector);
  }


  //<target| H |base>
  {
    SpinBlock newOverlapSystemBlock, overlapSystemBlock;
    systemDot = SpinBlock(systemDotStart, systemDotEnd, system.get_integralIndex(), true);
    if ( (system.get_sites().size() == 1 && dmrginp.spinAdapted()) || (!dmrginp.spinAdapted()&&system.get_sites().size() == 2) ) {
      int restartSize = 0; bool restart=false, warmUp = false;
      InitBlocks::InitStartingBlock(overlapSystemBlock, forward, targetState, baseState, 
				  sweepParams.get_forward_starting_size(), 
				  sweepParams.get_backward_starting_size(), restartSize, 
				    restart, warmUp, system.get_integralIndex());
    }
    else {
      overlapSystemBlock.set_integralIndex() = system.get_integralIndex();
      SpinBlock::restore(forward, systemsites, overlapSystemBlock, targetState, baseState);
    }
    newOverlapSystemBlock.set_integralIndex() = system.get_integralIndex();
    newOverlapSystemBlock.initialise_op_array(OVERLAP, false);
    newOverlapSystemBlock.setstoragetype(DISTRIBUTED_STORAGE);
    newOverlapSystemBlock.BuildSumBlock (NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, overlapSystemBlock, 
					 systemDot);
    
    LoadRotationMatrix (newSystem.get_sites(), ketrotateMatrix, baseState);

#ifndef SERIAL
    mpi::communicator world;
    broadcast(world, ketrotateMatrix, 0);
#endif

    dmrginp.operrotT -> start();
    newOverlapSystemBlock.transform_operators(brarotateMatrix, ketrotateMatrix);
    dmrginp.operrotT -> stop();
    SpinBlock::store(forward, newSystem.get_sites(), newOverlapSystemBlock, targetState, baseState);
  }




  if (dmrginp.outputlevel() > 0)
    mcheck("after rotation and transformation of block");

  if (dmrginp.outputlevel() > 0){
    pout << *dmrginp.guessgenT<<" "<<*dmrginp.multiplierT<<" "<<*dmrginp.operrotT<< "  "<<globaltimer.totalwalltime()<<" timer "<<endl;
    pout << *dmrginp.makeopsT<<" makeops "<<endl;
    pout << *dmrginp.datatransfer<<" datatransfer "<<endl;
    pout <<"oneindexopmult   twoindexopmult   Hc  couplingcoeff"<<endl;  
    pout << *dmrginp.oneelecT<<" "<<*dmrginp.twoelecT<<" "<<*dmrginp.hmultiply<<" "<<*dmrginp.couplingcoeff<<" hmult"<<endl;
    pout << *dmrginp.buildsumblock<<" "<<*dmrginp.buildblockops<<" build block"<<endl;
    pout << "addnoise  S_0_opxop  S_1_opxop   S_2_opxop"<<endl;
    pout << *dmrginp.addnoise<<" "<<*dmrginp.s0time<<" "<<*dmrginp.s1time<<" "<<*dmrginp.s2time<<endl;
  }

}


void SpinAdapted::SweepResponse::WavefunctionCanonicalize (SweepParams &sweepParams, SpinBlock& system, const bool &useSlater, const bool& dot_with_sys, int targetState, int correctionVector, int baseState)
{
  if (dmrginp.outputlevel() > 0) {
    mcheck("at the start of block and decimate");
    pout << "\t\t\t dot with system "<<dot_with_sys<<endl;
  }
  pout <<endl<< "\t\t\t Performing Blocking"<<endl;
  // figure out if we are going forward or backwards
  dmrginp.guessgenT -> start();
  
  SpinBlock newSystem;
  
  bool forward = (system.get_sites() [0] == 0);
  SpinBlock systemDot, environmentDot;
  int systemDotStart, systemDotEnd, environmentDotStart, environmentDotEnd;
  int systemDotSize = sweepParams.get_sys_add() - 1;
  int environmentDotSize = 0;
  if (forward)
    {
      systemDotStart = dmrginp.spinAdapted() ? 
	*system.get_sites().rbegin () + 1 : (*system.get_sites().rbegin ())/2 + 1 ;
      systemDotEnd = systemDotStart + systemDotSize;
      environmentDotStart = systemDotEnd + 1;
      environmentDotEnd = environmentDotStart + environmentDotSize;
    }
  else
    {
      systemDotStart = dmrginp.spinAdapted() ? 
	system.get_sites()[0] - 1 : (system.get_sites()[0])/2 - 1 ;
      systemDotEnd = systemDotStart - systemDotSize;
      environmentDotStart = systemDotEnd - 1;
      environmentDotEnd = environmentDotStart - environmentDotSize;
    }
  systemDot = SpinBlock(systemDotStart, systemDotEnd, system.get_integralIndex(), true);
  environmentDot = SpinBlock(environmentDotStart, environmentDotEnd, system.get_integralIndex(), true);
  SpinBlock environment, newEnvironment;
  
  SpinBlock big;  // new_sys = sys+sys_dot; new_env = env+env_dot; big = new_sys + new_env then renormalize to find new_sys(new)

  InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, targetState, targetState, sweepParams.get_sys_add(), dmrginp.direct(), system.get_integralIndex(), 
				 DISTRIBUTED_STORAGE, false, true);

  newSystem.set_loopblock(false);  environmentDot.set_loopblock(false); 
  InitBlocks::InitBigBlock(newSystem, environmentDot, big);

  
  
  //analyse_operator_distribution(big);
  dmrginp.guessgenT -> stop();
  dmrginp.multiplierT -> start();
  std::vector<Matrix> rotatematrix;
  
  if (dmrginp.outputlevel() > 0)
    mcheck(""); 
  if (dmrginp.outputlevel() == 0) {
    if (!dot_with_sys && sweepParams.get_onedot()) pout << "\t\t\t System  Block"<<system;    
    else pout << "\t\t\t System  Block"<<newSystem;
    pout << "\t\t\t Environment Block"<<newEnvironment<<endl;
    pout << "\t\t\t Solving wavefunction "<<endl;
  }
  
  std::vector<Wavefunction> lowerStates;
  
  
  //make the baseState
  int originalOutputlevel = dmrginp.outputlevel();
  dmrginp.setOutputlevel() = -1;
  lowerStates.resize(2); //only make the 
  DiagonalMatrix e;
  guessWaveTypes guesstype = sweepParams.get_block_iter() == 0 ? TRANSPOSE : TRANSFORM;
  
  std::vector<int> systemsites;
  if (dmrginp.spinAdapted())
    systemsites = system.get_sites();
  else {
    if (forward) {
      systemsites.push_back(0); systemsites.push_back(*system.get_sites().rbegin()/2);
    }
    else {
      systemsites.push_back(system.get_sites()[0]/2); systemsites.push_back(*system.get_sites().rbegin()/2);
    }
  }
  
  
  
  Wavefunction targetWave;
  GuessWave::transform_previous_twodot_to_onedot_wavefunction(targetWave, big, targetState);

  
  DensityMatrix bratracedMatrix;
  bratracedMatrix.allocate(newSystem.get_braStateInfo());
  operatorfunctions::MultiplyProduct(targetWave, Transpose(const_cast<Wavefunction&> (targetWave)), bratracedMatrix, 1.0);
  int largeNumber = 1000000;
  if (!mpigetrank())
    double error = makeRotateMatrix(bratracedMatrix, rotatematrix, largeNumber, sweepParams.get_keep_qstates());
  
#ifndef SERIAL
  mpi::communicator world;
  broadcast(world, rotatematrix, 0);
#endif

  
  targetWave.SaveWavefunctionInfo (big.get_braStateInfo(), big.get_leftBlock()->get_sites(), targetState);
  SaveRotationMatrix (newSystem.get_sites(), rotatematrix, targetState);
  
  newSystem.transform_operators(rotatematrix);
  SpinBlock::store(forward, newSystem.get_sites(), newSystem, targetState, targetState);
  
  
  dmrginp.setOutputlevel() = originalOutputlevel;
  
  
  
  
  
  //<target|O|cV>
  //save the updated overlap spinblock
  //dmrginp.setOutputlevel() = -1;
  {
    
    SpinBlock overlapBig;
    SpinBlock overlapsystem, overlapenvironment, overlapnewsystem, overlapnewenvironment;
    SpinBlock overlapsystemDot(systemDotStart, systemDotEnd, system.get_integralIndex(), true);
    SpinBlock overlapenvironmentDot(environmentDotStart, environmentDotEnd, system.get_integralIndex(), true);
    guessWaveTypes guesstype = sweepParams.get_block_iter() == 0 ? TRANSPOSE : TRANSFORM;
    
    DiagonalMatrix e;

    overlapsystem.set_integralIndex() = system.get_integralIndex();
    SpinBlock::restore(forward, systemsites, overlapsystem, targetState, correctionVector);
    overlapnewsystem.initialise_op_array(OVERLAP, false);
    overlapnewsystem.setstoragetype(DISTRIBUTED_STORAGE);
    overlapnewsystem.BuildSumBlock (NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, overlapsystem, overlapsystemDot);
    overlapnewsystem.set_loopblock(false);  overlapenvironmentDot.set_loopblock(false); 
    InitBlocks::InitBigBlock(overlapnewsystem, overlapenvironmentDot, overlapBig);

    
    Wavefunction iwave;
    GuessWave::transform_previous_twodot_to_onedot_wavefunction(iwave, overlapBig, correctionVector);


    std::vector<Matrix> ketrotatematrix;
    DensityMatrix tracedMatrix;
    tracedMatrix.allocate(overlapnewsystem.get_ketStateInfo());
    operatorfunctions::MultiplyProduct(iwave, Transpose(const_cast<Wavefunction&> (iwave)), tracedMatrix, 1.0);
    int largeNumber = 1000000;
    if (!mpigetrank())
      double error = makeRotateMatrix(tracedMatrix, ketrotatematrix, largeNumber, sweepParams.get_keep_qstates());
    
#ifndef SERIAL
    mpi::communicator world;
    broadcast(world, ketrotatematrix, 0);
#endif
    
    iwave.SaveWavefunctionInfo (overlapBig.get_ketStateInfo(), overlapBig.get_leftBlock()->get_sites(), correctionVector);
    SaveRotationMatrix (overlapnewsystem.get_sites(), ketrotatematrix, correctionVector);
    
    overlapnewsystem.transform_operators(rotatematrix, ketrotatematrix);
    SpinBlock::store(forward, overlapnewsystem.get_sites(), overlapnewsystem, targetState, correctionVector);
    overlapsystem.clear(); overlapenvironment.clear(); overlapnewsystem.clear(); overlapnewenvironment.clear();
    
  }
  
  
  //<target|O|base>
  //save the updated overlap spinblock
  {
    
    SpinBlock overlapBig;
    SpinBlock overlapsystem, overlapenvironment, overlapnewsystem, overlapnewenvironment;
    SpinBlock overlapsystemDot(systemDotStart, systemDotEnd, system.get_integralIndex(), true);
    SpinBlock overlapenvironmentDot(environmentDotStart, environmentDotEnd, system.get_integralIndex(), true);
    guessWaveTypes guesstype = sweepParams.get_block_iter() == 0 ? TRANSPOSE : TRANSFORM;
    
    DiagonalMatrix e;

    overlapsystem.set_integralIndex() = system.get_integralIndex();
    SpinBlock::restore(forward, systemsites, overlapsystem, targetState, baseState);
    overlapnewsystem.initialise_op_array(OVERLAP, false);
    overlapnewsystem.setstoragetype(DISTRIBUTED_STORAGE);
    overlapnewsystem.BuildSumBlock (NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, overlapsystem, overlapsystemDot);
    overlapnewsystem.set_loopblock(false);  overlapenvironmentDot.set_loopblock(false); 
    InitBlocks::InitBigBlock(overlapnewsystem, overlapenvironmentDot, overlapBig);

    
    Wavefunction iwave;
    GuessWave::transform_previous_twodot_to_onedot_wavefunction(iwave, overlapBig, baseState);


    std::vector<Matrix> ketrotatematrix;
    DensityMatrix tracedMatrix;
    tracedMatrix.allocate(overlapnewsystem.get_ketStateInfo());
    operatorfunctions::MultiplyProduct(iwave, Transpose(const_cast<Wavefunction&> (iwave)), tracedMatrix, 1.0);
    int largeNumber = 1000000;
    if (!mpigetrank())
      double error = makeRotateMatrix(tracedMatrix, ketrotatematrix, largeNumber, sweepParams.get_keep_qstates());
    
#ifndef SERIAL
    mpi::communicator world;
    broadcast(world, ketrotatematrix, 0);
#endif
    
    iwave.SaveWavefunctionInfo (overlapBig.get_ketStateInfo(), overlapBig.get_leftBlock()->get_sites(), baseState);
    SaveRotationMatrix (overlapnewsystem.get_sites(), ketrotatematrix, baseState);
    
    overlapnewsystem.transform_operators(rotatematrix, ketrotatematrix);
    SpinBlock::store(forward, overlapnewsystem.get_sites(), overlapnewsystem, targetState, baseState);
    overlapsystem.clear(); overlapenvironment.clear(); overlapnewsystem.clear(); overlapnewenvironment.clear();
    
  }
  dmrginp.setOutputlevel() = originalOutputlevel;
  
  
  if (dmrginp.outputlevel() > 0)
    mcheck("after rotation and transformation of block");
  
  if (dmrginp.outputlevel() > 0){
    pout << *dmrginp.guessgenT<<" "<<*dmrginp.multiplierT<<" "<<*dmrginp.operrotT<< "  "<<globaltimer.totalwalltime()<<" timer "<<endl;
    pout << *dmrginp.makeopsT<<" makeops "<<endl;
    pout << *dmrginp.datatransfer<<" datatransfer "<<endl;
    pout <<"oneindexopmult   twoindexopmult   Hc  couplingcoeff"<<endl;  
    pout << *dmrginp.oneelecT<<" "<<*dmrginp.twoelecT<<" "<<*dmrginp.hmultiply<<" "<<*dmrginp.couplingcoeff<<" hmult"<<endl;
    pout << *dmrginp.buildsumblock<<" "<<*dmrginp.buildblockops<<" build block"<<endl;
    pout << "addnoise  S_0_opxop  S_1_opxop   S_2_opxop"<<endl;
    pout << *dmrginp.addnoise<<" "<<*dmrginp.s0time<<" "<<*dmrginp.s1time<<" "<<*dmrginp.s2time<<endl;
  }
  
}
