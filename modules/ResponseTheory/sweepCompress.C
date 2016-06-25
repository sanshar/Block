/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "guess_wavefunction.h"
#include "sweepCompress.h"
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

using namespace boost;
using namespace std;


void SpinAdapted::SweepCompress::BlockDecimateAndCompress (SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, const bool &useSlater, const bool& dot_with_sys, int targetState, int baseState)
{
  int sweepiter = sweepParams.get_sweep_iter();
  if (dmrginp.outputlevel() > 0)
    mcheck("at the start of block and decimate");
  p2out << "\t\t\t dot with system "<<dot_with_sys<<endl;
  p1out <<endl<< "\t\t\t Performing Blocking"<<endl;
  // figure out if we are going forward or backwards
  dmrginp.guessgenT -> start();
  bool forward = (system.get_sites() [0] == 0);
  SpinBlock systemDot;
  SpinBlock environment, environmentDot, newEnvironment;
  SpinBlock envDot, big;
  int systemDotStart, systemDotEnd;
  int environmentDotStart, environmentDotEnd, environmentStart, environmentEnd;
  int systemDotSize = sweepParams.get_sys_add() - 1;
  int environmentDotSize = sweepParams.get_env_add() -1;
  if (forward)
  {
    systemDotStart = dmrginp.spinAdapted() ? *system.get_sites().rbegin () + 1 : (*system.get_sites().rbegin ())/2 + 1 ;
    systemDotEnd = systemDotStart + systemDotSize;
    environmentDotStart = systemDotEnd + 1;
    environmentDotEnd = environmentDotStart + environmentDotSize;
  }
  else
  {
    systemDotStart = dmrginp.spinAdapted() ? system.get_sites()[0] - 1 : (system.get_sites()[0])/2 - 1 ;
    systemDotEnd = systemDotStart - systemDotSize;
    environmentDotStart = systemDotEnd - 1;
    environmentDotEnd = environmentDotStart - environmentDotSize;
  }
  systemDot = SpinBlock(systemDotStart, systemDotEnd, system.get_integralIndex(), true);
  environmentDot = SpinBlock(environmentDotStart, environmentDotEnd, system.get_integralIndex(), true);

  Sweep::makeSystemEnvironmentBigBlocks(system, systemDot, newSystem, environment, environmentDot, newEnvironment, big, sweepParams, dot_with_sys, useSlater, system.get_integralIndex(), targetState, baseState);


  //analyse_operator_distribution(big);
  dmrginp.guessgenT -> stop();
  dmrginp.multiplierT -> start();
  std::vector<Matrix> rotatematrix;

  if (dmrginp.outputlevel() > 0)
    mcheck(""); 
  if (!dot_with_sys && sweepParams.get_onedot())
    { pout << "\t\t\t System  Block"<<system; }
  else
    { pout << "\t\t\t System  Block"<<newSystem; }
  pout << "\t\t\t Environment Block"<<newEnvironment<<endl;
  p1out << "\t\t\t Solving wavefunction "<<endl;

  std::vector<Wavefunction> solution; solution.resize(1);

  DiagonalMatrix e;


  //read the 0th wavefunction which we keep on the ket side because by default the ket stateinfo is used to initialize wavefunction
  //also when you use spinblock operators to multiply a state, it does so from the ket side i.e.  H|ket>
  GuessWave::guess_wavefunctions(solution, e, big, sweepParams.set_guesstype(), sweepParams.get_onedot(), dot_with_sys, 0.0, baseState); 

#ifndef SERIAL
  mpi::communicator world;
  broadcast(world, solution, 0);
#endif
  
  multiply_h davidson_f(big, sweepParams.get_onedot());
  vector<Wavefunction> outputState; outputState.resize(1);
  outputState[0].AllowQuantaFor(big.get_leftBlock()->get_braStateInfo(), big.get_rightBlock()->get_braStateInfo(), dmrginp.effective_molecule_quantum_vec());
  outputState[0].set_onedot(sweepParams.get_onedot());
  outputState[0].Clear();

  davidson_f(solution[0], outputState[0]);
  double overlap = 1.0;

  SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
  sweepParams.set_lowest_energy() = std::vector<double>(1,overlap);

  SpinBlock newbig;

  if (sweepParams.get_onedot() && !dot_with_sys)
  {
    InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, baseState, targetState, systemDot.size(), dmrginp.direct(), system.get_integralIndex(), DISTRIBUTED_STORAGE, false, true);
    InitBlocks::InitBigBlock(newSystem, environment, newbig); 

    Wavefunction tempwave = solution[0];
    GuessWave::onedot_shufflesysdot(big.get_ketStateInfo(), newbig.get_ketStateInfo(), solution[0], tempwave);  
    solution[0] = tempwave;

    tempwave = outputState[0];
    GuessWave::onedot_shufflesysdot(big.get_braStateInfo(), newbig.get_braStateInfo(), outputState[0], tempwave);  
    outputState[0] = tempwave;

    big.get_rightBlock()->clear();
    big.clear();
  }
  else
    newbig = big;

  std::vector<Matrix> brarotateMatrix, ketrotateMatrix;
  DensityMatrix bratracedMatrix(newSystem.get_braStateInfo()), kettracedMatrix(newSystem.get_ketStateInfo());
  bratracedMatrix.allocate(newSystem.get_braStateInfo()); kettracedMatrix.allocate(newSystem.get_ketStateInfo());

  bratracedMatrix.makedensitymatrix(outputState, newbig, dmrginp.weights(sweepiter), 0.0, 0.0, true);
  if (sweepParams.get_noise() > NUMERICAL_ZERO) {
    pout << "adding noise  "<<trace(bratracedMatrix)<<"  "<<sweepiter<<"  "<<dmrginp.weights(sweepiter)[0]<<endl;
    bratracedMatrix.add_onedot_noise_forCompression(solution[0], newbig, sweepParams.get_noise()*max(1.0,trace(bratracedMatrix)));
    if (trace(bratracedMatrix) <1e-14) 
      bratracedMatrix.SymmetricRandomise();
      
    pout << "after noise  "<<trace(bratracedMatrix)<<"  "<<sweepParams.get_noise()<<endl;
  }
  environment.clear();
  newEnvironment.clear();

  kettracedMatrix.makedensitymatrix(solution, newbig, dmrginp.weights(sweepiter), 0.0, 0.0, true);
  double braerror, keterror;
  if (!mpigetrank()) {
    keterror = makeRotateMatrix(kettracedMatrix, ketrotateMatrix, newbig.get_rightBlock()->get_ketStateInfo().totalStates, 0);
    braerror = makeRotateMatrix(bratracedMatrix, brarotateMatrix, sweepParams.get_keep_states(), sweepParams.get_keep_qstates());
  }

#ifndef SERIAL
  broadcast(world, ketrotateMatrix, 0);
  broadcast(world, brarotateMatrix, 0);
#endif

  //assert(keterror < NUMERICAL_ZERO);
  pout << "\t\t\t Total ket discarded weight "<<keterror<<endl<<endl;
  pout << "\t\t\t Total bra discarded weight "<<braerror<<endl<<endl;

  sweepParams.set_lowest_error() = braerror;

  SaveRotationMatrix (newbig.get_leftBlock()->get_sites(), ketrotateMatrix, baseState);
  SaveRotationMatrix (newbig.get_leftBlock()->get_sites(), brarotateMatrix, targetState);
  solution[0].SaveWavefunctionInfo (newbig.get_ketStateInfo(), newbig.get_leftBlock()->get_sites(), baseState);
  outputState[0].SaveWavefunctionInfo (newbig.get_braStateInfo(), newbig.get_leftBlock()->get_sites(), targetState);

  p1out <<"\t\t\t Performing Renormalization "<<endl;
  newSystem.transform_operators(brarotateMatrix, ketrotateMatrix);




  if (dmrginp.outputlevel() > 0)
    mcheck("after rotation and transformation of block");

  p2out << *dmrginp.guessgenT<<" "<<*dmrginp.multiplierT<<" "<<*dmrginp.operrotT<< "  "<<globaltimer.totalwalltime()<<" timer "<<endl;
  p2out << *dmrginp.makeopsT<<" makeops "<<endl;
  p2out << *dmrginp.datatransfer<<" datatransfer "<<endl;
  p2out <<"oneindexopmult   twoindexopmult   Hc  couplingcoeff"<<endl;  
  p2out << *dmrginp.oneelecT<<" "<<*dmrginp.twoelecT<<" "<<*dmrginp.hmultiply<<" "<<*dmrginp.couplingcoeff<<" hmult"<<endl;
  p2out << *dmrginp.buildsumblock<<" "<<*dmrginp.buildblockops<<" build block"<<endl;
  p2out << "addnoise  S_0_opxop  S_1_opxop   S_2_opxop"<<endl;
  p2out << *dmrginp.addnoise<<" "<<*dmrginp.s0time<<" "<<*dmrginp.s1time<<" "<<*dmrginp.s2time<<endl;

}

double SpinAdapted::SweepCompress::do_one(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize, int targetState, int baseState)
{
  int integralIndex = 0;
  SpinBlock system;
  const int nroots = dmrginp.nroots(sweepParams.get_sweep_iter());

  std::vector<double> finalEnergy(nroots,-1.0e10);
  std::vector<double> finalEnergy_spins(nroots,0.);
  double finalError = 0.;
  if (restart) {
    finalEnergy = sweepParams.get_lowest_energy();
    finalEnergy_spins = sweepParams.get_lowest_energy();
    finalError = sweepParams.get_lowest_error();
  }

  sweepParams.set_sweep_parameters();
  // a new renormalisation sweep routine
  pout << endl;
  if (forward)
    { pout << "\t\t\t Starting sweep "<< sweepParams.set_sweep_iter()<<" in forwards direction"<<endl; }
  else
    { pout << "\t\t\t Starting sweep "<< sweepParams.set_sweep_iter()<<" in backwards direction" << endl; }
  pout << "\t\t\t ============================================================================ " << endl;

  InitBlocks::InitStartingBlock (system,forward, baseState, targetState, sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), restartSize, restart, warmUp, integralIndex);
  if(!restart)
    sweepParams.set_block_iter() = 0;

 
  p2out << "\t\t\t Starting block is :: " << endl << system << endl;

  SpinBlock::store (forward, system.get_sites(), system, targetState, baseState); // if restart, just restoring an existing block --
  sweepParams.savestate(forward, system.get_sites().size());
  bool dot_with_sys = true;
  vector<int> syssites = system.get_sites();

  if (restart)
  {
    if (forward && system.get_complementary_sites()[0] >= dmrginp.last_site()/2)
      dot_with_sys = false;
    if (!forward && system.get_sites()[0]-1 < dmrginp.last_site()/2)
      dot_with_sys = false;
  }
  if (dmrginp.outputlevel() > 0)
    mcheck("at the very start of sweep");  // just timer

  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); ) // get_n_iters() returns the number of blocking iterations needed in one sweep
    {
      pout << "\n\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
      pout << "\t\t\t ----------------------------" << endl;
      if (forward)
	{ p1out << "\t\t\t Current direction is :: Forwards " << endl; }
      else
	{ p1out << "\t\t\t Current direction is :: Backwards " << endl; }

      if (sweepParams.get_block_iter() != 0) 
	sweepParams.set_guesstype() = TRANSFORM;
      else
        sweepParams.set_guesstype() = TRANSPOSE;


      
      p1out << "\t\t\t Blocking and Decimating " << endl;
	  
      SpinBlock newSystem; // new system after blocking and decimating

      //Need to substitute by:
      if (warmUp )
	Startup(sweepParams, system, newSystem, dot_with_sys, targetState, baseState);
      else {
	BlockDecimateAndCompress (sweepParams, system, newSystem, false, dot_with_sys, targetState, baseState);
      }
      
      //Need to substitute by?

      if (!warmUp ){

	//this criteria should work for state average or state specific because the lowest sweep energy is always the lowest of the average
	finalError = max(sweepParams.get_lowest_error(),finalError);
	finalEnergy[0] = max(sweepParams.get_lowest_energy()[0], finalEnergy[0]);
	pout << "final energy "<<finalEnergy[0]<<"  "<<sweepParams.get_lowest_energy()[0]<<endl;
      }
      
      system = newSystem;
      p2out << system<<endl;
      p2out << system.get_braStateInfo()<<endl;
      system.printOperatorSummary();
      
      //system size is going to be less than environment size
      if (forward && system.get_complementary_sites()[0] >= dmrginp.last_site()/2)
	    dot_with_sys = false;
      if (!forward && system.get_sites()[0]-1 < dmrginp.last_site()/2)
	    dot_with_sys = false;

      SpinBlock::store (forward, system.get_sites(), system, targetState, baseState);	 	
      syssites = system.get_sites();
      p1out << "\t\t\t saving state " << syssites.size() << endl;
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
      pout << "\n\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
      pout << "\t\t\t ----------------------------" << endl;
      if (forward)
	{ p1out << "\t\t\t Current direction is :: Forwards " << endl; }
      else
	{ p1out << "\t\t\t Current direction is :: Backwards " << endl; }
    sweepParams.set_onedot() = true;
    sweepParams.set_env_add() = 0;
    bool dot_with_sys = true;
    WavefunctionCanonicalize(sweepParams, system, warmUp, dot_with_sys, targetState, baseState);
    sweepParams.set_onedot() = false;
    sweepParams.set_env_add() = 1;
  }


  pout << "\t\t\t Largest Error for Sweep with " << sweepParams.get_keep_states() << " states is " << finalError << endl;
  pout << "\t\t\t Largest overlap for Sweep with " << sweepParams.get_keep_states() << " states is " << finalEnergy[0] << endl;
  sweepParams.set_largest_dw() = finalError;
  

  pout << "\t\t\t ============================================================================ " << endl;

  // update the static number of iterations

  ++sweepParams.set_sweep_iter();

  return finalError;
}


void SpinAdapted::SweepCompress::Startup (SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, const bool& dot_with_sys, int targetState, int baseState)
{
  bool useSlater = false;
  p1out <<endl<< "\t\t\t Performing Blocking"<<endl;
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
  systemDot = SpinBlock(systemDotStart, systemDotEnd, system.get_integralIndex(), false);
  SpinBlock environment, environmentDot, newEnvironment;


  int environmentDotStart, environmentDotEnd, environmentStart, environmentEnd;
  int environmentDotSize = sweepParams.get_env_add() -1;
  if (environmentDotSize <0) environmentDotSize = 0 ; 
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

  if (!sweepParams.get_onedot())
    environmentDot = SpinBlock(environmentDotStart, environmentDotEnd, system.get_integralIndex(), false);

  const int nexact = forward ? sweepParams.get_forward_starting_size() : sweepParams.get_backward_starting_size();

  //before halfway put the sysdot with system otherwise with environment
  if (!sweepParams.get_onedot()) {
      dmrginp.datatransfer -> start();
      system.addAdditionalCompOps();
      dmrginp.datatransfer -> stop();

      bool haveNormOps = dot_with_sys, haveCompOps = true;
      InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, baseState, targetState, sweepParams.get_sys_add(), dmrginp.direct(), 
				     system.get_integralIndex(), DISTRIBUTED_STORAGE, haveNormOps, haveCompOps);
      if (dmrginp.outputlevel() > 0)
         mcheck("");

      InitBlocks::InitNewEnvironmentBlock(environment, environmentDot, newEnvironment, system, systemDot, baseState, baseState,
					  sweepParams.get_sys_add(), sweepParams.get_env_add(), forward, dmrginp.direct(),
					  sweepParams.get_onedot(), nexact, useSlater, system.get_integralIndex(), !haveNormOps, haveCompOps, dot_with_sys);
      if (dmrginp.outputlevel() > 0)
         mcheck("");
  }
  else {
    dmrginp.datatransfer -> start();
    system.addAdditionalCompOps();
    dmrginp.datatransfer -> stop();

    bool haveNormOps = dot_with_sys, haveCompOps = true;
    if (dot_with_sys) {
      InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, targetState, baseState, sweepParams.get_sys_add(), dmrginp.direct(), system.get_integralIndex(), DISTRIBUTED_STORAGE, haveNormOps, haveCompOps);

    }
    InitBlocks::InitNewEnvironmentBlock(environment, systemDot, newEnvironment, system, systemDot, baseState, baseState,
					sweepParams.get_sys_add(), sweepParams.get_env_add(), forward, dmrginp.direct(),
					sweepParams.get_onedot(), nexact, useSlater, system.get_integralIndex(), !haveNormOps, haveCompOps, dot_with_sys);
  }
  SpinBlock big;  // new_sys = sys+sys_dot; new_env = env+env_dot; big = new_sys + new_env then renormalize to find new_sys(new)
  if (dot_with_sys) {
    newSystem.set_loopblock(true);
    system.set_loopblock(false);
    newEnvironment.set_loopblock(false);
    if (!sweepParams.get_onedot())
      environment.set_loopblock(false);
    InitBlocks::InitBigBlock(newSystem, newEnvironment, big); 
  }
  else{
    if (sweepParams.get_onedot()) {
      system.set_loopblock(false);
      newEnvironment.set_loopblock(true);
      environment.set_loopblock(true);
      InitBlocks::InitBigBlock(system, newEnvironment, big); 
    }
    else {
      newSystem.set_loopblock(false);
      system.set_loopblock(false);
      newEnvironment.set_loopblock(true);
      environment.set_loopblock(false);
      InitBlocks::InitBigBlock(newSystem, newEnvironment, big); 
    }
  }
  //analyse_operator_distribution(big);
  dmrginp.guessgenT -> stop();
  dmrginp.multiplierT -> start();
  std::vector<Matrix> rotatematrix;

  if (dmrginp.outputlevel() > 0)
    mcheck(""); 
  if (!dot_with_sys && sweepParams.get_onedot()) {
    pout << "\t\t\t System  Block"<<system;
    pout << "\t\t\t Environment Block"<<newEnvironment<<endl;
  }
  else {
    pout << "\t\t\t System  Block"<<newSystem;
    pout << "\t\t\t Environment Block"<<newEnvironment<<endl;
  }
  p1out << "\t\t\t Solving wavefunction "<<endl;

  std::vector<Wavefunction> solution; solution.resize(1);

  DiagonalMatrix e;
  e.ReSize(big.get_stateInfo().totalStates); e= 0;

  //read the 0th wavefunction which we keep on the ket side because by default the ket stateinfo is used to initialize wavefunction
  //also when you use spinblock operators to multiply a state, it does so from the ket side i.e.  H|ket>
  GuessWave::guess_wavefunctions(solution, e, big, sweepParams.set_guesstype(), sweepParams.get_onedot(), dot_with_sys, 0.0, baseState); 

  SpinBlock newbig;

  if (sweepParams.get_onedot() && !dot_with_sys)
  {
    InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, targetState, baseState, systemDot.size(), dmrginp.direct(), system.get_integralIndex(), DISTRIBUTED_STORAGE, false, true);
    InitBlocks::InitBigBlock(newSystem, environment, newbig); 

    Wavefunction tempwave = solution[0];
    tempwave.Clear();
    GuessWave::onedot_shufflesysdot(big.get_ketStateInfo(), newbig.get_ketStateInfo(), solution[0], tempwave);  
    solution[0] = tempwave;

    big.get_rightBlock()->clear();
    big.clear();
  }
  else
    newbig = big;


  SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));

  std::vector<Matrix> ketrotateMatrix, brarotateMatrix;
  DensityMatrix kettracedMatrix(newSystem.get_ketStateInfo());
  DensityMatrix bratracedMatrix(newSystem.get_braStateInfo());

  kettracedMatrix.allocate(newSystem.get_ketStateInfo());
  bratracedMatrix.allocate(newSystem.get_braStateInfo());


  bratracedMatrix.makedensitymatrix(solution, newbig, dmrginp.weights(0), 0.0, 
				    0.0, true);
  //bratracedMatrix.add_onedot_noise_forCompression(solution[0], newbig, sweepParams.get_noise()*max(1.0, trace(bratracedMatrix)));
  //bratracedMatrix.add_twodot_noise(newbig, sweepParams.get_noise());
  environment.clear();
  newEnvironment.clear();


  kettracedMatrix.makedensitymatrix(solution, newbig, dmrginp.weights(0), 0.0, 
				    0.0, true);
  //kettracedMatrix.add_twodot_noise(newbig, sweepParams.get_noise());
  double keterror, braerror;
  if (!mpigetrank()) {
    keterror = makeRotateMatrix(kettracedMatrix, ketrotateMatrix, newbig.get_rightBlock()->get_ketStateInfo().totalStates, 0);
    braerror = makeRotateMatrix(bratracedMatrix, brarotateMatrix, newbig.get_rightBlock()->get_braStateInfo().totalStates, 0);
  }
#ifndef SERIAL
  mpi::communicator world;
  broadcast(world, ketrotateMatrix, 0);
  broadcast(world, brarotateMatrix, 0);
#endif

  //assert(keterror < NUMERICAL_ZERO);
  p1out <<"\t\t\t Performing Renormalization "<<endl;
  pout << "\t\t\t Total ket discarded weight "<<keterror<<endl<<endl;
  pout << "\t\t\t Total bra discarded weight "<<braerror<<endl<<endl;
  sweepParams.set_lowest_error() = keterror;

  SaveRotationMatrix (newbig.get_leftBlock()->get_sites(), ketrotateMatrix, baseState);
  SaveRotationMatrix (newbig.get_leftBlock()->get_sites(), brarotateMatrix, targetState);
  solution[0].SaveWavefunctionInfo (newbig.get_ketStateInfo(), newbig.get_leftBlock()->get_sites(), baseState);
  solution[0].SaveWavefunctionInfo (newbig.get_braStateInfo(), newbig.get_leftBlock()->get_sites(), targetState);


  newSystem.transform_operators(brarotateMatrix, ketrotateMatrix);



  if (dmrginp.outputlevel() > 0)
    mcheck("after rotation and transformation of block");

  p2out << *dmrginp.guessgenT<<" "<<*dmrginp.multiplierT<<" "<<*dmrginp.operrotT<< "  "<<globaltimer.totalwalltime()<<" timer "<<endl;
  p2out << *dmrginp.makeopsT<<" makeops "<<endl;
  p2out << *dmrginp.datatransfer<<" datatransfer "<<endl;
  p2out <<"oneindexopmult   twoindexopmult   Hc  couplingcoeff"<<endl;  
  p2out << *dmrginp.oneelecT<<" "<<*dmrginp.twoelecT<<" "<<*dmrginp.hmultiply<<" "<<*dmrginp.couplingcoeff<<" hmult"<<endl;
  p2out << *dmrginp.buildsumblock<<" "<<*dmrginp.buildblockops<<" build block"<<endl;
  p2out << "addnoise  S_0_opxop  S_1_opxop   S_2_opxop"<<endl;
  p3out << *dmrginp.addnoise<<" "<<*dmrginp.s0time<<" "<<*dmrginp.s1time<<" "<<*dmrginp.s2time<<endl;

}

void SpinAdapted::SweepCompress::WavefunctionCanonicalize (SweepParams &sweepParams, SpinBlock& system, const bool &useSlater, const bool& dot_with_sys, int correctionVector, int baseState)
{
  if (dmrginp.outputlevel() > 0)
    mcheck("at the start of block and decimate");
  p2out << "\t\t\t dot with system "<<dot_with_sys<<endl;
  p1out <<endl<< "\t\t\t Performing Blocking"<<endl;
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
  vector<int> sitesenvdot(environmentDotSize+1, 0);
  int index = 0;
  for (int i=min(environmentDotStart, environmentDotEnd); i<max(environmentDotStart, environmentDotEnd)+1; i++) {
    sitesenvdot[index] = (i);
    index++;
  }

  SpinBlock::restore(!forward, sitesenvdot, environmentDot, correctionVector, baseState); 

  SpinBlock environment, newEnvironment;
  
  SpinBlock big;  // new_sys = sys+sys_dot; new_env = env+env_dot; big = new_sys + new_env then renormalize to find new_sys(new)

  system.addAdditionalCompOps();
  InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, correctionVector, baseState, sweepParams.get_sys_add(), dmrginp.direct(), 
				 system.get_integralIndex(), DISTRIBUTED_STORAGE, false, true);

  newSystem.set_loopblock(false);  environmentDot.set_loopblock(false); 
  InitBlocks::InitBigBlock(newSystem, environmentDot, big);

  
  
  //analyse_operator_distribution(big);
  dmrginp.guessgenT -> stop();
  dmrginp.multiplierT -> start();
  std::vector<Matrix> rotatematrix;
  
  if (dmrginp.outputlevel() > 0)
    mcheck(""); 
  if (!dot_with_sys && sweepParams.get_onedot())
    { pout << "\t\t\t System  Block"<<system; }
  else
    { pout << "\t\t\t System  Block"<<newSystem; }
  pout << "\t\t\t Environment Block"<<environmentDot<<endl;
  p1out << "\t\t\t Solving wavefunction "<<endl;
  
  
  
  //make the baseState
  int originalOutputlevel = dmrginp.outputlevel();
  //dmrginp.setOutputlevel() = -1;
  
  
  
  std::vector< Wavefunction> solution(1);
  if (!mpigetrank())
    GuessWave::transform_previous_twodot_to_onedot_wavefunction(solution[0], big, baseState);
#ifndef SERIAL
    mpi::communicator world;
    broadcast(world, solution, 0);
#endif

  multiply_h davidson_f(big, sweepParams.get_onedot());
  vector<Wavefunction> outputState; outputState.resize(1);
  outputState[0].AllowQuantaFor(big.get_leftBlock()->get_braStateInfo(), big.get_rightBlock()->get_braStateInfo(), dmrginp.effective_molecule_quantum_vec());
  outputState[0].set_onedot(sweepParams.get_onedot());
  outputState[0].Clear();

  davidson_f(solution[0], outputState[0]);

  
  std::vector<Matrix> ketrotateMatrix, brarotateMatrix;
  DensityMatrix bratracedMatrix, kettracedMatrix;
  bratracedMatrix.allocate(newSystem.get_braStateInfo());
  kettracedMatrix.allocate(newSystem.get_ketStateInfo());

  bratracedMatrix.makedensitymatrix(outputState, big, dmrginp.weights(0), 0.0, 0.0, true);

  kettracedMatrix.makedensitymatrix(solution, big, dmrginp.weights(0), 0.0, 0.0, true);
  double braerror, keterror;
  int largeNumber = 1000000;
  if (!mpigetrank()) {
    keterror = makeRotateMatrix(kettracedMatrix, ketrotateMatrix, largeNumber, 0);
    braerror = makeRotateMatrix(bratracedMatrix, brarotateMatrix, sweepParams.get_keep_states(), sweepParams.get_keep_qstates());
  }

#ifndef SERIAL
  broadcast(world, ketrotateMatrix, 0);
  broadcast(world, brarotateMatrix, 0);
#endif


  pout << "\t\t\t Total ket discarded weight "<<keterror<<endl<<endl;
  pout << "\t\t\t Total bra discarded weight "<<braerror<<endl<<endl;

  sweepParams.set_lowest_error() = braerror;

  SaveRotationMatrix (big.get_leftBlock()->get_sites(), ketrotateMatrix, baseState);
  SaveRotationMatrix (big.get_leftBlock()->get_sites(), brarotateMatrix, correctionVector);
  solution[0].SaveWavefunctionInfo (big.get_ketStateInfo(), big.get_leftBlock()->get_sites(), baseState);
  outputState[0].SaveWavefunctionInfo (big.get_braStateInfo(), big.get_leftBlock()->get_sites(), correctionVector);

  p1out <<"\t\t\t Performing Renormalization "<<endl;
  newSystem.transform_operators(brarotateMatrix, ketrotateMatrix);


  
  p2out << *dmrginp.guessgenT<<" "<<*dmrginp.multiplierT<<" "<<*dmrginp.operrotT<< "  "<<globaltimer.totalwalltime()<<" timer "<<endl;
  p2out << *dmrginp.makeopsT<<" makeops "<<endl;
  p2out << *dmrginp.datatransfer<<" datatransfer "<<endl;
  p2out <<"oneindexopmult   twoindexopmult   Hc  couplingcoeff"<<endl;  
  p2out << *dmrginp.oneelecT<<" "<<*dmrginp.twoelecT<<" "<<*dmrginp.hmultiply<<" "<<*dmrginp.couplingcoeff<<" hmult"<<endl;
  p2out << *dmrginp.buildsumblock<<" "<<*dmrginp.buildblockops<<" build block"<<endl;
  p2out << "addnoise  S_0_opxop  S_1_opxop   S_2_opxop"<<endl;
  p3out << *dmrginp.addnoise<<" "<<*dmrginp.s0time<<" "<<*dmrginp.s1time<<" "<<*dmrginp.s2time<<endl;
  
}
