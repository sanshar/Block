/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "sweep.h"
#include "global.h"
#include "solver.h"
#include "initblocks.h"
#include "MatrixBLAS.h"
#ifndef SERIAL
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif
#include "rotationmat.h"
#include "density.h"

#ifdef MOLPRO
#include "global/CxOutputStream.h"
#define pout if (dmrginp.outputlevel() < 0) xout
#endif


using namespace boost;
using namespace std;


void SpinAdapted::Sweep::BlockAndDecimate (SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, const bool &useSlater, const bool& dot_with_sys)
{
  if (dmrginp.outputlevel() > 0) {
    mcheck("at the start of block and decimate");
    pout << "\t\t\t dot with system "<<dot_with_sys<<endl;
  }
  pout <<endl<< "\t\t\t Performing Blocking"<<endl;
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
    environmentDot = SpinBlock(environmentDotStart, environmentDotEnd);
  
  const int nexact = forward ? sweepParams.get_forward_starting_size() : sweepParams.get_backward_starting_size();

  //before halfway put the sysdot with system otherwise with environment
  if (!sweepParams.get_onedot()) {
      dmrginp.datatransfer -> start();
      system.addAdditionalCompOps();
      dmrginp.datatransfer -> stop();
      InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, sweepParams.get_sys_add(), dmrginp.direct(), 
      			     DISTRIBUTED_STORAGE, dot_with_sys, true);
      if (dmrginp.outputlevel() > 0)
         mcheck("");

      InitBlocks::InitNewEnvironmentBlock(environment, environmentDot, newEnvironment, system, systemDot,
					  sweepParams.get_sys_add(), sweepParams.get_env_add(), forward, dmrginp.direct(),
					  sweepParams.get_onedot(), nexact, useSlater, !dot_with_sys, true, dot_with_sys);
      if (dmrginp.outputlevel() > 0)
         mcheck("");
  }
  else {
    dmrginp.datatransfer -> start();
    system.addAdditionalCompOps();
    dmrginp.datatransfer -> stop();
    if (dot_with_sys) {
      InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, sweepParams.get_sys_add(), dmrginp.direct(), DISTRIBUTED_STORAGE, dot_with_sys, true);

    }
    InitBlocks::InitNewEnvironmentBlock(environment, systemDot, newEnvironment, system, systemDot,
					sweepParams.get_sys_add(), sweepParams.get_env_add(), forward, dmrginp.direct(),
					sweepParams.get_onedot(), nexact, useSlater, !dot_with_sys, true, dot_with_sys);
  }
  SpinBlock big;
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
  if (dmrginp.outputlevel() == 0) {
    if (!dot_with_sys && sweepParams.get_onedot()) {
      pout << "\t\t\t System  Block"<<system;
      pout << "\t\t\t Environment Block"<<newEnvironment<<endl;
    }
    else {
      pout << "\t\t\t System  Block"<<newSystem;
      pout << "\t\t\t Environment Block"<<newEnvironment<<endl;
    }
    pout << "\t\t\t Solving wavefunction "<<endl;
  }

  newSystem.RenormaliseFrom (sweepParams.set_lowest_energy(), sweepParams.set_lowest_energy_spins(), sweepParams.set_lowest_error(), 
                             rotatematrix, sweepParams.get_keep_states(), 
                             sweepParams.get_keep_qstates(), sweepParams.get_davidson_tol(), big, sweepParams.get_guesstype(), sweepParams.get_noise(), 
                             sweepParams.get_additional_noise(), sweepParams.get_onedot(), system, systemDot, environmentDot, environment, 
			     dot_with_sys, useSlater, sweepParams.get_sweep_iter());

  std::vector<StateInfo> storeStates(3);
  storeStates[0] = newSystem.get_stateInfo();
  if(dot_with_sys)
    storeStates[1] = newEnvironment.get_stateInfo();
  else
    storeStates[1] = environment.get_stateInfo();

  if (dmrginp.outputlevel() > 0)
    mcheck("");
  environment.clear();
  newEnvironment.clear();

  pout <<"\t\t\t Performing Renormalization "<<endl;
  pout << "\t\t\t Total discarded weight "<<sweepParams.set_lowest_error()<<endl<<endl;

  dmrginp.multiplierT -> stop();
  dmrginp.operrotT -> start();

  newSystem.transform_operators(rotatematrix);
  storeStates[2] = newSystem.get_stateInfo();
  dmrginp.operrotT -> stop();
  if (dmrginp.outputlevel() > 0)
    mcheck("after rotation and transformation of block");

  if (dmrginp.outputlevel() > 0){
    pout << dmrginp.guessgenT<<" "<<dmrginp.multiplierT<<" "<<dmrginp.operrotT<< "  "<<globaltimer.totalwalltime()<<" timer "<<endl;
    pout << dmrginp.makeopsT<<" makeops "<<endl;
    pout << dmrginp.datatransfer<<" datatransfer "<<endl;
    pout <<"oneindexopmult   twoindexopmult   Hc  couplingcoeff"<<endl;  
    pout << dmrginp.oneelecT<<" "<<dmrginp.twoelecT<<" "<<dmrginp.hmultiply<<" "<<dmrginp.couplingcoeff<<" hmult"<<endl;
    pout << dmrginp.buildsumblock<<" "<<dmrginp.buildblockops<<" build block"<<endl;
    pout << "addnoise  S_0_opxop  S_1_opxop   S_2_opxop"<<endl;
    pout << dmrginp.addnoise<<" "<<dmrginp.s0time<<" "<<dmrginp.s1time<<" "<<dmrginp.s2time<<endl;
  }

}

double SpinAdapted::Sweep::do_one(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize)
{

  SpinBlock system;
  const int nroots = dmrginp.nroots(sweepParams.get_sweep_iter());
  std::vector<double> finalEnergy(nroots,1.0e10);
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
    pout << "\t\t\t Starting sweep "<< sweepParams.set_sweep_iter()<<" in forwards direction"<<endl;
  else
    pout << "\t\t\t Starting sweep "<< sweepParams.set_sweep_iter()<<" in backwards direction" << endl;
  pout << "\t\t\t ============================================================================ " << endl;

  InitBlocks::InitStartingBlock (system,forward, sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), restartSize, restart, warmUp);
  if(!restart)
    sweepParams.set_block_iter() = 0;

 
  if (dmrginp.outputlevel() > 0)
    pout << "\t\t\t Starting block is :: " << endl << system << endl;

  SpinBlock::store (forward, system.get_sites(), system); // if restart, just restoring an existing block --
  sweepParams.savestate(forward, system.get_sites().size());
  bool dot_with_sys = true;
  vector<int> syssites;
  {
    syssites = system.get_sites();
  }

  if (restart)
  {
    if (forward && system.get_complementary_sites()[0] >= dmrginp.last_site()/2)
      dot_with_sys = false;
    if (!forward && system.get_sites()[0]-1 < dmrginp.last_site()/2)
      dot_with_sys = false;
  }
  if (dmrginp.outputlevel() > 0)
    mcheck("at the very start of sweep");

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

      if (dmrginp.no_transform() || (sweepParams.get_sweep_iter()-sweepParams.get_restart_iter() == 0 && sweepParams.get_block_iter() == 0))
	      sweepParams.set_guesstype() = BASIC;
      else if (!warmUp && sweepParams.get_block_iter() != 0) 
  	    sweepParams.set_guesstype() = TRANSFORM;
      else if (!warmUp && sweepParams.get_block_iter() == 0 && 
                ((dmrginp.algorithm_method() == TWODOT_TO_ONEDOT && dmrginp.twodot_to_onedot_iter() != sweepParams.get_sweep_iter()) ||
                  dmrginp.algorithm_method() != TWODOT_TO_ONEDOT))
        sweepParams.set_guesstype() = TRANSPOSE;
      else
        sweepParams.set_guesstype() = BASIC;

      
      if (dmrginp.outputlevel() > 0)
         pout << "\t\t\t Blocking and Decimating " << endl;
	  
      SpinBlock newSystem;

      //Need to substitute by:
      if (warmUp && (sym=="dinfh"||sym=="trans"))
      //if (warmUp)// && (sym=="dinfh"||sym=="trans"))
         Startup(sweepParams, system, newSystem);
      else {
         if (sweepParams.set_sweep_iter() == 1 && sweepParams.get_block_iter() == 0)
           sweepParams.set_guesstype() = BASIC;
         BlockAndDecimate (sweepParams, system, newSystem, warmUp, dot_with_sys);
      }
      
      //Need to substitute by?
      if (!warmUp || !(sym == "dinfh"||sym=="trans") ){
      //if (!warmUp){// || !(sym == "dinfh"||sym=="trans") ){
	
	for(int j=0;j<nroots;++j)
   {
#ifndef MOLPRO
	  pout << "\t\t\t Total block energy for State [ " << j << 
	    " ] with " << sweepParams.get_keep_states()<<" States :: " << sweepParams.get_lowest_energy()[j]+dmrginp.get_coreenergy() <<endl;              
#else 
     //We might want to relax the output restrictions here, so it prints out with outputlevel=0
     if (dmrginp.outputlevel() < 0) {
        pout << "\t\t\t Total block energy for State [ " << j << 
          " ] with " << sweepParams.get_keep_states()<<" States :: " << fixed << setprecision(10) << sweepParams.get_lowest_energy()[j]+dmrginp.get_coreenergy() <<endl;              
      }
#endif
   }
	
	finalEnergy_spins = ((sweepParams.get_lowest_energy()[0] < finalEnergy[0]) ? sweepParams.get_lowest_energy_spins() : finalEnergy_spins);
	finalEnergy = ((sweepParams.get_lowest_energy()[0] < finalEnergy[0]) ? sweepParams.get_lowest_energy() : finalEnergy);
	finalError = max(sweepParams.get_lowest_error(),finalError);
	pout << endl;
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

      SpinBlock::store (forward, system.get_sites(), system);	 	
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
  for(int j=0;j<nroots;++j)
    pout << "\t\t\t Finished Sweep with " << sweepParams.get_keep_states() << " states and sweep energy for State [ " << j 
	 << " ] with Spin [ " << dmrginp.molecule_quantum().get_s()  << " ] :: " << finalEnergy[j]+dmrginp.get_coreenergy() << endl;
  pout << "\t\t\t Largest Error for Sweep with " << sweepParams.get_keep_states() << " states is " << finalError << endl;
  for(int j=0;j<nroots;++j){
    if (mpigetrank() == 0) {
#ifndef MOLPRO
      printf("\t\t\t M = %6i   Largest Discarded Weight = %8.3e  Sweep Energy = %20.10f \n",sweepParams.get_keep_states(), finalError, finalEnergy[j]+dmrginp.get_coreenergy());
#else 
      //printf("\t\t\t M = %6i   Largest Discarded Weight = %8.3e  Sweep Energy = %20.10f \n",sweepParams.get_keep_states(), finalError, finalEnergy[j]+dmrginp.get_coreenergy());
      xout << "\t\t\t M = " <<  sweepParams.get_keep_states() ; 
      xout << "\t Largest Discarded Weight = " << scientific << setprecision(8) << finalError ;
      xout << "\t Sweep Energy = " << fixed << setprecision(10) << finalEnergy[j]+dmrginp.get_coreenergy() << endl;
#endif
    }
  }
  pout << "\t\t\t ============================================================================ " << endl;

  // update the static number of iterations

  ++sweepParams.set_sweep_iter();

  return finalEnergy[0];
}

void SpinAdapted::Sweep::Startup (SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem)
{
  mcheck("at the start of block and decimate");
  dmrginp.guessgenT -> start();
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

  dmrginp.datatransfer -> start();
  system.addAdditionalCompOps();
  dmrginp.datatransfer -> stop();
  InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, sweepParams.get_sys_add(), dmrginp.direct(), 
				 DISTRIBUTED_STORAGE, true, true);


  int nquanta = newSystem.get_stateInfo().quanta.size();
  std::vector<DiagonalMatrix > energies(nquanta);
  std::vector<Matrix> rotateMatrix(nquanta);
  DensityMatrix transformmatrix; 
  transformmatrix.allocate(newSystem.get_stateInfo());
  SpinQuantum q(0,0,IrrepSpace(0));

  if (mpigetrank() == 0) {
    double minval = 1e12;
    for (int i=0; i<nquanta; i++) {
      diagonalise(newSystem.get_op_rep(HAM,q)->operator_element(i,i), energies[i], transformmatrix(i,i));
      for (int j=0; j<energies[i].Nrows(); j++) 
	if (minval > energies[i](j+1))
	  minval = energies[i](j+1);
    }
    for (int i=0; i<nquanta; i++) {
      for (int j=0; j<energies[i].Nrows(); j++) 
	energies[i](j+1) = 1.0/(energies[i](j+1)-minval+1);
    }


    vector<pair<int, int> > inorderwts;
    vector<vector<int> > wtsbyquanta;
    
    sort_weights(energies, inorderwts, wtsbyquanta);
    
    // make transformation matrix by various algorithms
    int keptstates = sweepParams.get_keep_states()/2, keptqstates = sweepParams.get_keep_states()-keptstates;
    int totalstatesbydm = min(static_cast<int>(inorderwts.size()), keptstates);
    int totalstatesbyquanta = min(static_cast<int>(inorderwts.size()), keptstates + keptqstates) - totalstatesbydm;
    if (totalstatesbyquanta < 0) totalstatesbyquanta = 0;
    
    pout << "\t\t\t total states using dm and quanta " << totalstatesbydm << " " << totalstatesbyquanta << endl;
    
    
    double error = assign_matrix_by_dm(rotateMatrix, energies, transformmatrix, inorderwts, wtsbyquanta, totalstatesbydm, totalstatesbyquanta, newSystem.size(), 2*totalstatesbydm);
    pout << "\t\t\t Total discarded weight "<<error<<endl;
  }

#ifndef SERIAL
  mpi::communicator world;
  broadcast(world, rotateMatrix, 0);
#endif

  dmrginp.operrotT -> start();
  newSystem.transform_operators(rotateMatrix);
  for (int i=0; i<dmrginp.nroots(); i++)
    SaveRotationMatrix (newSystem.get_sites(), rotateMatrix, i);
  dmrginp.operrotT -> stop();
  mcheck("after rotation and transformation of block");
  


  pout << dmrginp.guessgenT<<" "<<dmrginp.multiplierT<<" "<<dmrginp.operrotT<< "  "<<globaltimer.totalwalltime()<<" timer "<<endl;
  pout << dmrginp.makeopsT<<" makeops "<<endl;
  pout << dmrginp.datatransfer<<" datatransfer "<<endl;
  //cout << dmrginp.justmultiply<<" just multiply "<<endl;
  //cout << dmrginp.otherrotation<<" "<<dmrginp.spinrotation<<" "<<dmrginp.operrotT<<" rotations time "<<endl; 
  pout <<"oneindexopmult   twoindexopmult   Hc  couplingcoeff"<<endl;  
  pout << dmrginp.oneelecT<<" "<<dmrginp.twoelecT<<" "<<dmrginp.hmultiply<<" "<<dmrginp.couplingcoeff<<" hmult"<<endl;
  pout << dmrginp.buildsumblock<<" "<<dmrginp.buildblockops<<" build block"<<endl;
  pout << "addnoise  S_0_opxop  S_1_opxop   S_2_opxop"<<endl;
  pout << dmrginp.addnoise<<" "<<dmrginp.s0time<<" "<<dmrginp.s1time<<" "<<dmrginp.s2time<<endl;
  

  //mcheck("After renorm transform");
}
