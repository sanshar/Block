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
#include "operatorfunctions.h"

using namespace boost;
using namespace std;


void SpinAdapted::SweepResponse::BlockAndDecimate (SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, const bool &useSlater, const bool& dot_with_sys, int targetState, vector<int>& projectors, vector<int>& baseStates)
{
  if (dmrginp.outputlevel() > 0)
    mcheck("at the start of block and decimate");
  p2out << "\t\t\t dot with system "<<dot_with_sys<<endl;
  p1out <<endl<< "\t\t\t Performing Blocking"<<endl;
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
  if (!dot_with_sys && sweepParams.get_onedot())
    {
      pout << "\t\t\t System  Block"<<system;
    }
  else
    {
      pout << "\t\t\t System  Block"<<newSystem;
    }
  pout << "\t\t\t Environment Block"<<newEnvironment<<endl;
  p1out << "\t\t\t Solving wavefunction "<<endl;

  std::vector<Wavefunction> lowerStates;


  //make the baseState
  int originalOutputlevel = dmrginp.outputlevel();
  dmrginp.setOutputlevel() = -1;
  lowerStates.resize(projectors.size()+1); //only make the 
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

  //<target|O|firstOrderState>
  DensityMatrix branoiseMatrix(newSystem.get_braStateInfo());
  branoiseMatrix.allocate(newSystem.get_braStateInfo());
  
  if (dmrginp.outputlevel() > 0)
    mcheck("before making correction vector");
  for (int l=0; l<baseStates.size(); l++)
  {  
    //now one needs to make |phi_0> = O|psi_0> so that the |phi_0> has the same dimensions as our target state
    int firstOrderState = baseStates[l];

    int perturbationIntegral = l+1;
    SpinBlock perturbationBig, perturbationsystemdot, perturbationenvironmentdot;
    SpinBlock perturbationsystem, perturbationenvironment, perturbationnewsystem, perturbationnewenvironment;
    perturbationsystemdot = SpinBlock(systemDotStart, systemDotEnd, perturbationIntegral, true);
    perturbationenvironmentdot = SpinBlock(environmentDotStart, environmentDotEnd, perturbationIntegral, true);
    perturbationsystem.set_integralIndex() = perturbationIntegral;
    SpinBlock::restore(forward, systemsites, perturbationsystem, targetState, firstOrderState);

    Sweep::makeSystemEnvironmentBigBlocks(perturbationsystem, perturbationsystemdot, 
					  perturbationnewsystem, perturbationenvironment, 
					  perturbationenvironmentdot, perturbationnewenvironment, 
					  perturbationBig, sweepParams, dot_with_sys, useSlater,
					  perturbationIntegral, targetState, firstOrderState);
    
    Wavefunction iwave;
    GuessWave::guess_wavefunctions(iwave, e, perturbationBig, guesstype, 
				   sweepParams.get_onedot(), firstOrderState, dot_with_sys, 0.0);

    //dont add noise in the onedot algorithm
    if (!sweepParams.get_onedot() && sweepParams.get_noise() > NUMERICAL_ZERO && l == 0) { //only add noise using one basestate
      int sweepiter = sweepParams.get_sweep_iter();
      branoiseMatrix.add_onedot_noise_forCompression(iwave, perturbationBig, 1.0/DotProduct(iwave, iwave));
    }

#ifndef SERIAL
    mpi::communicator world;
    broadcast(world, iwave, 0);
#endif

    Wavefunction temp; 
    temp.set_onedot( sweepParams.get_onedot());
    temp.AllowQuantaFor(*perturbationBig.get_braStateInfo().leftStateInfo,*perturbationBig.get_braStateInfo().rightStateInfo, dmrginp.effective_molecule_quantum_vec());
    temp.Clear();


    perturbationBig.multiplyH(iwave, &temp, MAX_THRD);
    if(l==0)
      lowerStates[0] = temp;
    else
      ScaleAdd(1.0, temp, lowerStates[0]);


    perturbationsystem.clear(); perturbationenvironment.clear(); 
    perturbationnewsystem.clear(); perturbationnewenvironment.clear();  
  }
  if (mpigetrank() != 0) {
    lowerStates[0].CleanUp();
  }


  //<target|O|projectors>
  for (int l=0; l<projectors.size(); l++)
  {  
    
    //now one needs to make |phi_0> = O|psi_0> so that the |phi_0> has the same dimensions as our target state
    int perturbationIntegral = 0;
    SpinBlock overlapBig, overlapsystemdot, overlapenvironmentdot;
    SpinBlock overlapsystem, overlapenvironment, overlapnewsystem, overlapnewenvironment;
    overlapsystemdot = SpinBlock(systemDotStart, systemDotEnd, perturbationIntegral, true);
    overlapenvironmentdot = SpinBlock(environmentDotStart, environmentDotEnd, perturbationIntegral, true);
    overlapsystem.set_integralIndex() = perturbationIntegral;
    SpinBlock::restore(forward, systemsites, overlapsystem, targetState, projectors[l]);
    
    Sweep::makeSystemEnvironmentBigOverlapBlocks(systemsites, overlapsystemdot, overlapenvironmentdot, 
						 overlapsystem, overlapnewsystem, overlapenvironment, 
						 overlapnewenvironment, overlapBig, sweepParams, dot_with_sys, useSlater,
						 perturbationIntegral, targetState, projectors[l]);
    
    Wavefunction iwave;
    GuessWave::guess_wavefunctions(iwave, e, overlapBig, guesstype, 
				   sweepParams.get_onedot(), projectors[l], dot_with_sys, 0.0);
    
#ifndef SERIAL
    mpi::communicator world;
    broadcast(world, iwave, 0);
#endif

    Wavefunction temp; 
    temp.set_onedot( sweepParams.get_onedot());
    temp.AllowQuantaFor(*overlapBig.get_braStateInfo().leftStateInfo,*overlapBig.get_braStateInfo().rightStateInfo, dmrginp.effective_molecule_quantum_vec());
    temp.Clear();
    overlapBig.multiplyOverlap(iwave, &temp, MAX_THRD);

    lowerStates[l+1] = temp;
    int success = 0;
    if (mpigetrank() == 0) {
      for (int istate=1; istate<l+1; istate++)  {
	double overlap = pow(DotProduct(lowerStates[istate], lowerStates[istate]), 0.5);
	ScaleAdd(-DotProduct(lowerStates[istate], lowerStates[l+1])/overlap, lowerStates[istate], lowerStates[l+1]);
      }
    
      lowerStates[l+1].Normalise(&success);
    }
    if (mpigetrank() != 0) {
      lowerStates[l+1].CleanUp();
    }

    overlapsystem.clear(); overlapenvironment.clear(); 
    overlapnewsystem.clear(); overlapnewenvironment.clear();  
  }


  if (dmrginp.outputlevel() > 0)
    mcheck("Before renormalization");
  dmrginp.setOutputlevel() = originalOutputlevel;
  DensityMatrix bratracedMatrix;
  newSystem.RenormaliseFrom (sweepParams.set_lowest_energy(), sweepParams.set_lowest_energy_spins(),
			     sweepParams.set_lowest_error(), rotatematrix, 
			     sweepParams.get_keep_states(), sweepParams.get_keep_qstates(), 
			     sweepParams.get_davidson_tol(), big, sweepParams.get_guesstype(), 
			     0.0, 0.0, //noise 
			     sweepParams.get_onedot(), system, systemDot, environment, 
			     dot_with_sys, useSlater, sweepParams.get_sweep_iter(), targetState, 
			     lowerStates, &bratracedMatrix);



  if (dmrginp.outputlevel() > 0)
    mcheck("");
  environment.clear();
  newEnvironment.clear();

  p1out <<"\t\t\t Performing Renormalization "<<endl;

  rotatematrix.resize(0);

  if (dmrginp.outputlevel() > 0)
    mcheck("Before adding noise");
  ScaleAdd(sweepParams.get_noise()*(max(1.e-5, trace(branoiseMatrix))), branoiseMatrix, bratracedMatrix);
  if (!mpigetrank())
    sweepParams.set_lowest_error() = makeRotateMatrix(bratracedMatrix, rotatematrix, sweepParams.get_keep_states(), sweepParams.get_keep_qstates());

  p2out << "discarded weight "<<sweepParams.set_lowest_error()<<endl;

  Wavefunction targetWave; StateInfo braStateInfo;
  targetWave.LoadWavefunctionInfo (braStateInfo, newSystem.get_sites(), targetState);
#ifndef SERIAL
  mpi::communicator world;
  broadcast(world, rotatematrix, 0);
  broadcast(world, targetWave, 0);
#endif

  //<target|O|firstOrderState>
  dmrginp.setOutputlevel() = -1; 
  for (int l=0; l<baseStates.size(); l++) 
  {
    int firstOrderState = baseStates[l];
    int originalOutputlevel = dmrginp.outputlevel();
    //dmrginp.setOutputlevel() = -1;

    int perturbationIntegral = l+1;
    SpinBlock perturbationBig, perturbationsystemdot, perturbationenvironmentdot;
    SpinBlock perturbationsystem, perturbationenvironment, perturbationnewsystem, perturbationnewenvironment;
    perturbationsystemdot = SpinBlock(systemDotStart, systemDotEnd, perturbationIntegral, true);
    perturbationenvironmentdot = SpinBlock(environmentDotStart, environmentDotEnd, perturbationIntegral, true);
    perturbationsystem.set_integralIndex() = perturbationIntegral;
    SpinBlock::restore(forward, systemsites, perturbationsystem, targetState, firstOrderState);

    Sweep::makeSystemEnvironmentBigBlocks(perturbationsystem, perturbationsystemdot, 
					  perturbationnewsystem, perturbationenvironment, 
					  perturbationenvironmentdot, perturbationnewenvironment, 
					  perturbationBig, sweepParams, dot_with_sys, useSlater,
					  perturbationIntegral, targetState, firstOrderState);

    SpinBlock perturbationnewbig;
    if (sweepParams.get_onedot() && !dot_with_sys)
    {
      InitBlocks::InitNewSystemBlock(perturbationsystem, perturbationsystemdot, 
				     perturbationnewsystem, targetState, firstOrderState, 
				     perturbationsystemdot.size(), dmrginp.direct(), 
				     perturbationIntegral, DISTRIBUTED_STORAGE, false, true);
      InitBlocks::InitBigBlock(perturbationnewsystem, perturbationenvironment, perturbationnewbig); 

      perturbationBig.get_rightBlock()->clear();
      perturbationBig.clear();
    }
    else
      perturbationnewbig = perturbationBig;

    std::vector<Matrix> ketrotatematrix;
    if (mpigetrank() == 0) {
      Wavefunction iwave;
      GuessWave::guess_wavefunctions(iwave, e, perturbationnewbig, guesstype, 
				     sweepParams.get_onedot(), firstOrderState, true, 0.0);
      
      DensityMatrix tracedMatrix;
      tracedMatrix.allocate(perturbationnewbig.get_leftBlock()->get_ketStateInfo());
      operatorfunctions::MultiplyProduct(iwave, Transpose(const_cast<Wavefunction&> (iwave)), tracedMatrix, 1.0);
      int largeNumber = 1000000;
      if (!mpigetrank())
	double error = makeRotateMatrix(tracedMatrix, ketrotatematrix, largeNumber, sweepParams.get_keep_qstates());

      iwave.SaveWavefunctionInfo (perturbationnewbig.get_ketStateInfo(), perturbationnewbig.get_leftBlock()->get_sites(), firstOrderState);

      if ( l == 0)
	SaveRotationMatrix (perturbationnewsystem.get_sites(), ketrotatematrix);
      SaveRotationMatrix (perturbationnewsystem.get_sites(), ketrotatematrix, firstOrderState);
      
    }

#ifndef SERIAL
    broadcast(world, ketrotatematrix, 0);
#endif
        
    perturbationnewsystem.transform_operators(rotatematrix, ketrotatematrix);

    SpinBlock::store(forward, perturbationnewsystem.get_sites(), perturbationnewsystem, targetState, firstOrderState);


    perturbationsystem.clear(); perturbationenvironment.clear(); perturbationnewsystem.clear(); perturbationnewenvironment.clear();

  }


  //<target|O|projectors>
  for (int l=0; l<projectors.size(); l++)
  {
    int originalOutputlevel = dmrginp.outputlevel();
    //dmrginp.setOutputlevel() = -1;

    int perturbationIntegral = 0;
    SpinBlock perturbationBig, perturbationsystemdot, perturbationenvironmentdot;
    SpinBlock perturbationsystem, perturbationenvironment, perturbationnewsystem, perturbationnewenvironment;
    perturbationsystemdot = SpinBlock(systemDotStart, systemDotEnd, perturbationIntegral, true);
    perturbationenvironmentdot = SpinBlock(environmentDotStart, environmentDotEnd, perturbationIntegral, true);
    perturbationsystem.set_integralIndex() = perturbationIntegral;
    SpinBlock::restore(forward, systemsites, perturbationsystem, targetState, projectors[l]);

    Sweep::makeSystemEnvironmentBigOverlapBlocks(systemsites, perturbationsystemdot, perturbationenvironmentdot, 
						 perturbationsystem, perturbationnewsystem, perturbationenvironment, perturbationnewenvironment, 
						 perturbationBig, sweepParams, dot_with_sys, useSlater,
						 perturbationIntegral, targetState, projectors[l]);

    SpinBlock perturbationnewbig;
    if (sweepParams.get_onedot() && !dot_with_sys)
    {
      perturbationnewsystem.set_integralIndex() = perturbationsystem.get_integralIndex();
      perturbationnewsystem.initialise_op_array(OVERLAP, false);
      perturbationnewsystem.setstoragetype(DISTRIBUTED_STORAGE);
      perturbationnewsystem.BuildSumBlock(NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, perturbationsystem, perturbationsystemdot);
      InitBlocks::InitBigBlock(perturbationnewsystem, perturbationenvironment, perturbationnewbig); 

      perturbationBig.get_rightBlock()->clear();
      perturbationBig.clear();
    }
    else
      perturbationnewbig = perturbationBig;

    std::vector<Matrix> ketrotatematrix;
    if (mpigetrank() == 0) {
      Wavefunction iwave;
      GuessWave::guess_wavefunctions(iwave, e, perturbationnewbig, guesstype, 
				     sweepParams.get_onedot(), projectors[l], true, 0.0);
      
      dmrginp.setOutputlevel() = 10;
      DensityMatrix tracedMatrix;
      tracedMatrix.allocate(perturbationnewbig.get_leftBlock()->get_ketStateInfo());
      operatorfunctions::MultiplyProduct(iwave, Transpose(const_cast<Wavefunction&> (iwave)), tracedMatrix, 1.0);
      int largeNumber = 1000000;
      if (!mpigetrank())
	double error = makeRotateMatrix(tracedMatrix, ketrotatematrix, largeNumber, sweepParams.get_keep_qstates());

      iwave.SaveWavefunctionInfo (perturbationnewbig.get_ketStateInfo(), perturbationnewbig.get_leftBlock()->get_sites(), projectors[l]);
      SaveRotationMatrix (perturbationnewsystem.get_sites(), ketrotatematrix, projectors[l]);
    }

#ifndef SERIAL
    broadcast(world, ketrotatematrix, 0);
#endif
    
    
    if (dmrginp.outputlevel() > 0)
      mcheck("transform operators");
    perturbationnewsystem.transform_operators(rotatematrix, ketrotatematrix);

    SpinBlock::store(forward, perturbationnewsystem.get_sites(), perturbationnewsystem, targetState, projectors[l]);


    perturbationsystem.clear(); perturbationenvironment.clear(); perturbationnewsystem.clear(); perturbationnewenvironment.clear();

  }

#ifndef SERIAL
    broadcast(world, rotatematrix, 0);
#endif


  dmrginp.setOutputlevel() = originalOutputlevel;

  dmrginp.operrotT -> start();
  newSystem.transform_operators(rotatematrix);
  SaveRotationMatrix (newSystem.get_sites(), rotatematrix, targetState);
  dmrginp.operrotT -> stop();

  if (dmrginp.outputlevel() > 0)
    mcheck("after rotation and transformation of block");

  p2out << *dmrginp.guessgenT<<" "<<*dmrginp.multiplierT<<" "<<*dmrginp.operrotT<< "  "<<globaltimer.totalwalltime()<<" timer "<<endl;
  p2out << *dmrginp.makeopsT<<" makeops "<<endl;
  p2out << *dmrginp.datatransfer<<" datatransfer "<<endl;
  p2out <<"oneindexopmult   twoindexopmult   Hc  couplingcoeff"<<endl;  
  p2out << *dmrginp.oneelecT<<" "<<*dmrginp.twoelecT<<" "<<*dmrginp.hmultiply<<" "<<*dmrginp.couplingcoeff<<" hmult"<<endl;
  p2out << *dmrginp.buildsumblock<<" "<<*dmrginp.buildblockops<<" build block"<<endl;
  p2out << *dmrginp.blockintegrals<<"  "<<*dmrginp.blocksites<<"  "<<*dmrginp.statetensorproduct<<"  "<<*dmrginp.statecollectquanta<<"  "<<*dmrginp.buildsumblock<<" "<<*dmrginp.buildblockops<<" build sum block"<<endl;
  p2out << "addnoise  S_0_opxop  S_1_opxop   S_2_opxop"<<endl;
  p3out << *dmrginp.addnoise<<" "<<*dmrginp.s0time<<" "<<*dmrginp.s1time<<" "<<*dmrginp.s2time<<endl;
  

}



double SpinAdapted::SweepResponse::do_one(SweepParams &sweepParams, const bool &warmUp, const bool &forward, 
					    const bool &restart, const int &restartSize, int targetState, 
					  vector<int>& projectors, vector<int>& baseStates, int correctionVector)
{
  SpinBlock system;
  int activeSpaceIntegral = 0;
  std::vector<int> perturbationIntegral(dmrginp.getNumIntegrals()-1,0);
  for (int i=0; i<perturbationIntegral.size(); i++)
    perturbationIntegral[i] = i+1;
  std::vector<double> finalEnergy(1,0.0e100);
  double finalError = 0.;

  if (restart) {
    finalEnergy = sweepParams.get_lowest_energy();
    finalError = sweepParams.get_lowest_error();
  }

  sweepParams.set_sweep_parameters();
  // a new renormalisation sweep routine
  pout << endl;
  if (forward)
    {
      pout << "\t\t\t Starting sweep "<< sweepParams.set_sweep_iter()<<" in forwards direction"<<endl;
    }
  else
    {
      pout << "\t\t\t Starting sweep "<< sweepParams.set_sweep_iter()<<" in backwards direction" << endl;
    }
  pout << "\t\t\t ============================================================================ " << endl;

  InitBlocks::InitStartingBlock (system,forward, targetState, targetState,
				 sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), 
				 restartSize, restart, warmUp, activeSpaceIntegral);

  for (int l=0; l<projectors.size(); l++)
  {
    SpinBlock perturbationSystem;
    perturbationSystem.set_integralIndex() = 0;
    InitBlocks::InitStartingBlock (perturbationSystem,forward, targetState, projectors[l],
				   sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), 
				   restartSize, restart, warmUp, 0);
    SpinBlock::store (forward, system.get_sites(), perturbationSystem, targetState, projectors[l]);
  }

  for (int l=0; l<baseStates.size(); l++)
  {
    SpinBlock overlapSystem;
    overlapSystem.set_integralIndex() = l+1;
    InitBlocks::InitStartingBlock (overlapSystem,forward, targetState, baseStates[l],
				   sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), 
				   restartSize, restart, warmUp, perturbationIntegral[l]);
    SpinBlock::store (forward, system.get_sites(), overlapSystem, targetState, baseStates[l]);
  }

  if(!restart)
    sweepParams.set_block_iter() = 0;

 
  p2out << "\t\t\t Starting block is :: " << endl << system << endl;

  // if restart, just restoring an existing block --
  SpinBlock::store (forward, system.get_sites(), system, targetState, targetState);
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

 // get_n_iters() returns the number of blocking iterations needed in one sweep
  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); )
    {
      pout << "\n\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
      pout << "\t\t\t ----------------------------" << endl;
      if (forward)
	{
	  p1out << "\t\t\t Current direction is :: Forwards " << endl;
	}
      else
	{
	  p1out << "\t\t\t Current direction is :: Backwards " << endl;
	}


      if (sweepParams.get_block_iter() == 0 && sweepParams.get_sweep_iter() == 1)
	{
	  sweepParams.set_guesstype() = BASIC;
	}
      else if (sweepParams.get_block_iter() != 0)
	{
	  sweepParams.set_guesstype() = TRANSFORM;
	}
      else
	{
	  sweepParams.set_guesstype() = TRANSPOSE;
	}


      
      p1out << "\t\t\t Blocking and Decimating " << endl;
	  
      SpinBlock newSystem; // new system after blocking and decimating

      //Need to substitute by:
      if (warmUp ) {
	int correctionVector = dmrginp.guessState();
	p2out << "USING state "<<correctionVector<<" as initial guess"<<endl;

	StartUp(sweepParams, system, newSystem, dot_with_sys, 
		targetState, correctionVector, projectors, baseStates);
      }
      else {
	double cE = 0.0;
	if (dmrginp.calc_type() == RESPONSEBW) {
	  //When doing BW perturbation theory in the denominator we have (H0 - E) and not (H0-E0)
	  cE = coreEnergy[system.get_integralIndex()];
	  double e2 = BWPTenergy;
	  coreEnergy[system.get_integralIndex()] = cE - e2;
	  p2out << "\t\t\t BW perturbation theory  "<<cE<<"  "<<e2<<endl; 
	}

	BlockAndDecimate(sweepParams, system, newSystem, warmUp, 
			 dot_with_sys, targetState, projectors, baseStates);

	if (dmrginp.calc_type() == RESPONSEBW) 
	  coreEnergy[system.get_integralIndex()] = cE ;

      }
      
      //Need to substitute by?

      if (!warmUp ){

	//this criteria should work for state average or state specific because the lowest sweep energy is always the lowest of the average
	finalError = max(sweepParams.get_lowest_error(),finalError);
	finalEnergy[0] = min(sweepParams.get_lowest_energy()[0], finalEnergy[0]);
	BWPTenergy = finalEnergy[0];
	pout << "final energy "<<finalEnergy[0]<<"  "<<sweepParams.get_lowest_energy()[0]<<endl;
      }
      
      system = newSystem;
      p2out << system<<endl;
      system.printOperatorSummary();
      
      //system size is going to be less than environment size
      if (forward && system.get_complementary_sites()[0] >= dmrginp.last_site()/2)
	    dot_with_sys = false;
      if (!forward && system.get_sites()[0]-1 < dmrginp.last_site()/2)
	    dot_with_sys = false;

      SpinBlock::store (forward, system.get_sites(), system, targetState, targetState);
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
	{
	  p1out << "\t\t\t Current direction is :: Forwards " << endl;
	}
      else
	{
	  p1out << "\t\t\t Current direction is :: Backwards " << endl;
	}
    sweepParams.set_onedot() = true;
    sweepParams.set_env_add() = 0;
    bool dot_with_sys = true;
    WavefunctionCanonicalize(sweepParams, system, warmUp, dot_with_sys, targetState, projectors, baseStates);
    sweepParams.set_onedot() = false;
    sweepParams.set_env_add() = 1;
  }

  //pout << "\t\t\t Largest Error for Sweep with " << sweepParams.get_keep_states() << " states is " << finalError << endl;
  //pout << "\t\t\t Sweep Energy for Sweep with " << sweepParams.get_keep_states() << " states is " << finalEnergy[0] << endl;
  sweepParams.set_largest_dw() = finalError;
  
  if (mpigetrank() == 0)
    printf("\t\t\t M = %6i  state = %4i  Largest Discarded Weight = %8.3e  Sweep Energy = %20.10e \n",sweepParams.get_keep_states(), 0, finalError, finalEnergy[0]);

  pout << "\t\t\t ============================================================================ " << endl;

  // update the static number of iterations

  ++sweepParams.set_sweep_iter();

  if (mpigetrank()==0)
  {
    pout << "About to write dmrg energy"<<endl;
    std::string efile;
    efile = str(boost::format("%s%s") % dmrginp.load_prefix() % "/dmrg.e" );
    
    
    FILE* f = fopen(efile.c_str(), "wb");      
    double e = finalEnergy[0]; //sweepParams.get_lowest_energy()[0]; //instead of the lowest energy of the sweep, we record the last energy of the sweep
    fwrite( &e, 1, sizeof(double), f);
    fclose(f);
  }



  return finalError;
}


void SpinAdapted::SweepResponse::StartUp (SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, const bool& dot_with_sys, int targetState, int correctionVector, vector<int>& projectors, vector<int>& baseStates)
{
  if (dmrginp.outputlevel() > 0)
    mcheck("at the start of block and decimate");
  p2out << "\t\t\t dot with system "<<dot_with_sys<<endl;
  p1out <<endl<< "\t\t\t Performing Blocking"<<endl;
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
    system.addAdditionalCompOps();
    InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, correctionVector, correctionVector, 
				   sweepParams.get_sys_add(), dmrginp.direct(), system.get_integralIndex(), 
				   DISTRIBUTED_STORAGE, haveNormOps, haveCompOps);
    
    dmrginp.guessgenT -> stop();
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
  
  //<target| V |baseStates>
  for(int l=0; l<baseStates.size(); l++) 
  {
    dmrginp.guessgenT -> start();
    int perturbationIntegral = l+1;
    SpinBlock perturbationSystemDot, perturbationSystem, perturbationNewSystem;
    perturbationSystem.set_integralIndex() = perturbationIntegral;
    SpinBlock::restore(forward, systemsites, perturbationSystem, targetState, baseStates[l]);
    perturbationSystemDot = SpinBlock(systemDotStart, systemDotEnd, perturbationIntegral, true);

    perturbationSystem.addAdditionalCompOps();
    InitBlocks::InitNewSystemBlock(perturbationSystem, perturbationSystemDot, perturbationNewSystem,
				   targetState, baseStates[l], 
				   sweepParams.get_sys_add(), dmrginp.direct(), perturbationIntegral,
				   DISTRIBUTED_STORAGE, haveNormOps, haveCompOps);

    
    dmrginp.guessgenT -> stop();

    LoadRotationMatrix (newSystem.get_sites(), ketrotateMatrix, baseStates[l]);
    
#ifndef SERIAL
    mpi::communicator world;
    broadcast(world, ketrotateMatrix, 0);
#endif

    dmrginp.operrotT -> start();
    perturbationNewSystem.transform_operators(brarotateMatrix, ketrotateMatrix);
    dmrginp.operrotT -> stop();
    SpinBlock::store(forward, perturbationNewSystem.get_sites(), perturbationNewSystem, targetState, baseStates[l]);
  }


  //<target|O|baseState>
  //save the updated overlap spinblock
  for(int l=0; l<projectors.size(); l++) 
  {
    int perturbationIntegral = 0;
    
    SpinBlock overlapsystem, overlapsystemDot, overlapnewSystem;
    overlapsystem.set_integralIndex() = perturbationIntegral;
    SpinBlock::restore(forward, systemsites, overlapsystem, targetState, projectors[l]);
    overlapsystemDot = SpinBlock(systemDotStart, systemDotEnd, perturbationIntegral, true);

    overlapnewSystem.set_integralIndex() = perturbationIntegral;
    overlapnewSystem.initialise_op_array(OVERLAP, false);
    overlapnewSystem.setstoragetype(DISTRIBUTED_STORAGE);
    overlapnewSystem.BuildSumBlock(NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, overlapsystem, overlapsystemDot);

    LoadRotationMatrix (newSystem.get_sites(), ketrotateMatrix, projectors[l]);
    
#ifndef SERIAL
    mpi::communicator world;
    broadcast(world, ketrotateMatrix, 0);
#endif

    dmrginp.operrotT -> start();
    overlapnewSystem.transform_operators(brarotateMatrix, ketrotateMatrix);
    dmrginp.operrotT -> stop();
    SpinBlock::store(forward, overlapnewSystem.get_sites(), overlapnewSystem, targetState, projectors[l]);

  }



  if (dmrginp.outputlevel() > 0)
    mcheck("after rotation and transformation of block");

  p2out << *dmrginp.guessgenT<<" "<<*dmrginp.multiplierT<<" "<<*dmrginp.operrotT<< "  "<<globaltimer.totalwalltime()<<" timer "<<endl;
  p2out << *dmrginp.makeopsT<<" makeops "<<endl;
  p2out << *dmrginp.datatransfer<<" datatransfer "<<endl;
  p2out <<"oneindexopmult   twoindexopmult   Hc  couplingcoeff"<<endl;  
  p2out << *dmrginp.oneelecT<<" "<<*dmrginp.twoelecT<<" "<<*dmrginp.hmultiply<<" "<<*dmrginp.couplingcoeff<<" hmult"<<endl;
  p2out << *dmrginp.buildsumblock<<" "<<*dmrginp.buildblockops<<" build block"<<endl;
  p2out << *dmrginp.blockintegrals<<"  "<<*dmrginp.blocksites<<"  "<<*dmrginp.statetensorproduct<<"  "<<*dmrginp.statecollectquanta<<"  "<<*dmrginp.buildsumblock<<" "<<*dmrginp.buildblockops<<" build sum block"<<endl;
  p2out << "addnoise  S_0_opxop  S_1_opxop   S_2_opxop"<<endl;
  p3out << *dmrginp.addnoise<<" "<<*dmrginp.s0time<<" "<<*dmrginp.s1time<<" "<<*dmrginp.s2time<<endl;

}


void SpinAdapted::SweepResponse::WavefunctionCanonicalize (SweepParams &sweepParams, SpinBlock& system, const bool &useSlater, const bool& dot_with_sys, int targetState, vector<int>& projectors, vector<int>& baseStates)
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

  SpinBlock::restore(!forward, sitesenvdot, environmentDot, targetState, targetState); 

  SpinBlock environment, newEnvironment;
  
  SpinBlock big;  // new_sys = sys+sys_dot; new_env = env+env_dot; big = new_sys + new_env then renormalize to find new_sys(new)

  system.addAdditionalCompOps();
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
  if (!dot_with_sys && sweepParams.get_onedot())
    {
      pout << "\t\t\t System  Block"<<system;
    }
  else
    {
      pout << "\t\t\t System  Block"<<newSystem;
    }
  pout << "\t\t\t Environment Block"<<newEnvironment<<endl;
  p1out << "\t\t\t Solving wavefunction "<<endl;
  
  std::vector<Wavefunction> lowerStates;
  
  
  //make the baseState
  int originalOutputlevel = dmrginp.outputlevel();
  dmrginp.setOutputlevel() = -1;
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
  pout << "transform previous"<<endl;
  if (!mpigetrank())
    GuessWave::transform_previous_twodot_to_onedot_wavefunction(targetWave, big, targetState);

#ifndef SERIAL
    mpi::communicator world;
    broadcast(world, targetWave, 0);
#endif

  
  DensityMatrix bratracedMatrix;
  bratracedMatrix.allocate(newSystem.get_braStateInfo());
  operatorfunctions::MultiplyProduct(targetWave, Transpose(const_cast<Wavefunction&> (targetWave)), bratracedMatrix, 1.0);
  int largeNumber = 1000000;
  if (!mpigetrank())
    double error = makeRotateMatrix(bratracedMatrix, rotatematrix, largeNumber, sweepParams.get_keep_qstates());
  
  pout << "broadcast rotate"<<endl;
#ifndef SERIAL
  broadcast(world, rotatematrix, 0);
#endif

  
  targetWave.SaveWavefunctionInfo (big.get_braStateInfo(), big.get_leftBlock()->get_sites(), targetState);
  SaveRotationMatrix (newSystem.get_sites(), rotatematrix, targetState);
  
  newSystem.transform_operators(rotatematrix);
  
  
  dmrginp.setOutputlevel() = originalOutputlevel;
  
  
  
  
  
  
  //<target|O|firstOrderState>
  //save the updated overlap spinblock
  for (int l=0; l<baseStates.size(); l++)
  {
    int perturbationIntegral = l+1;
    
    SpinBlock overlapBig;
    SpinBlock overlapsystem, overlapenvironment, overlapnewsystem, overlapnewenvironment, overlapenvironmentDot;
    SpinBlock overlapsystemDot(systemDotStart, systemDotEnd, perturbationIntegral, true);

    overlapenvironmentDot.set_integralIndex() = perturbationIntegral;
    SpinBlock::restore(!forward, sitesenvdot, overlapenvironmentDot, targetState, baseStates[l]); 

    guessWaveTypes guesstype = sweepParams.get_block_iter() == 0 ? TRANSPOSE : TRANSFORM;
    
    DiagonalMatrix e;

    overlapsystem.set_integralIndex() = perturbationIntegral;
    SpinBlock::restore(forward, systemsites, overlapsystem, targetState, baseStates[l]);
    overlapsystem.addAdditionalCompOps();
    InitBlocks::InitNewSystemBlock(overlapsystem, overlapsystemDot, 
				   overlapnewsystem, targetState, baseStates[l], 
				   overlapsystemDot.size(), dmrginp.direct(), 
				   perturbationIntegral, DISTRIBUTED_STORAGE, false, true);

    overlapnewsystem.set_loopblock(false);  overlapenvironmentDot.set_loopblock(false); 
    InitBlocks::InitBigBlock(overlapnewsystem, overlapenvironmentDot, overlapBig);

    
    pout << "transform iwave"<<endl;
    Wavefunction iwave;
    if (!mpigetrank())
      GuessWave::transform_previous_twodot_to_onedot_wavefunction(iwave, overlapBig, baseStates[l]);

#ifndef SERIAL
    mpi::communicator world;
    broadcast(world, iwave, 0);
#endif

    std::vector<Matrix> ketrotatematrix;
    DensityMatrix tracedMatrix;
    tracedMatrix.allocate(overlapnewsystem.get_ketStateInfo());
    operatorfunctions::MultiplyProduct(iwave, Transpose(const_cast<Wavefunction&> (iwave)), tracedMatrix, 1.0);
    int largeNumber = 1000000;
    if (!mpigetrank())
      double error = makeRotateMatrix(tracedMatrix, ketrotatematrix, largeNumber, sweepParams.get_keep_qstates());
    
    pout << "broadcast ketrot"<<endl;
#ifndef SERIAL
    broadcast(world, ketrotatematrix, 0);
#endif
    
    iwave.SaveWavefunctionInfo (overlapBig.get_ketStateInfo(), overlapBig.get_leftBlock()->get_sites(), baseStates[l]);
    SaveRotationMatrix (overlapnewsystem.get_sites(), ketrotatematrix, baseStates[l]);
    
    overlapnewsystem.transform_operators(rotatematrix, ketrotatematrix);
    SpinBlock::store(forward, overlapnewsystem.get_sites(), overlapnewsystem, targetState, baseStates[l]);
    overlapsystem.clear(); overlapenvironment.clear(); overlapnewsystem.clear(); overlapnewenvironment.clear();
    
  }


  //<target|O|projectors>
  //save the updated overlap spinblock
  for (int l=0; l<projectors.size(); l++)
  {
    int perturbationIntegral = 0;
    
    SpinBlock overlapBig;
    SpinBlock overlapsystem, overlapenvironment, overlapnewsystem, overlapnewenvironment, overlapenvironmentDot;
    SpinBlock overlapsystemDot(systemDotStart, systemDotEnd, perturbationIntegral, true);

    overlapenvironmentDot.set_integralIndex() = perturbationIntegral;
    SpinBlock::restore(!forward, sitesenvdot, overlapenvironmentDot, targetState, projectors[l]); 

    guessWaveTypes guesstype = sweepParams.get_block_iter() == 0 ? TRANSPOSE : TRANSFORM;
    
    DiagonalMatrix e;

    overlapsystem.set_integralIndex() = perturbationIntegral;
    SpinBlock::restore(forward, systemsites, overlapsystem, targetState, projectors[l]);
    overlapnewsystem.set_integralIndex() = perturbationIntegral;
    overlapnewsystem.initialise_op_array(OVERLAP, false);
    overlapnewsystem.setstoragetype(DISTRIBUTED_STORAGE);
    overlapnewsystem.BuildSumBlock(NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, overlapsystem, overlapsystemDot);

    overlapnewsystem.set_loopblock(false);  overlapenvironmentDot.set_loopblock(false); 
    InitBlocks::InitBigBlock(overlapnewsystem, overlapenvironmentDot, overlapBig);

    
    pout << "transform iwave"<<endl;
    Wavefunction iwave;
    if (!mpigetrank())
      GuessWave::transform_previous_twodot_to_onedot_wavefunction(iwave, overlapBig, projectors[l]);

#ifndef SERIAL
    mpi::communicator world;
    broadcast(world, iwave, 0);
#endif

    std::vector<Matrix> ketrotatematrix;
    DensityMatrix tracedMatrix;
    tracedMatrix.allocate(overlapnewsystem.get_ketStateInfo());
    operatorfunctions::MultiplyProduct(iwave, Transpose(const_cast<Wavefunction&> (iwave)), tracedMatrix, 1.0);
    int largeNumber = 1000000;
    if (!mpigetrank())
      double error = makeRotateMatrix(tracedMatrix, ketrotatematrix, largeNumber, sweepParams.get_keep_qstates());
    
    pout << "broadcast ketrot"<<endl;
#ifndef SERIAL
    broadcast(world, ketrotatematrix, 0);
#endif
    
    iwave.SaveWavefunctionInfo (overlapBig.get_ketStateInfo(), overlapBig.get_leftBlock()->get_sites(), projectors[l]);
    SaveRotationMatrix (overlapnewsystem.get_sites(), ketrotatematrix, projectors[l]);
    
    overlapnewsystem.transform_operators(rotatematrix, ketrotatematrix);
    SpinBlock::store(forward, overlapnewsystem.get_sites(), overlapnewsystem, targetState, projectors[l]);
    overlapsystem.clear(); overlapenvironment.clear(); overlapnewsystem.clear(); overlapnewenvironment.clear();
    
  }

  dmrginp.setOutputlevel() = originalOutputlevel;
  
  pout << newSystem<<endl;
  SpinBlock::store(forward, newSystem.get_sites(), newSystem, targetState, targetState);
  
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
