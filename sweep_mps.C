/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "guess_wavefunction.h"
#include "sweepgenblock.h"
#include "sweep.h"
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

using namespace boost;
using namespace std;


void getComplementarySites(std::vector<int> &sites, std::vector<int> &complementarySites) 
{
  for (int i=0; i<dmrginp.last_site(); ++i)
    if (find(sites.begin(), sites.end(), i) == sites.end())
      complementarySites.push_back(i);
  return;
}

//Canonicalize wavefunction, takes the wavefunction and does a sweep to update all the roation matrices so that we get a consistent wavefunction along the whole sweep
void SpinAdapted::Sweep::CanonicalizeWavefunction(SweepParams &sweepParams, const bool &forward, int currentstate)
{

  sweepParams.set_sweep_parameters();
  sweepParams.set_block_iter() = 0;

  std::vector<int> sites;
  int new_site, wave_site;
  if (forward) {
    pout << "\t\t\t Starting sweep "<< sweepParams.set_sweep_iter()<<" in forwards direction"<<endl;
    new_site = 0;
  }
  else {
    pout << "\t\t\t Starting sweep "<< sweepParams.set_sweep_iter()<<" in backwards direction" << endl;
    new_site = dmrginp.spinAdapted() ? dmrginp.last_site()-1 : dmrginp.last_site()/2-1;
  }
  pout << "\t\t\t ============================================================================ " << 
    endl;

  if (dmrginp.spinAdapted())
    sites.push_back(new_site);
  else {
    sites.push_back(2*new_site);
    sites.push_back(2*new_site+1);
    std::sort(sites.begin(), sites.end());
  }
    
    
  //only need statinfos
  StateInfo stateInfo1; makeStateInfo(stateInfo1, new_site);
  
  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); ) {
      
    pout << "\n\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
    pout << "\t\t\t ----------------------------" << endl;
    
    if (forward) {
      new_site++;
      wave_site = new_site+1;
      p1out << "\t\t\t Current direction is :: Forwards " << endl;
    }
    else {
      new_site--;
      wave_site = new_site-1;
      p1out << "\t\t\t Current direction is :: Backwards " << endl;
    }
    std::vector<int> complementarySites, spindotsites(1, new_site), oldsites = sites, oldcomplement;

    if (dmrginp.spinAdapted())
      sites.push_back(new_site);
    else {
      sites.push_back(2*new_site);
      sites.push_back(2*new_site+1);
      std::sort(sites.begin(), sites.end());
    }

    getComplementarySites(sites, complementarySites);
    getComplementarySites(oldsites, oldcomplement);
    
    StateInfo siteState, newState1, bigstate, envstate; 
    makeStateInfo(siteState, new_site);
    TensorProduct(stateInfo1, siteState, newState1, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
    newState1.CollectQuanta();

    Wavefunction w; w.set_deltaQuantum() = dmrginp.effective_molecule_quantum_vec();
    w.set_onedot(true);

    if (!dmrginp.spinAdapted()) {
      std::vector<int> spinSites(complementarySites.size()/2, 0);
      for (int s=0; s<spinSites.size(); s++)
	spinSites[s] = complementarySites[2*s]/2;
      StateInfo::restore(!forward, spinSites, envstate, currentstate);
    }
    else
      StateInfo::restore(!forward, complementarySites, envstate, currentstate);

    TensorProduct(newState1, envstate, bigstate, PARTICLE_SPIN_NUMBER_CONSTRAINT);

    if (sweepParams.get_block_iter() == 0) 
      GuessWave::transpose_previous_wavefunction(w, bigstate, complementarySites, spindotsites, currentstate, true, true);
    else 
      GuessWave::transform_previous_wavefunction(w, bigstate, oldsites, oldcomplement, currentstate, true, true);
    
    w.SaveWavefunctionInfo(bigstate, sites, currentstate);

      
    //make the newstate
    std::vector<Matrix> rotation1; 
      
      
    DensityMatrix tracedMatrix;
    tracedMatrix.allocate(*bigstate.leftStateInfo);
    operatorfunctions::MultiplyProduct(w, Transpose(const_cast<Wavefunction&> (w)), tracedMatrix, 1.0);
    int largeNumber = 1000000;
    if (!mpigetrank())
      double error = makeRotateMatrix(tracedMatrix, rotation1, largeNumber, sweepParams.get_keep_qstates());
    SaveRotationMatrix (sites, rotation1, currentstate);
    
    StateInfo renormState1;
    SpinAdapted::StateInfo::transform_state(rotation1, newState1, renormState1);
    StateInfo::store(forward, sites, renormState1, currentstate);
    stateInfo1 = renormState1;
    ++sweepParams.set_block_iter();
  }
  
}


//Initialize stateinfo using the rotation matrices
void SpinAdapted::Sweep::InitializeStateInfo(SweepParams &sweepParams, const bool &forward, int currentstate)
{
  sweepParams.set_sweep_parameters();
  sweepParams.set_block_iter() = 0;

  std::vector<int> sites;
  int new_site, wave_site;
  if (forward) 
    new_site = 0;
  else 
    new_site = dmrginp.spinAdapted() ? dmrginp.last_site()-1 : dmrginp.last_site()/2-1;

  if (dmrginp.spinAdapted())
    sites.push_back(new_site);
  else {
    sites.push_back(2*new_site);
    sites.push_back(2*new_site+1);
    std::sort(sites.begin(), sites.end());
  }
    
  //only need statinfos
  StateInfo stateInfo1; makeStateInfo(stateInfo1, new_site);
  StateInfo::store(forward, sites, stateInfo1, currentstate);
  
  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); ) {
      
    if (forward) 
      new_site++;
    else 
      new_site--;

    if (dmrginp.spinAdapted())
      sites.push_back(new_site);
    else {
      sites.push_back(2*new_site);
      sites.push_back(2*new_site+1);
      std::sort(sites.begin(), sites.end());
    }

    
    StateInfo siteState, newState1; 
    makeStateInfo(siteState, new_site);
    TensorProduct(stateInfo1, siteState, newState1, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
    newState1.CollectQuanta();
  
      
    //make the newstate
    std::vector<Matrix> rotation1; 
      
    LoadRotationMatrix (sites, rotation1, currentstate);
    
    StateInfo renormState1;
    SpinAdapted::StateInfo::transform_state(rotation1, newState1, renormState1);
    StateInfo::store(forward, sites, renormState1, currentstate);
    stateInfo1 = renormState1;
    ++sweepParams.set_block_iter();
  }
  
}


//before you start optimizing each state you want to initalize all the overlap matrices
void Sweep::InitializeOverlapSpinBlocks(SweepParams &sweepParams, const bool &forward, int stateA, int stateB, int integralIndex)
{
  SpinBlock system;

  sweepParams.set_sweep_parameters();
  if (forward)
    { pout << "\t\t\t Starting sweep "<< sweepParams.set_sweep_iter()<<" in forwards direction"<<endl; }
  else
    { pout << "\t\t\t Starting sweep "<< sweepParams.set_sweep_iter()<<" in backwards direction" << endl; }
  pout << "\t\t\t ============================================================================ " << endl;

  int restartSize = 0; bool restart = false, warmUp = false;
  InitBlocks::InitStartingBlock (system,forward, stateA, stateB, sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), restartSize, restart, warmUp, integralIndex);

  sweepParams.set_block_iter() = 0;

 
  p2out << "\t\t\t Starting block is :: " << endl << system << endl;

  SpinBlock::store (forward, system.get_sites(), system, stateA, stateB); // if restart, just restoring an existing block --
  sweepParams.savestate(forward, system.get_sites().size());
  bool dot_with_sys = true;
  vector<int> syssites = system.get_sites();

  if (dmrginp.outputlevel() > 0)
    mcheck("at the very start of sweep");  // just timer

  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); ) // get_n_iters() returns the number of blocking iterations needed in one sweep
    {
      pout << "\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
      pout << "\t\t\t ----------------------------" << endl;
      if (forward)  { p1out << "\t\t\t Current direction is :: Forwards " << endl; }
      else  { p1out << "\t\t\t Current direction is :: Backwards " << endl; }

      SpinBlock systemDot, environmentDot;
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
      systemDot = SpinBlock(systemDotStart, systemDotEnd, integralIndex, true);

      SpinBlock newSystem; // new system after blocking and decimating
      newSystem.initialise_op_array(OVERLAP, false);
      newSystem.setstoragetype(DISTRIBUTED_STORAGE);
      newSystem.BuildSumBlock (NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, system, systemDot);

      std::vector<Matrix> brarotateMatrix, ketrotateMatrix;
      LoadRotationMatrix(newSystem.get_sites(), brarotateMatrix, stateA);
      LoadRotationMatrix(newSystem.get_sites(), ketrotateMatrix, stateB);
      newSystem.transform_operators(brarotateMatrix, ketrotateMatrix);

      
      system = newSystem;
      p2out << system<<endl;
      
      SpinBlock::store (forward, system.get_sites(), system, stateA, stateB);	 	
      ++sweepParams.set_block_iter();
      
      sweepParams.savestate(forward, syssites.size());
      if (dmrginp.outputlevel() > 0)
         mcheck("at the end of sweep iteration");
    }

  pout << "\t\t\t ============================================================================ " << endl;

  // update the static number of iterations
  return ;
  
}

