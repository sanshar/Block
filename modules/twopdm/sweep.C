/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "sweeptwopdm.h"
#include "global.h"
#include "solver.h"
#include "initblocks.h"
#include "rotationmat.h"
#include "davidson.h"
#include "linear.h"
#include "guess_wavefunction.h"
#include "twopdm.h"
#include "density.h"
#include "davidson.h"
#include "pario.h"

#ifndef SERIAL
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif
using namespace boost;
using namespace std;


void calcenergy(array_4d<double>& twopdm, int state)
{
  using namespace SpinAdapted;
  Matrix onepdm(2*dmrginp.last_site(), 2*dmrginp.last_site()); onepdm = 0.0;
  int nelec = dmrginp.real_particle_number();

  for (int i=0; i<dmrginp.last_site()*2; i++)
  for (int j=0; j<dmrginp.last_site()*2; j++)
  for (int k=0; k<dmrginp.last_site()*2; k++)
    onepdm(i+1, j+1) += twopdm(i, k, k, j);

  onepdm /= (nelec-1);
  double nel = 0.0, sz=0.0;
  for (int i=0; i<dmrginp.last_site(); i++) {
    nel += onepdm(2*i+1, 2*i+1)+onepdm(2*i+2, 2*i+2);
    sz += onepdm(2*i+1, 2*i+1)-onepdm(2*i+2, 2*i+2);
  }

  double energy = 0.0;
  for (int i=0; i<dmrginp.last_site()*2; i++)
  for (int j=0; j<dmrginp.last_site()*2; j++)
  for (int k=0; k<dmrginp.last_site()*2; k++)
  for (int l=0; l<dmrginp.last_site()*2; l++)
    energy += v_2[0](i,j,k,l)*twopdm(i,j,l,k);

  energy *= 0.5;

  for (int i=0; i<dmrginp.last_site()*2; i++)
  for (int j=0; j<dmrginp.last_site()*2; j++)
    energy += v_1[0](i,j) * onepdm(i+1,j+1);

  pout << "energy of state "<< state <<" = "<< energy<<endl;

  ofstream out("onepdm_fromtpdm");
  for (int i=0; i<dmrginp.last_site()*2; i++)
  for (int j=0; j<dmrginp.last_site()*2; j++)
    out<<i<<"  "<<j<<"  "<<onepdm(i+1,j+1)<<endl;
  out.close();
  
}

namespace SpinAdapted{
void SweepTwopdm::BlockAndDecimate (SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, const bool &useSlater, const bool& dot_with_sys, int state)
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
  //if (useSlater) {
  systemDot = SpinBlock(systemDotStart, systemDotEnd, system.get_integralIndex(), true);
    //SpinBlock::store(true, systemDot.get_sites(), systemDot);
    //}
    //else
    //SpinBlock::restore(true, spindotsites, systemDot);
  SpinBlock environment, environmentDot, newEnvironment;

  int environmentDotStart, environmentDotEnd, environmentStart, environmentEnd;

  const int nexact = forward ? sweepParams.get_forward_starting_size() : sweepParams.get_backward_starting_size();

  system.addAdditionalCompOps();
  InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, sweepParams.current_root(), sweepParams.current_root(), sweepParams.get_sys_add(), dmrginp.direct(), system.get_integralIndex(), DISTRIBUTED_STORAGE, true, true);
  
  InitBlocks::InitNewEnvironmentBlock(environment, systemDot, newEnvironment, system, systemDot, sweepParams.current_root(), sweepParams.current_root(), 
				      sweepParams.get_sys_add(), sweepParams.get_env_add(), forward, dmrginp.direct(),
				      sweepParams.get_onedot(), nexact, useSlater, system.get_integralIndex(), true, true, true);
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
  DensityMatrix tracedMatrix(newSystem.get_stateInfo());
  tracedMatrix.allocate(newSystem.get_stateInfo());
  tracedMatrix.makedensitymatrix(solution, big, std::vector<double>(1,1.0), 0.0, 0.0, false);
  rotateMatrix.clear();
  if (!mpigetrank())
    double error = makeRotateMatrix(tracedMatrix, rotateMatrix, sweepParams.get_keep_states(), sweepParams.get_keep_qstates());
  

#ifndef SERIAL
  mpi::broadcast(world,rotateMatrix,0);
#endif
#ifdef SERIAL
  const int numprocs = 1;
#endif
#ifndef SERIAL
  const int numprocs = world.size();
#endif
  if (sweepParams.get_block_iter() == 0)
    compute_twopdm_initial(solution, system, systemDot, newSystem, newEnvironment, big, numprocs, state);

  compute_twopdm_sweep(solution, system, systemDot, newSystem, newEnvironment, big, numprocs, state);

  if (sweepParams.get_block_iter()  == sweepParams.get_n_iters() - 1)
    compute_twopdm_final(solution, system, systemDot, newSystem, newEnvironment, big, numprocs, state);

  SaveRotationMatrix (newSystem.get_sites(), rotateMatrix, state);

  //for(int i=0;i<dmrginp.nroots();++i)
  solution[0].SaveWavefunctionInfo (big.get_stateInfo(), big.get_leftBlock()->get_sites(), state);

  newSystem.transform_operators(rotateMatrix);

}

double SweepTwopdm::do_one(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize, int state)
{
  Timer sweeptimer;
  int integralIndex = 0;
  if (dmrginp.hamiltonian() == BCS) {
    pout << "Two PDM with BCS calculations is not implemented" << endl;
    exit(0);
  }
  pout.precision(12);
  SpinBlock system;
  const int nroots = dmrginp.nroots();
  std::vector<double> finalEnergy(nroots,0.);
  std::vector<double> finalEnergy_spins(nroots,0.);
  double finalError = 0.;

  sweepParams.set_sweep_parameters();
  // a new renormalisation sweep routine
  pout << ((forward) ? "\t\t\t Starting renormalisation sweep in forwards direction" : "\t\t\t Starting renormalisation sweep in backwards direction") << endl;
  pout << "\t\t\t ============================================================================ " << endl;
  
  InitBlocks::InitStartingBlock (system,forward, sweepParams.current_root(), sweepParams.current_root(), sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), restartSize, restart, warmUp, integralIndex);
  if(!restart)
    sweepParams.set_block_iter() = 0;
 
  pout << "\t\t\t Starting block is :: " << endl << system << endl;
  if (!restart) 
    SpinBlock::store (forward, system.get_sites(), system, sweepParams.current_root(), sweepParams.current_root()); // if restart, just restoring an existing block --
  sweepParams.savestate(forward, system.get_sites().size());
  bool dot_with_sys = true;

  array_4d<double> twopdm(2*dmrginp.last_site(), 2*dmrginp.last_site(), 2*dmrginp.last_site(), 2*dmrginp.last_site());
  twopdm.Clear();

  save_twopdm_binary(twopdm, state, state); 


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

      BlockAndDecimate (sweepParams, system, newSystem, warmUp, dot_with_sys, state);

      for(int j=0;j<nroots;++j)
        pout << "\t\t\t Total block energy for State [ " << j << 
	  " ] with " << sweepParams.get_keep_states()<<" :: " << sweepParams.get_lowest_energy()[j] <<endl;              

      finalEnergy_spins = ((sweepParams.get_lowest_energy()[0] < finalEnergy[0]) ? sweepParams.get_lowest_energy_spins() : finalEnergy_spins);
      finalEnergy = ((sweepParams.get_lowest_energy()[0] < finalEnergy[0]) ? sweepParams.get_lowest_energy() : finalEnergy);
      finalError = max(sweepParams.get_lowest_error(),finalError);

      system = newSystem;

      pout << system<<endl;
      
      SpinBlock::store (forward, system.get_sites(), system, sweepParams.current_root(), sweepParams.current_root());	 	

      p1out << "\t\t\t saving state " << system.get_sites().size() << endl;
      ++sweepParams.set_block_iter();
      //sweepParams.savestate(forward, system.get_sites().size());
    }
  //for(int j=0;j<nroots;++j)
  {int j = state;
    pout << "\t\t\t Finished Sweep with " << sweepParams.get_keep_states() << " states and sweep energy for State [ " << j 
	 << " ] with Spin [ " << dmrginp.molecule_quantum().get_s()  << " ] :: " << finalEnergy[j] << endl;
  }
  pout << "\t\t\t Largest Error for Sweep with " << sweepParams.get_keep_states() << " states is " << finalError << endl;
  pout << "\t\t\t ============================================================================ " << endl;

  int i = state, j = state;
  //for (int j=0; j<=i; j++) {
  load_twopdm_binary(twopdm, i, j); 
  //calcenergy(twopdm, i);
  save_twopdm_text(twopdm, i, j);
  save_spatial_twopdm_text(twopdm, i, j);
  save_spatial_twopdm_binary(twopdm, i, j);
  

  // update the static number of iterations

  ++sweepParams.set_sweep_iter();

  ecpu = sweeptimer.elapsedcputime(); ewall = sweeptimer.elapsedwalltime();
  pout << "\t\t\t Elapsed Sweep CPU  Time (seconds): " << setprecision(3) << ecpu << endl;
  pout << "\t\t\t Elapsed Sweep Wall Time (seconds): " << setprecision(3) << ewall << endl;


  return finalEnergy[0];
}


}
