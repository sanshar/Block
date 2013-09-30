/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "npdm.h"
#include "sweeptwopdm.h"  // For old version of 2pdm

void dmrg(double sweep_tol);
void restart(double sweep_tol, bool reset_iter);
void ReadInput(char* conf);
void fullrestartGenblock();

namespace SpinAdapted {
namespace Npdm {

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void npdm_block_and_decimate( Npdm_driver& npdm_driver, SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, 
                              const bool &useSlater, const bool& dot_with_sys, const int state)
{
  Timer timer;
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
  InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, sweepParams.get_sys_add(), dmrginp.direct(), DISTRIBUTED_STORAGE, true, true);
  
  InitBlocks::InitNewEnvironmentBlock(environment, systemDot, newEnvironment, system, systemDot,
                                      sweepParams.get_sys_add(), sweepParams.get_env_add(), forward, dmrginp.direct(),
                                      sweepParams.get_onedot(), nexact, useSlater, true, true, true);
  SpinBlock big;
  newSystem.set_loopblock(true);
  system.set_loopblock(false);
  newEnvironment.set_loopblock(false);
//cout << "newSystem: LHS CRE size on p" << mpigetrank() << " = " << newSystem.get_leftBlock()->get_op_array(CRE).get_size() << " local ops" << endl;
//cout << "newSystem: RHS CRE size on p" << mpigetrank() << " = " << newSystem.get_rightBlock()->get_op_array(CRE).get_size() << " local ops" << endl;
//cout << "newEnvironment: CRE size on p" << mpigetrank() << " = " << newEnvironment.get_op_array(CRE).get_size() << " local ops" << endl;
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

//MAW
  int sweepPos = sweepParams.get_block_iter();
  int endPos = sweepParams.get_n_iters()-1;
//cout << "big: LHS CRE size on p" << mpigetrank() << " = " << big.get_leftBlock()->get_leftBlock()->get_op_array(CRE).get_size() << " local ops" << endl;
//cout << "big: DOT CRE size on p" << mpigetrank() << " = " << big.get_leftBlock()->get_rightBlock()->get_op_array(CRE).get_size() << " local ops" << endl;
//cout << "big: RHS CRE size on p" << mpigetrank() << " = " << big.get_rightBlock()->get_op_array(CRE).get_size() << " local ops" << endl;
  npdm_driver.compute_npdm_elements(solution, big, sweepPos, endPos);
//MAW

  SaveRotationMatrix (newSystem.get_sites(), rotateMatrix, state);

  //for(int i=0;i<dmrginp.nroots();++i)
  solution[0].SaveWavefunctionInfo (big.get_stateInfo(), big.get_leftBlock()->get_sites(), state);

  //FIXME Do we need to do this step for NPDM on the last sweep site? (It's not negligible cost...?)
  newSystem.transform_operators(rotateMatrix);
  pout << "NPDM block and decimate and compute elements " << timer.elapsedwalltime() << " " << timer.elapsedcputime() << endl;

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

double npdm_do_one_sweep(Npdm_driver& npdm_driver, SweepParams &sweepParams, const bool &warmUp, const bool &forward, 
                         const bool &restart, const int &restartSize, const int state)
{
  cout.precision(12);
  SpinBlock system;
  const int nroots = dmrginp.nroots();
  std::vector<double> finalEnergy(nroots,0.);
  std::vector<double> finalEnergy_spins(nroots,0.);
  double finalError = 0.;

  sweepParams.set_sweep_parameters();
  // a new renormalisation sweep routine
  pout << ((forward) ? "\t\t\t Starting renormalisation sweep in forwards direction" : "\t\t\t Starting renormalisation sweep in backwards direction") << endl;
  pout << "\t\t\t ============================================================================ " << endl;
  
  InitBlocks::InitStartingBlock( system, forward, sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), restartSize, restart, warmUp);
  pout << "\t\t\t Starting block is :: " << endl << system << endl;

  if (!restart) sweepParams.set_block_iter() = 0;
  if (!restart) SpinBlock::store (forward, system.get_sites(), system); // if restart, just restoring an existing block --

  sweepParams.savestate(forward, system.get_sites().size());
  bool dot_with_sys = true;

  // NPDM storage is either as full array_Nd<T>, which we fully allocate here, or in sparse format allocated dynamically
  if ( npdm_driver.use_full_array_ ) {
    npdm_driver.resize_npdm_array(2*dmrginp.last_site());
    npdm_driver.clear_npdm_array();
  }

  for (int i=0; i<nroots; i++) {
//FIXME only allows for one root at present (actually, should be trivial to change...)
assert(i==0);
//MAW    save_npdm_binary(i, i); 

    // Loop over all block sites
    for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); ) {
      Timer timer;

      pout << "\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
      pout << "\t\t\t ----------------------------" << endl;
      if (forward) pout << "\t\t\t Current direction is :: Forwards " << endl;
      else pout << "\t\t\t Current direction is :: Backwards " << endl;

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

      // Build Npdm elements
      npdm_block_and_decimate(npdm_driver, sweepParams, system, newSystem, warmUp, dot_with_sys, state);

      for(int j=0;j<nroots;++j)
        pout << "\t\t\t Total block energy for State [ " << j << 
	  " ] with " << sweepParams.get_keep_states()<<" :: " << sweepParams.get_lowest_energy()[j]+dmrginp.get_coreenergy() <<endl;              

      finalEnergy_spins = ((sweepParams.get_lowest_energy()[0] < finalEnergy[0]) ? sweepParams.get_lowest_energy_spins() : finalEnergy_spins);
      finalEnergy = ((sweepParams.get_lowest_energy()[0] < finalEnergy[0]) ? sweepParams.get_lowest_energy() : finalEnergy);
      finalError = max(sweepParams.get_lowest_error(),finalError);

      system = newSystem;

      pout << system<<endl;
      
      SpinBlock::store (forward, system.get_sites(), system);	 	

      pout << "\t\t\t saving state " << system.get_sites().size() << endl;
      ++sweepParams.set_block_iter();
      sweepParams.savestate(forward, system.get_sites().size());

      pout << "NPDM do one site time " << timer.elapsedwalltime() << " " << timer.elapsedcputime() << endl;
    }
  }

  //for(int j=0;j<nroots;++j)
  {int j = state;
    pout << "\t\t\t Finished Sweep with " << sweepParams.get_keep_states() << " states and sweep energy for State [ " << j 
	 << " ] with Spin [ " << dmrginp.molecule_quantum().get_s()  << " ] :: " << finalEnergy[j]+dmrginp.get_coreenergy() << endl;
  }
  pout << "\t\t\t Largest Error for Sweep with " << sweepParams.get_keep_states() << " states is " << finalError << endl;
  pout << "\t\t\t ============================================================================ " << endl;


  // Combine NPDM elements from all mpi ranks and dump files
  assert(state==0);
  int i = state, j = state;
  npdm_driver.save_sparse_array(i,j);
  if ( npdm_driver.use_full_array_ ) npdm_driver.save_full_array(i,j);

  // Update the static number of iterations
  ++sweepParams.set_sweep_iter();

  return finalEnergy[0];

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void npdm( int npdm_order )
{

  double sweep_tol = 1e-7;
  sweep_tol = dmrginp.get_sweep_tol();
  bool direction;
  int restartsize;
  bool direction_copy; int restartsize_copy;
  SweepParams sweepParams;
  SweepParams sweep_copy;

  if (dmrginp.algorithm_method() == TWODOT) {
    pout << "Npdm not allowed with twodot algorithm" << endl;
    abort();
  }

  if (RESTART && !FULLRESTART)
    restart(sweep_tol, reset_iter);
  else if (FULLRESTART) {
    fullrestartGenblock();
    reset_iter = true;
    sweepParams.restorestate(direction, restartsize);
    sweepParams.calc_niter();
    sweepParams.savestate(direction, restartsize);
    restart(sweep_tol, reset_iter);
  }
  else {
    dmrg(sweep_tol);
  }

  // Screening can break things for NPDM (e.g. smaller operators won't be available from which to build larger ones etc...?)
  dmrginp.screen_tol() = 0.0;
  dmrginp.Sz() = dmrginp.total_spin_number();
  dmrginp.do_npdm_ops() = true;

  // Prepare NPDM operators
  Timer timer;
  sweep_copy.restorestate(direction_copy, restartsize_copy);
  dmrginp.set_fullrestart() = true;
  sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
  SweepGenblock::do_one(sweepParams, false, !direction, false, 0, 0);
  dmrginp.set_fullrestart() = false;
  pout << "NPDM SweepGenblock time " << timer.elapsedwalltime() << " " << timer.elapsedcputime() << endl;

  switch (npdm_order) {
  case (1):
    // Compute onepdm elements
    SweepOnepdm::do_one(sweepParams, false, direction, false, 0);
    sweep_copy.savestate(direction_copy, restartsize_copy);
    break;
  case (2):
    // Compute twopdm elements
    for (int state=0; state<dmrginp.nroots(); state++) {
      Twopdm_driver twopdm_driver;
      sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
      if (false) {
        // Compute twopdm with the original code
        SweepTwopdm::do_one(sweepParams, false, direction, false, 0, state);
      } 
      else {
        // Compute twopdm with general npdm code
        npdm_do_one_sweep(twopdm_driver, sweepParams, false, direction, false, 0, state);
      }
    }
    break;
  case (3):
    // Compute threepdm elements
    for (int state=0; state<dmrginp.nroots(); state++) {
      Threepdm_driver threepdm_driver;
      sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
      npdm_do_one_sweep(threepdm_driver, sweepParams, false, direction, false, 0, state);
    }
    break;
  case (4):
    // Compute fourpdm elements
    for (int state=0; state<dmrginp.nroots(); state++) {
      Fourpdm_driver fourpdm_driver;
      sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
      // Not all 4-index ops are implemented yet!!
      assert(false);
      npdm_do_one_sweep(fourpdm_driver, sweepParams, false, direction, false, 0, state);
    }
    break;
  }
  sweep_copy.savestate(direction_copy, restartsize_copy);

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void npdm_restart( int npdm_order )
{
  bool direction;
  int restartsize;
  bool direction_copy; int restartsize_copy;
  SweepParams sweepParams;
  SweepParams sweep_copy;

  if(sym == "dinfh") {
    pout << "Npdm not implemented with dinfh symmetry"<<endl;
    abort();
  }

  sweepParams.restorestate(direction, restartsize);
  if(!sweepParams.get_onedot() || dmrginp.algorithm_method() == TWODOT) {
    pout << "Npdm only runs for the onedot algorithm" << endl;
    abort();
  }

  dmrginp.screen_tol() = 0.0; //need to turn screening off for onepdm
  dmrginp.Sz() = dmrginp.total_spin_number();
  dmrginp.do_npdm_ops() = true;
  dmrginp.screen_tol() = 0.0;

  sweep_copy.restorestate(direction_copy, restartsize_copy);
  dmrginp.set_fullrestart() = true;
  sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
  SweepGenblock::do_one(sweepParams, false, !direction, false, 0, 0); //this will generate the cd operators
  dmrginp.set_fullrestart() = false;

  switch (npdm_order) {
  case (1):
    SweepOnepdm::do_one(sweepParams, false, direction, false, 0);
    sweep_copy.savestate(direction_copy, restartsize_copy);
    break;
  case (2):
    for (int state=0; state<dmrginp.nroots(); state++) {
      Twopdm_driver twopdm_driver;
      sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
      npdm_do_one_sweep(twopdm_driver, sweepParams, false, direction, false, 0, state);
    }
    break;
  }
  sweep_copy.savestate(direction_copy, restartsize_copy);

}

//===========================================================================================================================================================

}
}

