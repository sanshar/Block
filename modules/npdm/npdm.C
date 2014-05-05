/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "npdm.h"
#include "sweep.h"
#include "sweepgenblock.h"
#include "density.h"
#include "sweeponepdm.h"  // For legacy version of 1pdm
#include "sweeptwopdm.h"  // For legacy version of 2pdm
#include "npdm_driver.h"
#include "nevpt2_npdm_driver.h"
#include "transition_sweeponepdm.h"
#include "transition_sweeptwopdm.h"

void dmrg(double sweep_tol);
void restart(double sweep_tol, bool reset_iter);
void ReadInput(char* conf);
void fullrestartGenblock();

namespace SpinAdapted {
namespace Npdm {

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void npdm_block_and_decimate( Npdm_driver_base& npdm_driver, SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, 
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
    systemDot = SpinBlock(systemDotStart, systemDotEnd, true);
    //SpinBlock::store(true, systemDot.get_sites(), systemDot);
    //}
    //else
    //SpinBlock::restore(true, spindotsites, systemDot);
  SpinBlock environment, environmentDot, newEnvironment;

  int environmentDotStart, environmentDotEnd, environmentStart, environmentEnd;

  const int nexact = forward ? sweepParams.get_forward_starting_size() : sweepParams.get_backward_starting_size();

  system.addAdditionalCompOps();
  InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, sweepParams.current_root(), sweepParams.current_root(),
                                 sweepParams.get_sys_add(), dmrginp.direct(), DISTRIBUTED_STORAGE, true, true);
  
  InitBlocks::InitNewEnvironmentBlock(environment, systemDot, newEnvironment, system, systemDot,
                                      sweepParams.current_root(), sweepParams.current_root(),
                                      sweepParams.get_sys_add(), sweepParams.get_env_add(), forward, dmrginp.direct(),
                                      sweepParams.get_onedot(), nexact, useSlater, true, true, true);
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

  int sweepPos = sweepParams.get_block_iter();
  int endPos = sweepParams.get_n_iters()-1;

  npdm_driver.compute_npdm_elements(solution, big, sweepPos, endPos);

  SaveRotationMatrix (newSystem.get_sites(), rotateMatrix, state);

  //FIXME multiple roots
  //for(int i=0;i<dmrginp.nroots();++i)
  solution[0].SaveWavefunctionInfo (big.get_stateInfo(), big.get_leftBlock()->get_sites(), state);

  //FIXME Do we need to do this step for NPDM on the last sweep site? (It's not negligible cost...?)
  newSystem.transform_operators(rotateMatrix);
  pout << "NPDM block and decimate and compute elements " << timer.elapsedwalltime() << " " << timer.elapsedcputime() << endl;

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

double npdm_do_one_sweep(Npdm_driver_base& npdm_driver, SweepParams &sweepParams, const bool &warmUp, const bool &forward, 
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
  
  InitBlocks::InitStartingBlock( system, forward, sweepParams.current_root(), sweepParams.current_root(), 
                                 sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), restartSize, restart, warmUp);
  pout << "\t\t\t Starting block is :: " << endl << system << endl;

  if (!restart) sweepParams.set_block_iter() = 0;
  if (!restart) SpinBlock::store (forward, system.get_sites(), system, sweepParams.current_root(), sweepParams.current_root() ); // if restart, just restoring an existing block --

  sweepParams.savestate(forward, system.get_sites().size());
  bool dot_with_sys = true;

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

    // Build npdm elements
    npdm_block_and_decimate(npdm_driver, sweepParams, system, newSystem, warmUp, dot_with_sys, state);

    for(int j=0;j<nroots;++j)
      pout << "\t\t\t Total block energy for State [ " << j << 
 " ] with " << sweepParams.get_keep_states()<<" :: " << sweepParams.get_lowest_energy()[j]+dmrginp.get_coreenergy() <<endl;              

    finalEnergy_spins = ((sweepParams.get_lowest_energy()[0] < finalEnergy[0]) ? sweepParams.get_lowest_energy_spins() : finalEnergy_spins);
    finalEnergy = ((sweepParams.get_lowest_energy()[0] < finalEnergy[0]) ? sweepParams.get_lowest_energy() : finalEnergy);
    finalError = max(sweepParams.get_lowest_error(),finalError);

    system = newSystem;

    pout << system<<endl;
    
    SpinBlock::store (forward, system.get_sites(), system, sweepParams.current_root(), sweepParams.current_root() );

    pout << "\t\t\t saving state " << system.get_sites().size() << endl;
    ++sweepParams.set_block_iter();
    //sweepParams.savestate(forward, system.get_sites().size());

    pout << "NPDM do one site time " << timer.elapsedwalltime() << " " << timer.elapsedcputime() << endl;
  }

  //for(int j=0;j<nroots;++j)
  {int j = state;
    pout << "\t\t\t Finished Sweep with " << sweepParams.get_keep_states() << " states and sweep energy for State [ " << j 
	 << " ] with Spin [ " << dmrginp.molecule_quantum().get_s()  << " ] :: " << finalEnergy[j]+dmrginp.get_coreenergy() << endl;
  }
  pout << "\t\t\t Largest Error for Sweep with " << sweepParams.get_keep_states() << " states is " << finalError << endl;
  pout << "\t\t\t ============================================================================ " << endl;

  // Dump NPDM to disk if necessary
  npdm_driver.save_data( state, state );

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

  if(sym == "dinfh") {
    pout << "Npdm not implemented with dinfh symmetry"<<endl;
    abort();
  }

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
  dmrginp.oneindex_screen_tol() = 0.0; //need to turn screening off for one index ops
  dmrginp.twoindex_screen_tol() = 0.0; //need to turn screening off for two index ops
  dmrginp.Sz() = dmrginp.total_spin_number().getirrep();
  dmrginp.do_npdm_ops() = true;
  sweep_copy.restorestate(direction_copy, restartsize_copy);

  // Initialize npdm_driver
  boost::shared_ptr<Npdm_driver_base> npdm_driver;
  if ( (dmrginp.hamiltonian() == QUANTUM_CHEMISTRY) && dmrginp.spinAdapted() ) {
    dmrginp.new_npdm_code() = true;
    if      (npdm_order == 1) npdm_driver = boost::shared_ptr<Npdm_driver_base>( new Onepdm_driver( dmrginp.last_site() ) );
    else if (npdm_order == 2) npdm_driver = boost::shared_ptr<Npdm_driver_base>( new Twopdm_driver( dmrginp.last_site() ) );
    else if (npdm_order == 3) npdm_driver = boost::shared_ptr<Npdm_driver_base>( new Threepdm_driver( dmrginp.last_site() ) );
    else if (npdm_order == 4) npdm_driver = boost::shared_ptr<Npdm_driver_base>( new Fourpdm_driver( dmrginp.last_site() ) );
    else if (npdm_order == 0) npdm_driver = boost::shared_ptr<Npdm_driver_base>( new Nevpt2_npdm_driver( dmrginp.last_site() ) );
    else abort();
  }

  // Not state-specific
  //--------------------
  if ( !dmrginp.setStateSpecific() ) {
    // Prepare NPDM operators
    Timer timer;
    dmrginp.set_fullrestart() = true;
    sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
    SweepGenblock::do_one(sweepParams, false, !direction, false, 0, -1, -1); //this will generate the cd operators                               
    pout << "NPDM SweepGenblock time " << timer.elapsedwalltime() << " " << timer.elapsedcputime() << endl;
    dmrginp.set_fullrestart() = false;
  
    for (int state=0; state<dmrginp.nroots(); state++) {
      sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
      // Do NPDM sweep
      if ( dmrginp.new_npdm_code() ) {
        Timer timerX;
        npdm_driver->clear();
        npdm_do_one_sweep(*npdm_driver, sweepParams, false, direction, false, 0, state);
        pout << "NPDM sweep time " << timerX.elapsedwalltime() << " " << timerX.elapsedcputime() << endl;
      } 
      else {
        if (npdm_order == 1) SweepOnepdm::do_one(sweepParams, false, direction, false, 0, state);      // Compute onepdm with the original code
        else if (npdm_order == 2) SweepTwopdm::do_one(sweepParams, false, direction, false, 0, state); // Compute twopdm with the original code
        else abort();
      }
    }
  }

  // State-specific
  //----------------
  else {
    // this only generates onpem between the same wavefunctions and cannot generate transition pdms. atleast not for now
    for (int state=0; state<dmrginp.nroots(); state++) {
      dmrginp.set_fullrestart() = true;
      sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;

      if (mpigetrank() == 0) {
        Sweep::InitializeStateInfo(sweepParams, direction, state);
        Sweep::InitializeStateInfo(sweepParams, !direction, state);
        Sweep::CanonicalizeWavefunction(sweepParams, direction, state);
        Sweep::CanonicalizeWavefunction(sweepParams, !direction, state);
        Sweep::CanonicalizeWavefunction(sweepParams, direction, state);
      }
      // Prepare NPDM operators
      SweepGenblock::do_one(sweepParams, false, !direction, false, 0, state, state); //this will generate the cd operators
      dmrginp.set_fullrestart() = false;
      // Do NPDM sweep
      if ( dmrginp.new_npdm_code() ) {
        npdm_driver->clear();
        npdm_do_one_sweep(*npdm_driver, sweepParams, false, direction, false, 0, state);
      }
      else {
        if (npdm_order == 1) SweepOnepdm::do_one(sweepParams, false, direction, false, 0, state);      // Compute onepdm with the original code
        else if (npdm_order == 2) SweepTwopdm::do_one(sweepParams, false, direction, false, 0, state); // Compute twopdm with the original code
        else abort();
      }
    }
  }

  sweep_copy.savestate(direction_copy, restartsize_copy);

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Note this routine calls the old onepdm and twopdm code
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

  if (dmrginp.algorithm_method() == TWODOT) {
    pout << "Npdm not allowed with twodot algorithm" << endl;
    abort();
  }

  dmrginp.oneindex_screen_tol() = 0.0; //need to turn screening off for one index ops
  dmrginp.twoindex_screen_tol() = 0.0; //need to turn screening off for two index ops
  dmrginp.Sz() = dmrginp.total_spin_number().getirrep();
  dmrginp.do_npdm_ops() = true;

  sweep_copy.restorestate(direction_copy, restartsize_copy);

  // Not state-specific
  //--------------------
  if ( !dmrginp.setStateSpecific() ) {
    dmrginp.set_fullrestart() = true;
    sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
    SweepGenblock::do_one(sweepParams, false, !direction, false, 0, -1, -1); //this will generate the cd operators                               
    dmrginp.set_fullrestart() = false;
  
    for (int state=0; state<dmrginp.nroots(); state++) {
      sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
      // 1PDM
      if (npdm_order == 1) SweepOnepdm::do_one(sweepParams, false, direction, false, 0, state);
      // 2PDM
      else if (npdm_order == 2) SweepTwopdm::do_one(sweepParams, false, direction, false, 0, state);
      else abort();
    }
  } 

  // State-specific
  //----------------
  else {
    // this only generates onpem between the same wavefunctions and cannot generate
    // transition pdms. atleast not for now
    for (int state=0; state<dmrginp.nroots(); state++) {
      dmrginp.set_fullrestart() = true;
      sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
      
      if (mpigetrank() == 0) {
        Sweep::InitializeStateInfo(sweepParams, direction, state);
        Sweep::InitializeStateInfo(sweepParams, !direction, state);
        Sweep::CanonicalizeWavefunction(sweepParams, direction, state);
        Sweep::CanonicalizeWavefunction(sweepParams, !direction, state);
        Sweep::CanonicalizeWavefunction(sweepParams, direction, state);
      }
      SweepGenblock::do_one(sweepParams, false, !direction, false, 0, state, state); //this will generate the cd operators
      dmrginp.set_fullrestart() = false;
      // 1PDM
      if (npdm_order == 1) SweepOnepdm::do_one(sweepParams, false, direction, false, 0, state);
      // 2PDM
      else if (npdm_order == 2) SweepTwopdm::do_one(sweepParams, false, direction, false, 0, state);
      else abort();
    }
  }

  sweep_copy.savestate(direction_copy, restartsize_copy);
}


//===========================================================================================================================================================

void transition_pdm( int npdm_order )
{
  double sweep_tol = 1e-7;
  sweep_tol = dmrginp.get_sweep_tol();
  bool direction;
  int restartsize;
  bool direction_copy; int restartsize_copy;
  SweepParams sweepParams;
  SweepParams sweep_copy;
  dmrginp.setimplicitTranspose() = false;

  if(sym == "dinfh") {
    pout << "Npdm not implemented with dinfh symmetry"<<endl;
    abort();
  }

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
  dmrginp.oneindex_screen_tol() = 0.0; //need to turn screening off for one index ops
  dmrginp.twoindex_screen_tol() = 0.0; //need to turn screening off for two index ops
  dmrginp.Sz() = dmrginp.total_spin_number().getirrep();
  dmrginp.do_npdm_ops() = true;
  sweep_copy.restorestate(direction_copy, restartsize_copy);

  // Initialize npdm_driver
  boost::shared_ptr<Npdm_driver_base> npdm_driver;
  if ( (dmrginp.hamiltonian() == QUANTUM_CHEMISTRY) && dmrginp.spinAdapted() ) {
    dmrginp.new_npdm_code() = false;
    if      (npdm_order == 1) npdm_driver = boost::shared_ptr<Npdm_driver_base>( new Onepdm_driver( dmrginp.last_site() ) );
    else if (npdm_order == 2) npdm_driver = boost::shared_ptr<Npdm_driver_base>( new Twopdm_driver( dmrginp.last_site() ) );
    else if (npdm_order == 3) npdm_driver = boost::shared_ptr<Npdm_driver_base>( new Threepdm_driver( dmrginp.last_site() ) );
    else if (npdm_order == 4) npdm_driver = boost::shared_ptr<Npdm_driver_base>( new Fourpdm_driver( dmrginp.last_site() ) );
    else if (npdm_order == 0) npdm_driver = boost::shared_ptr<Npdm_driver_base>( new Nevpt2_npdm_driver( dmrginp.last_site() ) );
    else abort();
  }

  //--------------------
    // Prepare NPDM operators
  Timer timer;
  dmrginp.set_fullrestart() = true;
  sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
  SweepGenblock::do_one(sweepParams, false, !direction, false, 0, -1, -1); //this will generate the cd operators                               
  pout << "NPDM SweepGenblock time " << timer.elapsedwalltime() << " " << timer.elapsedcputime() << endl;
  dmrginp.set_fullrestart() = false;
  
  for (int state=0; state<dmrginp.nroots(); state++) {
     for(int stateB=0; stateB<= state; stateB++){
   //  Usually, interacting operator is hermitian, <\Phi_1|V|\Phi_2> = <\Phi_2|V|\Phi_1>
   //  for(int stateB=0; stateB<= dmrginp.nroots(); stateB++){
        sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
        if (npdm_order == 1) transition_onepdm::do_one(sweepParams, false, direction, false, 0, state,stateB);     
        else if (npdm_order == 2) transition_twopdm::do_one(sweepParams, false, direction, false, 0, state,stateB);
        else abort();
  }
  }

  sweep_copy.savestate(direction_copy, restartsize_copy);

}

void transition_pdm_restart( int npdm_order )
{
  bool direction;
  int restartsize;
  bool direction_copy; int restartsize_copy;
  SweepParams sweepParams;
  SweepParams sweep_copy;
  dmrginp.setimplicitTranspose() = false;

  if(sym == "dinfh") {
    pout << "Npdm not implemented with dinfh symmetry"<<endl;
    abort();
  }

  if (dmrginp.algorithm_method() == TWODOT) {
    pout << "Npdm not allowed with twodot algorithm" << endl;
    abort();
  }

  dmrginp.oneindex_screen_tol() = 0.0; //need to turn screening off for one index ops
  dmrginp.twoindex_screen_tol() = 0.0; //need to turn screening off for two index ops
  dmrginp.Sz() = dmrginp.total_spin_number().getirrep();
  dmrginp.do_npdm_ops() = true;

  sweep_copy.restorestate(direction_copy, restartsize_copy);

  // Not state-specific
  //--------------------
  if ( !dmrginp.setStateSpecific() ) {
    dmrginp.set_fullrestart() = true;
    sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
    SweepGenblock::do_one(sweepParams, false, !direction, false, 0, -1, -1); //this will generate the cd operators                               
    dmrginp.set_fullrestart() = false;
  
    for (int state=0; state<dmrginp.nroots(); state++) {
     for(int stateB=0; stateB<= state; stateB++){
   //  Usually, interacting operator is hermitian, <\Phi_1|V|\Phi_2> = <\Phi_2|V|\Phi_1>
  //    for(int stateB=0; stateB < dmrginp.nroots(); stateB++){
        sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
        // 1PDM
        if (npdm_order == 1) transition_onepdm::do_one(sweepParams, false, direction, false, 0, state,stateB);
        // 2PDM
        else if (npdm_order == 2) transition_twopdm::do_one(sweepParams, false, direction, false, 0, state,stateB);
        else abort();
      }
    }
  } 

  sweep_copy.savestate(direction_copy, restartsize_copy);
}
}
}

