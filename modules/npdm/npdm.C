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
#include "pario.h"
#include "ds1_sweeponepdm.h"  //EL
#include "ds0_sweeponepdm.h"  //EL


void dmrg(double sweep_tol);
void restart(double sweep_tol, bool reset_iter);
void ReadInput(char* conf);
void fullrestartGenblock();

namespace SpinAdapted {
namespace Npdm {

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void npdm_block_and_decimate( Npdm_driver_base& npdm_driver, SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, 
                              const bool &useSlater, const bool& dot_with_sys, const int state, const int stateB)
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
  if(dmrginp.setStateSpecific() || dmrginp.transition_diff_irrep()){
    InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, state, stateB,
                                   sweepParams.get_sys_add(), dmrginp.direct(), system.get_integralIndex(), DISTRIBUTED_STORAGE, true, true);
    
    InitBlocks::InitNewEnvironmentBlock(environment, systemDot, newEnvironment, system, systemDot,
                                      state, stateB,
                                      sweepParams.get_sys_add(), sweepParams.get_env_add(), forward, dmrginp.direct(),
                                      sweepParams.get_onedot(), nexact, useSlater, environment.get_integralIndex(), true, true, true);
  }
  else{
    InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, sweepParams.current_root(), sweepParams.current_root(),
                                   sweepParams.get_sys_add(), dmrginp.direct(), system.get_integralIndex(), DISTRIBUTED_STORAGE, true, true);
    
    InitBlocks::InitNewEnvironmentBlock(environment, systemDot, newEnvironment, system, systemDot,
                                      sweepParams.current_root(), sweepParams.current_root(),
                                      sweepParams.get_sys_add(), sweepParams.get_env_add(), forward, dmrginp.direct(),
                                      sweepParams.get_onedot(), nexact, useSlater, environment.get_integralIndex(), true, true, true);

  }
  SpinBlock big;
  newSystem.set_loopblock(true);
  system.set_loopblock(false);
  newEnvironment.set_loopblock(false);
  InitBlocks::InitBigBlock(newSystem, newEnvironment, big); 

  const int nroots = dmrginp.nroots();
  std::vector<Wavefunction> solution;
  if(state==stateB){
    solution.resize(1);
    DiagonalMatrix e;
    GuessWave::guess_wavefunctions(solution[0], e, big, sweepParams.get_guesstype(), true, state, true, 0.0); 

  }
  else{
    solution.resize(2);
    DiagonalMatrix e;
    GuessWave::guess_wavefunctions(solution[0], e, big, sweepParams.get_guesstype(), true, state, true, 0.0,false); 
    GuessWave::guess_wavefunctions(solution[1], e, big, sweepParams.get_guesstype(), true, stateB, true, 0.0,true); 
  }

#ifndef SERIAL
  mpi::communicator world;
  mpi::broadcast(world, solution, 0);
#endif


  //GuessWave::guess_wavefunctions(solution[0], e, big, sweepParams.get_guesstype(), true, state, false, 0.0); 
  //GuessWave::guess_wavefunctions(solution[1], e, big, sweepParams.get_guesstype(), true, stateB, true, 0.0); 



  //bra and ket rotation matrices are calculated from different density matrices.


  std::vector<Matrix> rotateMatrix;
  std::vector<Matrix> rotateMatrixB;

  if(state!=stateB){

    DensityMatrix tracedMatrix(newSystem.get_braStateInfo());
    tracedMatrix.allocate(newSystem.get_braStateInfo());
    tracedMatrix.makedensitymatrix(std::vector<Wavefunction>(1,solution[0]), big, std::vector<double>(1,1.0), 0.0, 0.0, false);
    rotateMatrix.clear();

    DensityMatrix tracedMatrixB(newSystem.get_ketStateInfo());
    tracedMatrixB.allocate(newSystem.get_ketStateInfo());
    tracedMatrixB.makedensitymatrix(std::vector<Wavefunction>(1,solution[1]), big, std::vector<double>(1,1.0), 0.0, 0.0, false);
    rotateMatrixB.clear();
    if (!mpigetrank()){
      double error = makeRotateMatrix(tracedMatrix, rotateMatrix, sweepParams.get_keep_states(), sweepParams.get_keep_qstates());
      error = makeRotateMatrix(tracedMatrixB, rotateMatrixB, sweepParams.get_keep_states(), sweepParams.get_keep_qstates());
    }
  }
  else{
    DensityMatrix tracedMatrix(newSystem.get_stateInfo());
    tracedMatrix.allocate(newSystem.get_stateInfo());
    tracedMatrix.makedensitymatrix(std::vector<Wavefunction>(1,solution[0]), big, std::vector<double>(1,1.0), 0.0, 0.0, false);
    rotateMatrix.clear();
    if (!mpigetrank()){
      double error = makeRotateMatrix(tracedMatrix, rotateMatrix, sweepParams.get_keep_states(), sweepParams.get_keep_qstates());
    }

  }


  


  int sweepPos = sweepParams.get_block_iter();
  int endPos = sweepParams.get_n_iters()-1;
  npdm_driver.compute_npdm_elements(solution, big, sweepPos, endPos);
  SaveRotationMatrix (newSystem.get_sites(), rotateMatrix, state);
  solution[0].SaveWavefunctionInfo (big.get_braStateInfo(), big.get_leftBlock()->get_sites(), state);

  if(state!=stateB){
    SaveRotationMatrix (newSystem.get_sites(), rotateMatrixB, stateB);
    solution[1].SaveWavefunctionInfo (big.get_ketStateInfo(), big.get_leftBlock()->get_sites(), stateB);
  }


  //FIXME
  //Maybe, for StateSpecific calculations, we can load rotation matrices, wavefuntions from the disk. 
  //There is no longer need to transform wavefuntions and to make rotation matrices from the density matrices.
  
  //FIXME
  //If in the state-average pdm, different states do not share the same rotation matrices as they do in energy calculations. Making rotation matrices from 
  //density matrices of different states is neccessary. 
  
  //if(newSystem.get_sites().size()>1)
  //if (!mpigetrank()){
  //LoadRotationMatrix (newSystem.get_sites(), rotateMatrix, state);
  //LoadRotationMatrix (newSystem.get_sites(), rotateMatrixB, stateB);
  //}
  #ifndef SERIAL
    mpi::broadcast(world,rotateMatrix,0);
    if(state!=stateB)
      mpi::broadcast(world,rotateMatrixB,0);
  #endif

  // Do we need to do this step for NPDM on the last sweep site? (It's not negligible cost...?)
  //It crashes at the last sweep site. 
  //Since it is useless, Just omit it at the last sweep site.
 // if( sweepParams.get_block_iter()  != sweepParams.get_n_iters() - 1)
  {
    if(state!=stateB)
      newSystem.transform_operators(rotateMatrix,rotateMatrixB);
    else
      newSystem.transform_operators(rotateMatrix);
  }

  //newSystem.transform_operators(rotateMatrix,rotateMatrixB);
  ecpu = timer.elapsedcputime();ewall=timer.elapsedwalltime();
  p3out << "NPDM block and decimate and compute elements " << ewall << " " << ecpu << endl;

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

double npdm_do_one_sweep(Npdm_driver_base &npdm_driver, SweepParams &sweepParams, const bool &warmUp, const bool &forward, 
                         const bool &restart, const int &restartSize, const int state, const int stateB)
{
  Timer sweeptimer;
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
  
  int integralIndex = 0;
  if(dmrginp.setStateSpecific() || dmrginp.transition_diff_irrep())
    InitBlocks::InitStartingBlock( system, forward, state, stateB, sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), restartSize, restart, warmUp, integralIndex);
  else 
    InitBlocks::InitStartingBlock( system, forward, sweepParams.current_root(), sweepParams.current_root(), sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), restartSize, restart, warmUp, integralIndex);

  pout << "\t\t\t Starting block is :: " << endl << system << endl;

  if (!restart) sweepParams.set_block_iter() = 0;
  if(dmrginp.setStateSpecific() || dmrginp.transition_diff_irrep()){
    if (!restart) SpinBlock::store (forward, system.get_sites(), system, state, stateB ); // if restart, just restoring an existing block --
  }
  else{
    if (!restart) SpinBlock::store (forward, system.get_sites(), system, sweepParams.current_root(), sweepParams.current_root() ); // if restart, just restoring an existing block --
  }

  sweepParams.savestate(forward, system.get_sites().size());
  bool dot_with_sys = true;

  // Loop over all block sites
  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); ) {
    Timer timer;

    pout << "\n\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
    pout << "\t\t\t ----------------------------" << endl;
    if (forward) { p1out << "\t\t\t Current direction is :: Forwards " << endl; }
    else { p1out << "\t\t\t Current direction is :: Backwards " << endl; }

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
    
    //pout << "guess wave funtion type: " << sweepParams.get_guesstype()<<endl;
    p1out << "\t\t\t Blocking and Decimating " << endl;
 
    SpinBlock newSystem;

    // Build npdm elements
    npdm_block_and_decimate(npdm_driver, sweepParams, system, newSystem, warmUp, dot_with_sys, state, stateB);

//    for(int j=0;j<nroots;++j)
//      pout << "\t\t\t Total block energy for State [ " << j << 
// " ] with " << sweepParams.get_keep_states()<<" :: " << sweepParams.get_lowest_energy()[j]+dmrginp.get_coreenergy() <<endl;              
//
//    finalEnergy_spins = ((sweepParams.get_lowest_energy()[0] < finalEnergy[0]) ? sweepParams.get_lowest_energy_spins() : finalEnergy_spins);
//    finalEnergy = ((sweepParams.get_lowest_energy()[0] < finalEnergy[0]) ? sweepParams.get_lowest_energy() : finalEnergy);
//    finalError = max(sweepParams.get_lowest_error(),finalError);

    system = newSystem;

    pout << system<<endl;
    
    if(dmrginp.setStateSpecific() || dmrginp.transition_diff_irrep())
      SpinBlock::store (forward, system.get_sites(), system, state, stateB);
    else
      SpinBlock::store (forward, system.get_sites(), system, sweepParams.current_root(), sweepParams.current_root() );

    p1out << "\t\t\t saving state " << system.get_sites().size() << endl;
    ++sweepParams.set_block_iter();
    //sweepParams.savestate(forward, system.get_sites().size());

    ecpu = timer.elapsedcputime();ewall=timer.elapsedwalltime();
    p3out << "NPDM do one site time " << ewall << " " << ecpu << endl;
  }

  //for(int j=0;j<nroots;++j)
//  {int j = state;
//    pout << "\t\t\t Finished Sweep with " << sweepParams.get_keep_states() << " states and sweep energy for State [ " << j 
//	 << " ] with Spin [ " << dmrginp.molecule_quantum().get_s()  << " ] :: " << finalEnergy[j]+dmrginp.get_coreenergy() << endl;
//  }
//  {int j = stateB;
//    pout << "\t\t\t Finished Sweep with " << sweepParams.get_keep_states() << " states and sweep energy for State [ " << j 
//	 << " ] with Spin [ " << dmrginp.molecule_quantum().get_s()  << " ] :: " << finalEnergy[j]+dmrginp.get_coreenergy() << endl;
//  }
//  pout << "\t\t\t Largest Error for Sweep with " << sweepParams.get_keep_states() << " states is " << finalError << endl;
//  pout << "\t\t\t ============================================================================ " << endl;

  // Dump NPDM to disk if necessary
  npdm_driver.save_data( state, stateB );

  // Update the static number of iterations
  ++sweepParams.set_sweep_iter();

  ecpu = sweeptimer.elapsedcputime();ewall=sweeptimer.elapsedwalltime();
  pout << "\t\t\t Elapsed Sweep CPU  Time (seconds): " << setprecision(3) << ecpu << endl;
  pout << "\t\t\t Elapsed Sweep Wall Time (seconds): " << setprecision(3) << ewall << endl;

  return finalEnergy[0];

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void npdm(NpdmOrder npdm_order, bool restartpdm, bool transitionpdm, bool dS)
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

  if(!restartpdm){
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
  }
  if(transitionpdm)
    dmrginp.setimplicitTranspose() = false;

  dmrginp.do_pdm() = true;

  // Screening can break things for NPDM (e.g. smaller operators won't be available from which to build larger ones etc...?)
  dmrginp.oneindex_screen_tol() = 0.0; //need to turn screening off for one index ops
  dmrginp.twoindex_screen_tol() = 0.0; //need to turn screening off for two index ops
  dmrginp.Sz() = dmrginp.total_spin_number().getirrep();
  dmrginp.do_npdm_ops() = true;
  sweep_copy.restorestate(direction_copy, restartsize_copy);

  // Initialize npdm_driver
  boost::shared_ptr<Npdm_driver_base> npdm_driver;
  //if ( (dmrginp.hamiltonian() == QUANTUM_CHEMISTRY) && dmrginp.spinAdapted() ) {
  if ( (dmrginp.hamiltonian() == QUANTUM_CHEMISTRY) || (dmrginp.hamiltonian() == BCS)) {
    //By default, new_npdm_code is false.
    //For npdm_order 1 or 2. new_npdm_code is determined by default or manual setting.
    //For the other situation, only old or new code is suitable.
    if(npdm_order == NPDM_PAIRMATRIX || npdm_order == NPDM_THREEPDM || npdm_order == NPDM_FOURPDM || npdm_order == NPDM_NEVPT2 ||  transitionpdm == true  || dmrginp.spinAdapted() == false || dmrginp.setStateSpecific())
      dmrginp.set_new_npdm_code();

    if(dmrginp.new_npdm_code()){
    if      (npdm_order == NPDM_ONEPDM) npdm_driver = boost::shared_ptr<Npdm_driver_base>( new Onepdm_driver( dmrginp.last_site() ) );
    else if ((npdm_order == NPDM_DS0)||(npdm_order == NPDM_DS1)) npdm_driver = boost::shared_ptr<Npdm_driver_base>( new Onepdm_driver( dmrginp.last_site() ) );
    else if (npdm_order == NPDM_TWOPDM) npdm_driver = boost::shared_ptr<Npdm_driver_base>( new Twopdm_driver( dmrginp.last_site() ) );
    else if (npdm_order == NPDM_THREEPDM) npdm_driver = boost::shared_ptr<Npdm_driver_base>( new Threepdm_driver( dmrginp.last_site() ) );
    else if (npdm_order == NPDM_FOURPDM) npdm_driver = boost::shared_ptr<Npdm_driver_base>( new Fourpdm_driver( dmrginp.last_site() ) );
    else if (npdm_order == NPDM_NEVPT2) npdm_driver = boost::shared_ptr<Npdm_driver_base>( new Nevpt2_npdm_driver( dmrginp.last_site() ) );
    else if (npdm_order == NPDM_PAIRMATRIX) npdm_driver = boost::shared_ptr<Npdm_driver_base>( new Pairpdm_driver( dmrginp.last_site() ) );
    else abort();
    }
  }

 
if (dS) {
    for (int state=0; state<dmrginp.nroots(); state++) {
      for(int stateB=0; stateB<state; stateB++){
       dmrginp.set_fullrestart() = true;
       sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
       SweepGenblock::do_one(sweepParams, false, direction, false, 0, state, stateB); //this will generate the cd operators
       dmrginp.set_fullrestart() = false;
       npdm_driver->clear();
          sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
          if (npdm_order == NPDM_DS1) ds1_onepdm::do_one(sweepParams, false, !direction, false, 0, state, stateB);
          else if (npdm_order == NPDM_DS0) ds0_onepdm::do_one(sweepParams, false, !direction, false, 0, state, stateB);
          else abort();
      }
  }
}
else {
  if (dmrginp.specificpdm().size()!=0)
  {
    Timer timer;
    dmrginp.set_fullrestart() = true;
    sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
	  dmrginp.npdm_generate() = true;
    if ( !dmrginp.setStateSpecific())
      SweepGenblock::do_one(sweepParams, false, !direction, false, 0, -1, -1); //this will generate the cd operators                               
    else if (dmrginp.specificpdm().size()==1)
      SweepGenblock::do_one(sweepParams, false, !direction, false, 0, dmrginp.specificpdm()[0], dmrginp.specificpdm()[0]); //this will generate the cd operators                               
    else if (dmrginp.specificpdm().size()==2)
      SweepGenblock::do_one(sweepParams, false, !direction, false, 0, dmrginp.specificpdm()[0], dmrginp.specificpdm()[1]); //this will generate the cd operators                               
    else abort();
		dmrginp.npdm_generate() = false;
    ecpu = timer.elapsedcputime();ewall=timer.elapsedwalltime();
    p3out << "\t\t\t NPDM SweepGenblock time " << ewall << " " << ecpu << endl;
    dmrginp.set_fullrestart() = false;

    sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
    Timer timerX;
    npdm_driver->clear();
    if (dmrginp.specificpdm().size()==1)
      npdm_do_one_sweep(*npdm_driver, sweepParams, false, direction, false, 0,dmrginp.specificpdm()[0],dmrginp.specificpdm()[0]);
    else if (dmrginp.specificpdm().size()==2)
      npdm_do_one_sweep(*npdm_driver, sweepParams, false, direction, false, 0,dmrginp.specificpdm()[0],dmrginp.specificpdm()[1]);
    else abort();
    ecpu = timerX.elapsedcputime();ewall=timerX.elapsedwalltime();
    p3out << "\t\t\t NPDM sweep time " << ewall << " " << ecpu << endl;
    return;
  }

  if(dmrginp.transition_diff_irrep()){
    // It is used when bra and ket has different spatial irrep
    // For now, only the transtion pdm between two wavefuntions( 1 as bra and 0 as ket) are calculation
    // If the spatial irrep information is stored in wavefuntions, transition pdm among  several wavefunctions( i for bra and j for ket, there are n(n-1)/2 kinds of situations.) are possible.
    for (int state=0; state<dmrginp.nroots(); state++) {
      for(int stateB=0; stateB<state; stateB++){
        Timer timer;
        dmrginp.set_fullrestart() = true;
        sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
				dmrginp.npdm_generate() = true;
        SweepGenblock::do_one(sweepParams, false, !direction, false, 0, state, stateB); //this will generate the cd operators                               
				dmrginp.npdm_generate() = false;
        ecpu = timer.elapsedcputime();ewall=timer.elapsedwalltime();
        p3out << "\t\t\t NPDM SweepGenblock time " << ewall << " " << ecpu << endl;
        dmrginp.set_fullrestart() = false;

        sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
        Timer timerX;
        npdm_driver->clear();
        npdm_do_one_sweep(*npdm_driver, sweepParams, false, direction, false, 0, state,stateB);
        ecpu = timerX.elapsedcputime();ewall=timerX.elapsedwalltime();
        p3out << "\t\t\t NPDM sweep time " << ewall << " " << ecpu << endl;
       }
       }
    }

    else if( !dmrginp.setStateSpecific()){
    Timer timer;
    dmrginp.set_fullrestart() = true;
		dmrginp.npdm_generate() = true;
    sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
    SweepGenblock::do_one(sweepParams, false, !direction, false, 0, -1, -1); //this will generate the cd operators                               
		dmrginp.npdm_generate() = false;
    ecpu = timer.elapsedcputime();ewall=timer.elapsedwalltime();
    p3out << "\t\t\t NPDM SweepGenblock time " << ewall << " " << ecpu << endl;
    dmrginp.set_fullrestart() = false;

    
    if(transitionpdm){
      //  <\Phi_k|a^+_ia_j|\Phi_l> = <\Phi_l|a^+_ja_i|\Phi_k>*
      //  Therefore, only calculate the situations with k >= l.
      //  for(int stateB=0; stateB<= dmrginp.nroots(); stateB++){
      for (int state=0; state<dmrginp.nroots(); state++) {
         for(int stateB=0; stateB<=state; stateB++){
          sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
          Timer timerX;
          npdm_driver->clear();
          npdm_do_one_sweep(*npdm_driver, sweepParams, false, direction, false, 0, state,stateB);
          ecpu = timerX.elapsedcputime();ewall=timerX.elapsedwalltime();
          p3out << "\t\t\t NPDM sweep time " << ewall << " " << ecpu << endl;
          }
      }

    }
    else {
      for (int state=0; state<dmrginp.nroots(); state++) {
        sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
        if ( dmrginp.new_npdm_code() ) {
          Timer timerX;
          npdm_driver->clear();
          npdm_do_one_sweep(*npdm_driver, sweepParams, false, direction, false, 0, state,state);
          ecpu = timerX.elapsedcputime();ewall=timerX.elapsedwalltime();
          p3out << "\t\t\t NPDM sweep time " << ewall << " " << ecpu << endl;
        } 
        else{
          if (npdm_order == NPDM_ONEPDM) SweepOnepdm::do_one(sweepParams, false, direction, false, 0, state);     
          else if (npdm_order == NPDM_TWOPDM) SweepTwopdm::do_one(sweepParams, false, direction, false, 0, state);
          else abort();
        }
       }
       }

    }

  else {
    // state-specific
    if(transitionpdm){
      for (int state=0; state<dmrginp.nroots(); state++) {
        for (int stateB=0; stateB<=state; stateB++) {
  
        dmrginp.set_fullrestart() = true;
        sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
  
        if (mpigetrank() == 0) {
          Sweep::InitializeStateInfo(sweepParams, direction, state);
          Sweep::InitializeStateInfo(sweepParams, !direction, state);
          Sweep::CanonicalizeWavefunction(sweepParams, direction, state);
          Sweep::CanonicalizeWavefunction(sweepParams, !direction, state);
          Sweep::CanonicalizeWavefunction(sweepParams, direction, state);
        }
  
        sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
        if (mpigetrank() == 0) {
          Sweep::InitializeStateInfo(sweepParams, direction, stateB);
          Sweep::InitializeStateInfo(sweepParams, !direction, stateB);
          Sweep::CanonicalizeWavefunction(sweepParams, direction, stateB);
          Sweep::CanonicalizeWavefunction(sweepParams, !direction, stateB);
          Sweep::CanonicalizeWavefunction(sweepParams, direction, stateB);
        }
        // Prepare NPDM operators
        Timer timer;
        sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
				dmrginp.npdm_generate() = true;
        SweepGenblock::do_one(sweepParams, false, !direction, false, 0, state, stateB); //this will generate the cd operators
				dmrginp.npdm_generate() = false;
        dmrginp.set_fullrestart() = false;
        ecpu = timer.elapsedcputime();ewall=timer.elapsedwalltime();
        p3out << "\t\t\t NPDM SweepGenblock time " << ewall << " " << ecpu << endl;
  
        Timer timerX;
        npdm_driver->clear();
        sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
        npdm_do_one_sweep(*npdm_driver, sweepParams, false, direction, false, 0, state,stateB);
        ecpu = timerX.elapsedcputime();ewall=timerX.elapsedwalltime();
        p3out << "\t\t\t NPDM sweep time " << ewall << " " << ecpu << endl;
      }
    }
  }
  else{
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
			dmrginp.npdm_generate() = true;
      SweepGenblock::do_one(sweepParams, false, !direction, false, 0, state, state); //this will generate the cd operators
			dmrginp.npdm_generate() = false;
      dmrginp.set_fullrestart() = false;
      // Do NPDM sweep
      npdm_driver->clear();
      npdm_do_one_sweep(*npdm_driver, sweepParams, false, direction, false, 0, state,state);
   }
  }
  
}

  sweep_copy.savestate(direction_copy, restartsize_copy);

}
}
}
}

