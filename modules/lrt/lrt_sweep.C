/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "modules/lrt/lrt_sweep.h"
#include "global.h"
#include "initblocks.h"
#include "modules/lrt/lrt_initblocks.h"
#include "modules/lrt/lrt_transform_gauge.h"
#include "modules/lrt/lrt_davidson.h"
#include "MatrixBLAS.h"
#include <boost/format.hpp>

#ifndef SERIAL
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif

//#include "rotationmat.h"
//#include "density.h"
#include "pario.h"

using namespace boost;
using namespace std;


void SpinAdapted::Sweep::LRT::BlockAndDecimate
(SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem,
 const vector<double>& eigenvalues, vector<double>& rnorm, Matrix& h_subspace, Matrix& s_subspace, const Matrix& alpha,
 const bool &useSlater, const bool& dot_with_sys, int nroots, int mroots, int kroots, const bool& deflation_sweep)
{
  const int lroots = mroots + nroots - kroots;

pout << "DEBUG @ Sweep::LRT::BlockAndDecimate: called, nroots = " << nroots << ", mroots = " << mroots << ", kroots = " << kroots << ", lroots = " << lroots << ", " << (deflation_sweep ? "deflation" : "normal") << endl;
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

  const int nexact = forward ? sweepParams.get_forward_starting_size() : sweepParams.get_backward_starting_size();

  //before halfway put the sysdot with system otherwise with environment
  if (sweepParams.get_onedot()) {
    dmrginp.datatransfer -> start();
    system.addAdditionalCompOps();
    dmrginp.datatransfer -> stop();
    if (dot_with_sys) {
      InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, sweepParams.get_sys_add(),
                                     dmrginp.direct(), DISTRIBUTED_STORAGE, dot_with_sys, true, lroots);

    }
    InitBlocks::LRT::InitNewEnvironmentBlock(environment, systemDot, newEnvironment, system, systemDot, alpha,
                                             sweepParams.get_sys_add(), sweepParams.get_env_add(), forward, dmrginp.direct(),
                                             sweepParams.get_onedot(), nexact, useSlater, !dot_with_sys, true, dot_with_sys, mroots);
  }
  else {
    pout << "\t\t\t DMRG-LRT calculation can only be performed upon onedot wavefunction" << endl;
    abort();
  }

  SpinBlock big;
  if (dot_with_sys) {
    newSystem.set_loopblock(true);
    system.set_loopblock(false);
    newEnvironment.set_loopblock(false);
    InitBlocks::InitBigBlock(newSystem, newEnvironment, big); 
  }
  else{
    system.set_loopblock(false);
    newEnvironment.set_loopblock(true);
    environment.set_loopblock(true);
    InitBlocks::InitBigBlock(system, newEnvironment, big); 
  }
  //analyse_operator_distribution(big);
  dmrginp.guessgenT -> stop();
  dmrginp.multiplierT -> start();
  std::vector< std::vector<Matrix> > rotatematrices;

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

//pout << "DEBUG @ Sweep::LRT::BlockAndDecimate: print built operators" << endl;
//  if (!dot_with_sys && sweepParams.get_onedot()) {
//pout << "DEBUG @ Sweep::LRT::BlockAndDecimate: system operators" << endl;
//    system.printOperatorSummary();
//pout << "DEBUG @ Sweep::LRT::BlockAndDecimate: environment operators" << endl;
//    newEnvironment.printOperatorSummary();
//  }
//  else {
//pout << "DEBUG @ Sweep::LRT::BlockAndDecimate: system operators" << endl;
//    newSystem.printOperatorSummary();
//pout << "DEBUG @ Sweep::LRT::BlockAndDecimate: environment operators" << endl;
//    newEnvironment.printOperatorSummary();
//  }

  GuessWave::LRT::rotate_previous_wavefunction(big, alpha, mroots, sweepParams.get_guesstype(),
                                               sweepParams.get_onedot(), dot_with_sys);
  if(deflation_sweep) mroots = nroots;

  bool last_site = false;
  if(sweepParams.get_block_iter() == sweepParams.get_n_iters() - 1)
       last_site = true;

  newSystem.RenormaliseFrom_lrt(eigenvalues, rnorm, rotatematrices, nroots, mroots, kroots, h_subspace, s_subspace,
                                sweepParams.get_keep_states(), sweepParams.get_keep_qstates(), big,
                                sweepParams.get_guesstype(), sweepParams.get_onedot(), last_site,
                                system, systemDot, environmentDot, environment, dot_with_sys, sweepParams.get_sweep_iter());

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

  newSystem.transform_operators_lrt(rotatematrices);
  newSystem.transform_operators(rotatematrices[0]);

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

double SpinAdapted::Sweep::LRT::do_one
(SweepParams &sweepParams, const bool& warmUp, const bool& forward, const bool& restart, const int& restartSize)
{
  SpinBlock system;
  const int nroots = dmrginp.nroots(sweepParams.get_sweep_iter());
  std::vector<double> finalEnergy(nroots, 0.);
  std::vector<double> finalEnergy_spins(nroots, 0.);
  double finalError = 0.;
  if (restart) {
    finalEnergy[0]      = sweepParams.get_lowest_energy()[0];
    finalEnergy_spins[0]= sweepParams.get_lowest_energy()[0];
    finalError = sweepParams.get_lowest_error();
  }

  // FIXME: add here, initializations of DMRG-LRT calculation?
#ifndef SERIAL
  mpi::communicator world;
#endif

  bool deflation_sweep;
  int mroots, i_conv_root;
  Matrix h_subspace, s_subspace;
  std::vector<double> eigenvalues(nroots);
  Matrix alpha;

  // here, warmUp is used to detect if it's the initial sweep
  // guess wavefunctions will be computed from Krylov subspace
  if(mpigetrank() == 0) {
    eigenvalues[0] = 0.0; // <-- 0-th energy, but not necessary (TODO)
    if(!warmUp) {
      SpinAdapted::LRT::LoadDavidsonInfo(h_subspace, s_subspace, mroots, i_conv_root, deflation_sweep);
//pout << "DEBUG @ Sweep::LRT::do_one: printing subspace H" << endl;
//      for(int i = 1; i < mroots; ++i) {
//        pout << "\t\t\t ";
//        for(int j = 1; j < mroots; ++j) {
//          pout << setw(16) << fixed << setprecision(8) << h_subspace(i, j);
//        }
//        pout << endl;
//      }
//pout << "DEBUG @ Sweep::LRT::do_one: printing subspace S" << endl;
//      for(int i = 1; i < mroots; ++i) {
//        pout << "\t\t\t ";
//        for(int j = 1; j < mroots; ++j) {
//          pout << setw(16) << fixed << setprecision(8) << s_subspace(i, j);
//        }
//        pout << endl;
//      }

      DiagonalMatrix ritzval;
      diagonalise(h_subspace, s_subspace, ritzval, alpha);
      for(int i = 1; i < nroots; ++i) {
        eigenvalues[i] = ritzval(i, i);
      }
//pout << "DEBUG @ Sweep::LRT::do_one: Ritz vectors" << endl;
//      for(int i = 0; i < alpha.Nrows(); ++i) {
//        pout << "\t";
//        for(int j = 0; j < alpha.Ncols(); ++j) {
//          pout << setw(16) << fixed << setprecision(8) << alpha(i+1,j+1);
//        }
//        pout << endl;
//      }
// DEBUG: make no correction vectors
//    mroots = nroots;
//    i_conv_root = nroots;
//    deflation_sweep = false;
// DEBUG: end
    }
    else {
      // compute guesses from Krylov subspace
      mroots = 1;
      i_conv_root = 1;
      deflation_sweep = false;
    }
  }

#ifndef SERIAL
  mpi::broadcast(world, deflation_sweep, 0);
  mpi::broadcast(world, mroots, 0);
  mpi::broadcast(world, i_conv_root, 0);
  mpi::broadcast(world, eigenvalues, 0);
  mpi::broadcast(world, alpha, 0);
#endif

  for(int i = 1; i < nroots; ++i) finalEnergy[i] = eigenvalues[i];
  std::vector<double> rnorm(nroots, 0.0);

  if (deflation_sweep) mroots = nroots;

  int lroots = mroots + nroots - i_conv_root;
  h_subspace.ReSize(lroots-1, lroots-1); h_subspace = 0.0;
  s_subspace.ReSize(lroots-1, lroots-1); s_subspace = 0.0;

  sweepParams.set_sweep_parameters();
  // a new renormalisation sweep routine
  pout << endl;
  if (forward)
    pout << "\t\t\t Starting sweep "<< sweepParams.set_sweep_iter()<<" in forwards direction"<<endl;
  else
    pout << "\t\t\t Starting sweep "<< sweepParams.set_sweep_iter()<<" in backwards direction" << endl;
  pout << "\t\t\t ============================================================================ " << endl;

  InitBlocks::InitStartingBlock (system, forward, sweepParams.get_forward_starting_size(),
                                                  sweepParams.get_backward_starting_size(), restartSize, restart, warmUp, 1);
// FIXME: not sure whether starting block needs to have state indices
//                                                sweepParams.get_backward_starting_size(), restartSize, restart, warmUp, lroots);
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

//if (restart)
//{
//  if (forward && system.get_complementary_sites()[0] >= dmrginp.last_site()/2)
//    dot_with_sys = false;
//  if (!forward && system.get_sites()[0]-1 < dmrginp.last_site()/2)
//    dot_with_sys = false;
//}
  if (dmrginp.outputlevel() > 0)
    mcheck("at the very start of sweep");

  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); ) {

    pout << "\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
    pout << "\t\t\t ----------------------------" << endl;
    if (dmrginp.outputlevel() > 0) {
      if (forward)
        pout << "\t\t\t Current direction is :: Forwards " << endl;
      else
        pout << "\t\t\t Current direction is :: Backwards " << endl;
    }

//  if (dmrginp.no_transform() || (sweepParams.get_sweep_iter()-sweepParams.get_restart_iter() == 0 && sweepParams.get_block_iter() == 0))
//    sweepParams.set_guesstype() = BASIC;
//  else if (!warmUp && sweepParams.get_block_iter() != 0) 
//    sweepParams.set_guesstype() = TRANSFORM;
//  else if (!warmUp && sweepParams.get_block_iter() == 0 && 
//          ((dmrginp.algorithm_method() == TWODOT_TO_ONEDOT && dmrginp.twodot_to_onedot_iter() != sweepParams.get_sweep_iter()) ||
//            dmrginp.algorithm_method() != TWODOT_TO_ONEDOT))
//    sweepParams.set_guesstype() = TRANSPOSE;
//  else
//    sweepParams.set_guesstype() = BASIC;

//  if (!warmUp && sweepParams.get_block_iter() == 0) {
    if (sweepParams.get_block_iter() == 0) {
pout << "DEBUG @ Sweep::LRT::do_one: check point 8-1 (guess = transpose)" << endl;
      sweepParams.set_guesstype() = TRANSPOSE;
    }
//  else if (!warmUp && sweepParams.get_block_iter() != 0) {
    else {
pout << "DEBUG @ Sweep::LRT::do_one: check point 8-2 (guess = transform)" << endl;
      sweepParams.set_guesstype() = TRANSFORM;
    }
//  else {
//    pout << "\t\t\t ERROR: guess type BASIC was requested, this makes a problem in the gauge-transformation" << endl;
//    abort();
//  }
    
    if (dmrginp.outputlevel() > 0)
       pout << "\t\t\t Blocking and Decimating " << endl;
        
    SpinBlock newSystem;

    BlockAndDecimate (sweepParams, system, newSystem, eigenvalues, rnorm, h_subspace, s_subspace, alpha,
                      false, dot_with_sys, nroots, mroots, i_conv_root, deflation_sweep);

//pout << "DEBUG @ Sweep::LRT::do_one: printing subspace H" << endl;
//    for(int i = 1; i < lroots; ++i) {
//      pout << "\t\t\t ";
//      for(int j = 1; j < lroots; ++j) {
//        pout << setw(16) << fixed << setprecision(8) << h_subspace(i, j);
//      }
//      pout << endl;
//    }
//pout << "DEBUG @ Sweep::LRT::do_one: printing subspace S" << endl;
//    for(int i = 1; i < lroots; ++i) {
//      pout << "\t\t\t ";
//      for(int j = 1; j < lroots; ++j) {
//        pout << setw(16) << fixed << setprecision(8) << s_subspace(i, j);
//      }
//      pout << endl;
//    }

    system = newSystem;
//  if (dmrginp.outputlevel() > 0){
//     pout << system<<endl;
//     system.printOperatorSummary();
//  }

//  //system size is going to be less than environment size
//  if (forward && system.get_complementary_sites()[0] >= dmrginp.last_site()/2)
//    dot_with_sys = false;
//  if (!forward && system.get_sites()[0]-1 < dmrginp.last_site()/2)
//    dot_with_sys = false;

    SpinBlock::store (forward, system.get_sites(), system);                 
    syssites = system.get_sites();
    if (dmrginp.outputlevel() > 0)
       pout << "\t\t\t saving state " << syssites.size() << endl;
    ++sweepParams.set_block_iter();
      
#ifndef SERIAL
    world.barrier();
#endif
    sweepParams.savestate(forward, syssites.size());
    if (dmrginp.outputlevel() > 0)
       mcheck("at the end of sweep iteration");

  } // end block iter

#ifndef SERIAL
  mpi::broadcast(world, rnorm, 0);
#endif

  // check convergence and save iteration info
  if (mpigetrank() == 0) {
    i_conv_root = 1;
    if(!warmUp) {
      for(; i_conv_root < nroots; ++i_conv_root) {
        if(rnorm[i_conv_root] > sweepParams.get_davidson_tol()) break;
      }
    }

    mroots = lroots;
    deflation_sweep = mroots > dmrginp.deflation_max_size() ? true : false;

    SpinAdapted::LRT::SaveDavidsonInfo(h_subspace, s_subspace, mroots, i_conv_root, deflation_sweep);
  }

  for(int j = 0; j < nroots; ++j) {
    pout << "\t\t\t Finished Sweep with " << sweepParams.get_keep_states() << " states and sweep energy for State [ " << j 
         << " ] with Spin [ " << dmrginp.molecule_quantum().get_s()  << " ] :: " << fixed << setprecision(8) << finalEnergy[j]+dmrginp.get_coreenergy()
         << " ( R-norm  " << scientific << setprecision(3) << rnorm[j] << " ) " << endl;
  }

  pout << "\t\t\t Largest Error for Sweep with " << sweepParams.get_keep_states() << " states is " << finalError << endl;
  sweepParams.set_largest_dw() = finalError;
  pout << "\t\t\t ============================================================================ " << endl;

  // update the static number of iterations

  ++sweepParams.set_sweep_iter();
  if (!mpigetrank()) {
#ifndef MOLPRO
    FILE* f = fopen("dmrg.e", "wb");
#else
    std::string efile;
    efile = str(boost::format("%s%s") % dmrginp.load_prefix() % "/dmrg.e" );
    FILE* f = fopen(efile.c_str(), "wb");
#endif
    
    for(int j = 0; j < nroots; ++j) {
      double e = finalEnergy[j]+dmrginp.get_coreenergy(); 
      fwrite( &e, 1, sizeof(double), f);
    }
    fclose(f);
  }

  return *max_element(rnorm.begin(), rnorm.end());
}

