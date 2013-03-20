/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "modules/lrt/lrt_sweep.h"
#include "global.h"
#include "initblocks.h"
#include "MatrixBLAS.h"
#include "modules/lrt/lrt_initblocks.h"
#include "modules/lrt/lrt_transform_gauge.h"
#include "modules/lrt/lrt_davidson.h"
#include "modules/lrt/lrt_solver.h"
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
 const vector<double>& eigenvalues, vector<double>& rnorm, vector<double>& ynorm,
 Matrix& a_subspace, Matrix& b_subspace, Matrix& s_subspace, Matrix& d_subspace, const Matrix& alpha,
 const bool &useSlater, const bool& dot_with_sys, const bool& rpa_sweep, const bool& rpa_sweep_2nd, int nroots, int mroots, int kroots)
{
  const int lroots = mroots + nroots - kroots;

  int Nroots = (rpa_sweep ? 2 * nroots - 1 : nroots);
  int Mroots = (rpa_sweep ? 2 * mroots - 1 : mroots);
  int Lroots = (rpa_sweep ? 2 * lroots - 1 : lroots);

  pout << "\t\t\t Davidson info: nroots = " << nroots << ", trial vecs = " << mroots << ", converged = " << kroots << endl;

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
    // broadcast additional compops (0-th)
    system.addAdditionalCompOps();
    // broadcast additional compops (1-st)
    for(int i = 1; i < Lroots; ++i) {
      system.addAdditionalCompOps(0, i);
      system.addAdditionalCompOps(i, 0);
    }
    dmrginp.datatransfer -> stop();
    if (dot_with_sys) {
      InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, sweepParams.get_sys_add(),
                                     dmrginp.direct(), DISTRIBUTED_STORAGE, dot_with_sys, true, Lroots);

    }
    InitBlocks::LRT::InitNewEnvironmentBlock(environment, systemDot, newEnvironment, system, systemDot, alpha,
                                             sweepParams.get_sys_add(), sweepParams.get_env_add(), forward, dmrginp.direct(),
                                             sweepParams.get_onedot(), nexact, useSlater, !dot_with_sys, true, dot_with_sys, rpa_sweep_2nd, Mroots);
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

  if(!rpa_sweep_2nd)
    GuessWave::LRT::rotate_previous_wavefunction(big, alpha, Mroots, sweepParams.get_guesstype(), sweepParams.get_onedot(), dot_with_sys);

  bool last_site = (sweepParams.get_block_iter() == sweepParams.get_n_iters() - 1) ? true : false;

  newSystem.RenormaliseFrom_lrt(eigenvalues, rnorm, ynorm, rotatematrices, nroots, mroots, kroots,
                                a_subspace, b_subspace, s_subspace, d_subspace,
                                sweepParams.get_keep_states(), sweepParams.get_keep_qstates(), big,
                                sweepParams.get_guesstype(), sweepParams.get_onedot(), last_site, rpa_sweep, rpa_sweep_2nd,
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
(SweepParams &sweepParams, const bool& warmUp, const bool& forward, const bool& rpa_sweep_2nd, const bool& restart, const int& restartSize)
{
  SpinBlock system;
  const int nroots = dmrginp.nroots(sweepParams.get_sweep_iter());
  double zerothEnergy;
  std::vector<double> finalEnergy(nroots, 0.);
  std::vector<double> finalEnergy_spins(nroots, 0.);
//double finalError = 0.;
  zerothEnergy = sweepParams.get_lowest_energy()[0];
  finalEnergy_spins.resize(nroots, sweepParams.get_lowest_energy()[0]);
//finalError = sweepParams.get_lowest_error();

  //====================================================================================================
  // Load and solve generalized eigenvalue problem: H x = E S x
  //====================================================================================================

#ifndef SERIAL
  mpi::communicator world;
#endif

  bool rpa_sweep = (dmrginp.lrt_type() == RPA);
  bool deflation_sweep;
  int mroots, i_conv_root;
  Matrix a_subspace, b_subspace, s_subspace, d_subspace;
  std::vector<double> eigenvalues;
  std::vector<double> rnorm;
  std::vector<double> ynorm;
  Matrix alpha;

  // here, warmUp is used to detect if it's the initial sweep
  // guess wavefunctions will be computed from Krylov subspace
  if(mpigetrank() == 0) {
    if(!warmUp) {
      SpinAdapted::LRT::LoadDavidsonInfo(a_subspace, b_subspace, s_subspace, d_subspace,
                                         eigenvalues, rnorm, ynorm, mroots, i_conv_root, deflation_sweep);

      if(dmrginp.outputlevel() > 0 && !rpa_sweep_2nd) {
        pout << "\t\t\t printing reduced A-matrix: " << a_subspace.Nrows() << " x " << a_subspace.Ncols() << endl;
        for(int i = 1; i < mroots; ++i) {
          pout << "\t ";
          for(int j = 1; j < mroots; ++j) {
            pout << setw(16) << fixed << setprecision(8) << a_subspace(i, j);
          }
          pout << endl;
        }
        pout << endl;

        if(rpa_sweep) {

        pout << "\t\t\t printing reduced B-matrix: " << b_subspace.Nrows() << " x " << b_subspace.Ncols() << endl;
        for(int i = 1; i < mroots; ++i) {
          pout << "\t ";
          for(int j = 1; j < mroots; ++j) {
            pout << setw(16) << fixed << setprecision(8) << b_subspace(i, j);
          }
          pout << endl;
        }
        pout << endl;

        }

        pout << "\t\t\t printing reduced S-matrix: " << s_subspace.Nrows() << " x " << s_subspace.Ncols() << endl;
        for(int i = 1; i < mroots; ++i) {
          pout << "\t ";
          for(int j = 1; j < mroots; ++j) {
            pout << setw(16) << fixed << setprecision(8) << s_subspace(i, j);
          }
          pout << endl;
        }
        pout << endl;

        if(rpa_sweep) {

        pout << "\t\t\t printing reduced D-matrix: " << d_subspace.Nrows() << " x " << d_subspace.Ncols() << endl;
        for(int i = 1; i < mroots; ++i) {
          pout << "\t ";
          for(int j = 1; j < mroots; ++j) {
            pout << setw(16) << fixed << setprecision(8) << d_subspace(i, j);
          }
          pout << endl;
        }
        pout << endl;

        }
      }

      if(!rpa_sweep_2nd) {
        DiagonalMatrix freq;
        if(rpa_sweep)
          mroots = 1 + SpinAdapted::LRT::compute_eigenvalues(a_subspace, b_subspace, s_subspace, d_subspace, freq, alpha);
        else
          mroots = 1 + SpinAdapted::LRT::compute_eigenvalues(a_subspace, s_subspace, freq, alpha);

        for(int i = 1; i < nroots; ++i) {
          eigenvalues[i] = freq(i, i); // excitation frequencies
        }

        if(dmrginp.outputlevel() > 0) {
          pout << "\t\t\t printing frequencies" << endl;
          pout << "\t ";
          for(int i = 1; i < eigenvalues.size(); ++i) {
            pout << setw(16) << fixed << setprecision(8) << eigenvalues[i];
          }
          pout << endl << endl;
        }
        if(dmrginp.outputlevel() > 1) {
          pout << "\t\t\t printing state rotation vectors" << endl;
          for(int i = 0; i < alpha.Nrows(); ++i) {
            pout << "\t ";
            for(int j = 0; j < alpha.Ncols(); ++j) {
              pout << setw(16) << fixed << setprecision(8) << alpha(i+1,j+1);
            }
            pout << endl;
          }
        }
      }
// DEBUG: make no correction vectors
//    mroots = nroots;
//    i_conv_root = nroots;
//    deflation_sweep = false;
// DEBUG: end
    }
    else {
      if(rpa_sweep_2nd) {
        SpinAdapted::LRT::LoadDavidsonInfo(a_subspace, b_subspace, s_subspace, d_subspace,
                                           eigenvalues, rnorm, ynorm, mroots, i_conv_root, deflation_sweep);
      }
      else {
        // compute guesses from Krylov subspace
        eigenvalues.resize(nroots, 0.0);
        mroots = 1;
        i_conv_root = 1;
        deflation_sweep = false;
      }
    }
    eigenvalues[0] = zerothEnergy; // <-- 0-th energy
  }

#ifndef SERIAL
  mpi::broadcast(world, eigenvalues, 0);
  mpi::broadcast(world, rnorm, 0);
  mpi::broadcast(world, ynorm, 0);
  mpi::broadcast(world, mroots, 0);
  mpi::broadcast(world, i_conv_root, 0);
  mpi::broadcast(world, deflation_sweep, 0);
  mpi::broadcast(world, alpha, 0);
#endif

  if (deflation_sweep) mroots = nroots;

  int lroots = mroots + nroots - i_conv_root;

  if(!rpa_sweep_2nd) {
    a_subspace.ReSize(lroots-1, lroots-1); a_subspace = 0.0;
    b_subspace.ReSize(lroots-1, lroots-1); b_subspace = 0.0;
    s_subspace.ReSize(lroots-1, lroots-1); s_subspace = 0.0;
    d_subspace.ReSize(lroots-1, lroots-1); d_subspace = 0.0;
    rnorm.resize(nroots); fill(rnorm.begin(), rnorm.end(), 0.0);
    ynorm.resize(lroots); fill(ynorm.begin(), ynorm.end(), 0.0);
  }

  //====================================================================================================
  // Sweep Iteration for Davidson update
  //====================================================================================================

  sweepParams.set_sweep_parameters();
  // a new renormalisation sweep routine
  pout << endl;
  if (forward) {
    pout << "\t\t\t Starting sweep " << sweepParams.set_sweep_iter() << " in forwards direction" << endl;
  }
  else {
    pout << "\t\t\t Starting sweep " << sweepParams.set_sweep_iter() << " in backwards direction" << endl;
  }
  pout << "\t\t\t ============================================================================ " << endl;

  InitBlocks::InitStartingBlock (system, forward, sweepParams.get_forward_starting_size(),
                                                  sweepParams.get_backward_starting_size(), restartSize, restart, warmUp);
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

  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); ) {

    pout << "\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
    pout << "\t\t\t ----------------------------" << endl;
    if (dmrginp.outputlevel() > 0) {
      if (forward) {
        pout << "\t\t\t Current direction is :: Forwards " << endl;
      }
      else {
        pout << "\t\t\t Current direction is :: Backwards " << endl;
      }
    }

    if (sweepParams.get_block_iter() == 0) {
      sweepParams.set_guesstype() = TRANSPOSE;
    }
    else {
      sweepParams.set_guesstype() = TRANSFORM;
    }
    
//pout << "DEBUG @ Sweep::LRT::do_one: before BlockAndDecimate, printing system operators" << endl;
//  if (dmrginp.outputlevel() > 0){
//     pout << system<<endl;
//     system.printOperatorSummary();
//  }

    if (dmrginp.outputlevel() > 0)
       pout << "\t\t\t Blocking and Decimating " << endl;
        
    SpinBlock newSystem;

    BlockAndDecimate (sweepParams, system, newSystem, eigenvalues, rnorm, ynorm, a_subspace, b_subspace, s_subspace, d_subspace, alpha,
                      false, dot_with_sys, rpa_sweep, rpa_sweep_2nd, nroots, mroots, i_conv_root);

    if(dmrginp.outputlevel() > 2) { // this might give verbose output
      pout << "\t\t\t printing reduced A-matrix: " << a_subspace.Nrows() << " x " << a_subspace.Ncols() << endl;
      for(int i = 1; i < lroots; ++i) {
        pout << "\t ";
        for(int j = 1; j < lroots; ++j) {
          pout << setw(16) << fixed << setprecision(8) << a_subspace(i, j);
        }
        pout << endl;
      }
      pout << endl;

      if(rpa_sweep) {

      pout << "\t\t\t printing reduced B-matrix: " << b_subspace.Nrows() << " x " << b_subspace.Ncols() << endl;
      for(int i = 1; i < lroots; ++i) {
        pout << "\t ";
        for(int j = 1; j < lroots; ++j) {
          pout << setw(16) << fixed << setprecision(8) << b_subspace(i, j);
        }
        pout << endl;
      }
      pout << endl;

      }

      pout << "\t\t\t printing reduced S-matrix: " << s_subspace.Nrows() << " x " << s_subspace.Ncols() << endl;
      for(int i = 1; i < lroots; ++i) {
        pout << "\t ";
        for(int j = 1; j < lroots; ++j) {
          pout << setw(16) << fixed << setprecision(8) << s_subspace(i, j);
        }
        pout << endl;
      }
      pout << endl;

      if(rpa_sweep) {

      pout << "\t\t\t printing reduced D-matrix: " << d_subspace.Nrows() << " x " << d_subspace.Ncols() << endl;
      for(int i = 1; i < lroots; ++i) {
        pout << "\t ";
        for(int j = 1; j < lroots; ++j) {
          pout << setw(16) << fixed << setprecision(8) << d_subspace(i, j);
        }
        pout << endl;
      }
      pout << endl;

      }
    }

    system = newSystem;

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
    world.barrier();
#endif
    sweepParams.savestate(forward, syssites.size());
    if (dmrginp.outputlevel() > 0)
       mcheck("at the end of sweep iteration");

  } // end block iter

  // compute RPA energy
  double rpaEcorr = 0.0;
  if(mpigetrank() == 0) {
    for(int i = 1; i < nroots; ++i)
      rpaEcorr -= eigenvalues[i] * ynorm[i];
  }

#ifndef SERIAL
  mpi::broadcast(world, rpaEcorr, 0);
  mpi::broadcast(world, rnorm, 0);
  mpi::broadcast(world, ynorm, 0);
  mpi::broadcast(world, eigenvalues, 0);
#endif

  finalEnergy[0] = eigenvalues[0] - rpaEcorr;
  for(int i = 1; i < nroots; ++i)
    finalEnergy[i] = finalEnergy[0] + eigenvalues[i];

  // check convergence and save iteration info
  if (mpigetrank() == 0) {
    if(!warmUp && (!rpa_sweep || rpa_sweep_2nd)) {
      i_conv_root = 1;
      for(; i_conv_root < nroots; ++i_conv_root) {
        if(rnorm[i_conv_root] > sweepParams.get_davidson_tol()) break;
      }
    }

    if(!rpa_sweep || rpa_sweep_2nd) mroots = lroots;

    deflation_sweep = (mroots > dmrginp.deflation_max_size()) ? true : false;

    SpinAdapted::LRT::SaveDavidsonInfo(a_subspace, b_subspace, s_subspace, d_subspace,
                                       eigenvalues, rnorm, ynorm, mroots, i_conv_root, deflation_sweep);
  }

  if(!warmUp && (!rpa_sweep || rpa_sweep_2nd)) {
    if(dmrginp.outputlevel() > 0) {
      pout << "\t\t\t printing eigenvalues of effective hamiltonian" << endl;
      for(int j = 1; j < nroots; ++j) {
        pout << "\t\t\t excitation energy for State [ " << j << " ] with Spin [ " << dmrginp.molecule_quantum().get_s()
             << " ] :: " << setw(16) << fixed << setprecision(12) << eigenvalues[j]
             << " ( Y-norm  " << setw(20) << scientific << setprecision(12) << ynorm[j] << " ) " << endl;
      }
    }

    pout << "\t\t\t ============================================================================ " << endl;
    pout << "\t\t\t energy for State [ 0 ] with Spin [ "
         << dmrginp.molecule_quantum().get_s()  << " ] :: "
         << setw(16) << fixed << setprecision(12) << finalEnergy[0]+dmrginp.get_coreenergy() << endl;
    if(rpa_sweep)
      pout << "\t\t\t ( RPA correlation energy :: " << setw(16) << fixed << setprecision(12) << rpaEcorr << " ) " << endl;

    for(int j = 1; j < nroots; ++j) {
      pout << "\t\t\t Finished Sweep with " << sweepParams.get_keep_states() << " states and energy for State [ " << j 
           << " ] with Spin [ " << dmrginp.molecule_quantum().get_s()  << " ] :: "
           << setw(16) << fixed << setprecision(12) << finalEnergy[j]+dmrginp.get_coreenergy()
           << " ( R-norm  " << scientific << setprecision(3) << rnorm[j] << " ) " << endl;
    }
  }

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

  double max_rnorm;
  if(mpigetrank() == 0) {
    max_rnorm = *max_element(rnorm.begin()+1, rnorm.end());
  }
#ifndef SERIAL
  mpi::broadcast(world, max_rnorm, 0);
#endif

  return max_rnorm;
}

