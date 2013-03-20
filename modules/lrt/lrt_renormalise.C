/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "spinblock.h"
#include <boost/bind.hpp>
#include <boost/functional.hpp>
#include <boost/function.hpp>
#include <boost/make_shared.hpp>
#include "solver.h"
#include "modules/lrt/lrt_solver.h"
#include "operatorloops.h"
#include <numeric>
#include "Operators.h"
#include "wavefunction.h"
#include "rotationmat.h"
#include "modules/lrt/lrt_rotationmat.h"
#include "density.h"
#include "operatorfunctions.h"
#include "initblocks.h"
#include "guess_wavefunction.h"
#include "linear.h"
#include "davidson.h"
#include "modules/lrt/lrt_davidson.h"
#include <stdlib.h>
//#include "diis.h"

#ifndef SERIAL
#include <boost/mpi.hpp>
#endif

#include "pario.h"

using namespace boost;
using namespace std;

namespace SpinAdapted {
using namespace operatorfunctions;

void SpinBlock::RenormaliseFrom_lrt
(const vector<double> &energies, vector<double>& rnorm, vector<double>& ynorm, vector< vector<Matrix> >& rotateMatrices, int nroots, int mroots, int kroots,
 Matrix& a_subspace, Matrix& b_subspace, Matrix& s_subspace, Matrix& d_subspace, const int keptstates, const int keptqstates,
 SpinBlock& big, const guessWaveTypes &guesswavetype, const bool &onedot, const bool &last_site, const bool &rpa_sweep, const bool &rpa_sweep_2nd,
 SpinBlock& System, SpinBlock& sysDot, SpinBlock& envDot, SpinBlock& environment, const bool& dot_with_sys, int sweepiter)
{
  int lroots = mroots + nroots - kroots;

  int Nroots = (rpa_sweep ? 2 * nroots - 1 : nroots);
  int Mroots = (rpa_sweep ? 2 * mroots - 1 : mroots);
  int Lroots = (rpa_sweep ? 2 * lroots - 1 : lroots);

  vector<Wavefunction> wave_solutions_1st;
  vector<Wavefunction> wave_solutions_2nd;

  dmrginp.davidsonT -> start();
  if (dmrginp.outputlevel() > 0)
    mcheck("before davidson but after all blocks are built");

  dmrginp.solvewf -> start();

  //
  // NOTE: wave_solutions_1st[0] = wave_solutions_2nd[0] : 0-th wavefunction
  //

  LRT::solve_wavefunction(wave_solutions_1st, wave_solutions_2nd,
                          energies, rnorm, big, guesswavetype, onedot, dot_with_sys, rpa_sweep, rpa_sweep_2nd, nroots, mroots, kroots);

  if(rpa_sweep) {
    // compute 1-st order components for RPA
    // specialize to use dot with environment at last site (necessary for parallel run)
    if(!rpa_sweep_2nd && last_site) {
      LRT::multiply_h_left davidson_f(big, onedot);
      LRT::RPA::compute_matrix_elements(wave_solutions_1st, energies, ynorm, davidson_f,
                                        a_subspace, b_subspace, s_subspace, d_subspace, rpa_sweep_2nd, lroots);
    }
    // compute 2-nd order components for RPA
    else if(rpa_sweep_2nd && !last_site) {
      LRT::multiply_h_left davidson_f(big, onedot);
      LRT::RPA::compute_matrix_elements(wave_solutions_2nd, energies, ynorm, davidson_f,
                                        a_subspace, b_subspace, s_subspace, d_subspace, rpa_sweep_2nd, lroots);
    }
  }
  else {
    // compute 1-st order components for TDA
    // specialize to use dot with environment at last site (necessary for parallel run)
    if(last_site) {
      LRT::multiply_h_left davidson_f(big, onedot);
      LRT::TDA::compute_matrix_elements(wave_solutions_1st, energies, davidson_f, a_subspace, s_subspace, lroots);
    }
  }

  dmrginp.solvewf -> stop();
  SpinBlock newsystem;
  SpinBlock newenvironment;
  SpinBlock newbig;
  dmrginp.postwfrearrange -> start();

  // maybe better to turn 'dot with env' off for computing (1 - L(0)L(0)')C(I) later
  if (onedot && !dot_with_sys) {
    InitBlocks::InitNewSystemBlock(System, sysDot, newsystem, sysDot.size(), dmrginp.direct(), DISTRIBUTED_STORAGE, false, true, Lroots);
    InitBlocks::InitBigBlock(newsystem, environment, newbig); 
    for (int i = 0; i < Lroots && mpigetrank() == 0; i++) {
      Wavefunction tempwave = wave_solutions_1st[i];
      GuessWave::onedot_shufflesysdot(big.get_stateInfo(), newbig.get_stateInfo(), wave_solutions_1st[i], tempwave);  
      wave_solutions_1st[i] = tempwave;
    }
    *this = newsystem;
    if (dmrginp.outputlevel() > 0)
       cout << newsystem.get_twoInt().get()<<"  "<<get_twoInt().get()<<"  Ints "<<endl;
    envDot.clear();
    big.get_rightBlock()->clear();
    big.clear();
  }
  else
    newbig = big;

  dmrginp.postwfrearrange -> stop();

  if (dmrginp.outputlevel() > 0)
    mcheck("after davidson before noise");

  dmrginp.davidsonT -> stop();

  dmrginp.rotmatrixT -> start();

  DensityMatrix tracedMatrix;
  tracedMatrix.allocate(stateInfo);

  rotateMatrices.resize(Lroots);
  vector<Wavefunction> projected_wave_solutions(Lroots);

  vector< vector<double> > selectedwts;
  vector< vector<double> > rejectedwts;
  vector< Matrix > rejectedbasis;

  if(mpigetrank() == 0) {
    // FIXME: re-computing 0-th rotation matrix is wasteful, but computing 1-st rotation matrices needs its eigenvalues
    MultiplyProduct (wave_solutions_1st[0], Transpose(wave_solutions_1st[0]), tracedMatrix, 1.0);

    // find and sort weight info
    DensityMatrix transformmatrix;
    transformmatrix.allocate(stateInfo);
    std::vector<DiagonalMatrix> eigenMatrix;
    diagonalise_dm(tracedMatrix, transformmatrix, eigenMatrix);

    vector<pair<int, int> > inorderwts;
    vector<vector<int> > wtsbyquanta;
  
    int sys_dot_size = *get_sites().rbegin ()+1 ;
    sort_weights(eigenMatrix, inorderwts, wtsbyquanta);

    // make transformation matrix by various algorithms
    int totalstatesbydm = min(static_cast<int>(inorderwts.size()), keptstates);
    int totalstatesbyquanta = min(static_cast<int>(inorderwts.size()), keptstates + keptqstates) - totalstatesbydm;
    if (totalstatesbyquanta < 0) totalstatesbyquanta = 0;
  
    if (dmrginp.outputlevel() > 0)
      pout << "\t\t\t total states using dm and quanta " << totalstatesbydm << " " << totalstatesbyquanta << endl;
  
// FIXME: rejected basis might not be necessary, after fixed this should be removed
//  vector< vector<double> > selectedwts;
    LRT::assign_matrix_by_dm(eigenMatrix, rotateMatrices[0], selectedwts, rejectedbasis, rejectedwts, transformmatrix,
                             inorderwts, wtsbyquanta, totalstatesbydm, totalstatesbyquanta, size(), dmrginp.last_site()-size());

    projected_wave_solutions[0] = wave_solutions_1st[0];
    for(int i = 1; i < Lroots; ++i) {
      if(last_site)
        projected_wave_solutions[i] = wave_solutions_1st[i];
      else
        LRT::project_onto_rejectedspace(wave_solutions_1st[i], rotateMatrices[0], dot_with_sys, projected_wave_solutions[i]);
    }
  }

  // compute 1-st order components for RPA
  if(rpa_sweep) {
    if(!rpa_sweep_2nd && !last_site) {
      LRT::multiply_h_left davidson_f(newbig, onedot);
      LRT::RPA::compute_matrix_elements(projected_wave_solutions, energies, ynorm, davidson_f,
                                        a_subspace, b_subspace, s_subspace, d_subspace, rpa_sweep_2nd, lroots);
    }
  }
  // compute 1-st order components for TDA
  else {
    if(!last_site) {
      LRT::multiply_h_left davidson_f(newbig, onedot);
      LRT::TDA::compute_matrix_elements(projected_wave_solutions, energies, davidson_f, a_subspace, s_subspace, lroots);
    }
  }

  if(mpigetrank() == 0) {
    for(int i = 1; i < Lroots; ++i) {
      DensityMatrix tracedMatrix_deriv;
      tracedMatrix_deriv.allocate(stateInfo);
      MultiplyProduct (projected_wave_solutions[i], Transpose(wave_solutions_1st[0]), tracedMatrix_deriv, 1.0);
      LRT::assign_matrix_by_dm_deriv(rotateMatrices[0], selectedwts, rejectedwts, tracedMatrix_deriv, rotateMatrices[i], false);
    }
  }

#ifndef SERIAL
  mpi::communicator world;
  broadcast(world, rotateMatrices, 0);
#endif

  for (int i = 0; i < Lroots; ++i)
    SaveRotationMatrix (newbig.leftBlock->sites, rotateMatrices[i], i);
// FOR DEBUG TEST
//for (int i = 0; i < Lroots; ++i)
  for (int i = 1; i < Lroots; ++i)
    wave_solutions_1st[i].SaveWavefunctionInfo (newbig.stateInfo, newbig.leftBlock->sites, i);

  dmrginp.rotmatrixT -> stop();
}

}
