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
(const vector<double> &energies, vector<double>& rnorm, 
 vector< vector<Matrix> >& rotateMatrices, int nroots, int mroots, int kroots,
 Matrix& h_subspace, Matrix& s_subspace, const int keptstates, const int keptqstates,
 SpinBlock& big, const guessWaveTypes &guesswavetype, const bool &onedot, const bool &last_site, SpinBlock& System, 
 SpinBlock& sysDot, SpinBlock& envDot, SpinBlock& environment, const bool& dot_with_sys, int sweepiter)
{
  const int lroots = mroots + nroots - kroots;
  vector<Wavefunction> wave_solutions;
  dmrginp.davidsonT -> start();
  if (dmrginp.outputlevel() > 0)
    mcheck("before davidson but after all blocks are built");

  dmrginp.solvewf -> start();
  LRT::solve_wavefunction(wave_solutions, energies, rnorm, big, guesswavetype, onedot, dot_with_sys, nroots, mroots, kroots);

  // specialize to use dot with environment at last site (necessary for parallel run)
  if(last_site) {
    LRT::multiply_h_left davidson_f(big, onedot);
    LRT::TDA::compute_matrix_elements(wave_solutions, davidson_f, h_subspace, s_subspace, lroots);
  }

  dmrginp.solvewf -> stop();
  SpinBlock newsystem;
  SpinBlock newenvironment;
  SpinBlock newbig;
  dmrginp.postwfrearrange -> start();

  // maybe better to turn 'dot with env' off for computing (1 - L(0)L(0)')C(I) later
  if (onedot && !dot_with_sys) {
    InitBlocks::InitNewSystemBlock(System, sysDot, newsystem, sysDot.size(), dmrginp.direct(), DISTRIBUTED_STORAGE, false, true, lroots);
    InitBlocks::InitBigBlock(newsystem, environment, newbig); 
    for (int i = 0; i < lroots && mpigetrank() == 0; i++) {
      Wavefunction tempwave = wave_solutions[i];
      GuessWave::onedot_shufflesysdot(big.get_stateInfo(), newbig.get_stateInfo(), wave_solutions[i], tempwave);  
      wave_solutions[i] = tempwave;
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

  rotateMatrices.resize(lroots);
  vector<Wavefunction> projected_wave_solutions(lroots);

  vector< vector<double> > selectedwts;
  vector< vector<double> > rejectedwts;
  vector< Matrix > rejectedbasis;

  if(mpigetrank() == 0) {
    // FIXME: re-computing 0-th rotation matrix is wasteful, but computing 1-st rotation matrices needs its eigenvalues
    MultiplyProduct (wave_solutions[0], Transpose(wave_solutions[0]), tracedMatrix, 1.0);

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

    projected_wave_solutions[0] = wave_solutions[0];
    for(int i = 1; i < lroots; ++i) {
      if(last_site)
        projected_wave_solutions[i] = wave_solutions[i];
      else
        LRT::project_onto_rejectedspace(wave_solutions[i], rotateMatrices[0], dot_with_sys, projected_wave_solutions[i]);
    }
  }

  // might be used the same subroutine for both TDA and RPA ?
  if(!last_site) {
    LRT::multiply_h_left davidson_f(newbig, onedot);
    LRT::TDA::compute_matrix_elements(projected_wave_solutions, davidson_f, h_subspace, s_subspace, lroots);
  }

  if(mpigetrank() == 0) {
    for(int i = 1; i < lroots; ++i) {
      DensityMatrix tracedMatrix_deriv;
      tracedMatrix_deriv.allocate(stateInfo);
      // FIXME: one of them doesn't contribute
      // (maybe second one, but not sure, since it depends on storage structure of BLOCK code)
      MultiplyProduct (projected_wave_solutions[i], Transpose(wave_solutions[0]), tracedMatrix_deriv, 1.0);
//pout << "DEBUG @ SpinBlock::RenormaliseFrom_lrt: check point dm_deriv - 3" << endl;
//      MultiplyProduct (wave_solutions[0], Transpose(wave_solutions[i]), tracedMatrix_deriv, 1.0);

//    LRT::assign_matrix_by_dm_deriv(rotateMatrices[0], selectedwts, rejectedbasis, rejectedwts, tracedMatrix_deriv, rotateMatrices[i]);
//    LRT::project_onto_rejectedspace(wave_solutions[i], rejectedbasis, dot_with_sys, projected_wave_solutions[i]);
//    LRT::assign_matrix_by_dm_deriv(rotateMatrices[0], selectedwts, tracedMatrix_deriv, rotateMatrices[i], true);

//    LRT::assign_matrix_by_dm_deriv(rotateMatrices[0], selectedwts, tracedMatrix_deriv, rotateMatrices[i], false);
      LRT::assign_matrix_by_dm_deriv(rotateMatrices[0], selectedwts, rejectedwts, tracedMatrix_deriv, rotateMatrices[i], false);

//    LRT::project_onto_rejectedspace(wave_solutions[i], rotateMatrices[0], dot_with_sys, projected_wave_solutions[i]);
    }
  }

#ifndef SERIAL
  mpi::communicator world;
  broadcast(world, rotateMatrices, 0);
#endif

  for (int i = 0; i < lroots; ++i)
    SaveRotationMatrix (newbig.leftBlock->sites, rotateMatrices[i], i);
  for (int i = 0; i < lroots; ++i)
//for (int i = 1; i < lroots; ++i) // keep 0-th wavefunction
    wave_solutions[i].SaveWavefunctionInfo (newbig.stateInfo, newbig.leftBlock->sites, i);
  dmrginp.rotmatrixT -> stop();
}

}
