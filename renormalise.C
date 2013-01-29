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
#include "operatorloops.h"
#include <numeric>
#include "Operators.h"
#include "wavefunction.h"
#include "rotationmat.h"
#include "density.h"
#include "initblocks.h"
#include "guess_wavefunction.h"
#include "linear.h"
#include "davidson.h"
#include <stdlib.h>
//#include "diis.h"

#ifndef SERIAL
#include <boost/mpi.hpp>
#endif

#include "pario.h"

using namespace boost;
using namespace std;

namespace SpinAdapted{
void SpinBlock::RenormaliseFrom(vector<double> &energies, vector<double> &spins, double& error, vector<Matrix>& rotateMatrix, 
				const int keptstates, const int keptqstates, const double tol, SpinBlock& big, 
				const guessWaveTypes &guesswavetype, const double noise, const double additional_noise, const bool &onedot, SpinBlock& System, 
				SpinBlock& sysDot, SpinBlock& envDot, SpinBlock& environment, const bool& dot_with_sys,
				const bool& warmUp, int sweepiter)
{
  int nroots = dmrginp.nroots(sweepiter);
  vector<Wavefunction> wave_solutions(nroots);
  dmrginp.davidsonT -> start();
  if (dmrginp.outputlevel() > 0)
    mcheck("before davidson but after all blocks are built");

  dmrginp.solvewf -> start();
  Solver::solve_wavefunction(wave_solutions, energies, big, tol, guesswavetype, onedot, dot_with_sys, warmUp, additional_noise);

  dmrginp.solvewf -> stop();
  SpinBlock newsystem;
  SpinBlock newenvironment;
  SpinBlock newbig;
  dmrginp.postwfrearrange -> start();

  if (onedot && !dot_with_sys)
  {
    InitBlocks::InitNewSystemBlock(System, sysDot, newsystem, sysDot.size(), dmrginp.direct(), DISTRIBUTED_STORAGE, false, true);
    InitBlocks::InitBigBlock(newsystem, environment, newbig); 
    for (int i=0; i<nroots&& mpigetrank()==0; i++) 
    {
      Wavefunction tempwave = wave_solutions[i];
      GuessWave::onedot_shufflesysdot(big.get_stateInfo(), newbig.get_stateInfo(),wave_solutions[i], tempwave);  
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

  bool normalnoise = warmUp;
  if (newbig.get_rightBlock()->size() < 2)
    normalnoise = true;
  
  dmrginp.addnoise -> start();
  double twodotnoise = 0.0;
  if (dmrginp.noise_type() == RANDOM)
    twodotnoise = additional_noise;
  
  tracedMatrix.makedensitymatrix(wave_solutions, newbig, dmrginp.weights(sweepiter), noise, twodotnoise, normalnoise);
  dmrginp.addnoise -> stop();
  if (dmrginp.outputlevel() > 0)
    mcheck("after density matrix before rotation matrix");
  if (!mpigetrank())
    error = makeRotateMatrix(tracedMatrix, rotateMatrix, keptstates, keptqstates);

#ifndef SERIAL
  mpi::communicator world;
  broadcast(world, rotateMatrix, 0);
#endif

  for (int i=0; i<nroots; i++)
    SaveRotationMatrix (newbig.leftBlock->sites, rotateMatrix, i);
  for(int i=0;i<nroots;++i)
    wave_solutions[i].SaveWavefunctionInfo (newbig.stateInfo, newbig.leftBlock->sites, i);
  dmrginp.rotmatrixT -> stop();
  if (dmrginp.outputlevel() > 0)
    mcheck("after noise and calculation of density matrix");
}

double SpinBlock::makeRotateMatrix(DensityMatrix& tracedMatrix, vector<Matrix>& rotateMatrix, const int& keptstates, const int& keptqstates)
{
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
  
  return assign_matrix_by_dm(rotateMatrix, eigenMatrix, transformmatrix, inorderwts, wtsbyquanta, totalstatesbydm, 
			      totalstatesbyquanta, size(), dmrginp.last_site()-size());
}

}
