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
#include "diis.h"

#ifndef SERIAL
#include <boost/mpi.hpp>
#endif

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
  dmrginp.davidsonT.start();
  if (dmrginp.outputlevel() != 0)
    mcheck("before davidson but after all blocks are built");

  DiagonalMatrix e;
  dmrginp.solvewf.start();
  Solver::solve_wavefunction(wave_solutions, energies, big, tol, guesswavetype, onedot, dot_with_sys, warmUp, e, additional_noise);

  dmrginp.solvewf.stop();
  SpinBlock newsystem;
  SpinBlock newenvironment;
  SpinBlock newbig;
  dmrginp.postwfrearrange.start();

  if (onedot && !dot_with_sys)
  {
    InitBlocks::InitNewSystemBlock(System, sysDot, newsystem, sysDot.size(), dmrginp.direct(), DISTRIBUTED_STORAGE, false, true);
    InitBlocks::InitBigBlock(newsystem, environment, newbig); 
    for (int i=0; i<nroots; i++) 
    {
      //wave_solutions[i].clearTensors();
      Wavefunction tempwave = wave_solutions[i];
      GuessWave::onedot_shufflesysdot(big.get_stateInfo(), newbig.get_stateInfo(),wave_solutions[i], tempwave);  
      wave_solutions[i] = tempwave;
    }
    *this = newsystem;
    envDot.clear();
    big.get_rightBlock()->clear();
    big.clear();
  }
  else
    newbig = big;
  dmrginp.postwfrearrange.stop();

  if (dmrginp.outputlevel() != 0)
    mcheck("after davidson before noise");
  //if (mpigetrank() == 0) system("free -m");
  dmrginp.davidsonT.stop();

  dmrginp.rotmatrixT.start();
  DensityMatrix tracedMatrix;
  tracedMatrix.allocate(stateInfo);

  //pout <<"\t\t\t System state Info\n"<< stateInfo<<endl;
  // do post-solution density matrix modifications
  if (dmrginp.outputlevel() != 0)
    pout <<"\t\t\t Total left, right and keptstates quanta: "<<newbig.get_stateInfo().leftStateInfo->totalStates<<" "<<newbig.get_stateInfo().rightStateInfo->totalStates<<" "<<keptstates<<endl;
  int leftQuantaUsed = 0;
  for (int i=0; i<wave_solutions[0].nrows(); i++) {
    int maxUsed = 0;
    for (int j=0; j<wave_solutions[0].ncols(); j++)
      if( wave_solutions[0].allowed(i, j))
	if (maxUsed < min(wave_solutions[0](i,j).Nrows(), wave_solutions[0](i,j).Ncols()))
	  maxUsed = min(wave_solutions[0](i,j).Nrows(), wave_solutions[0](i,j).Ncols());
    leftQuantaUsed += maxUsed;
  }
  if (dmrginp.outputlevel() != 0)
    pout <<"\t\t\t Left Quanta used in state 0: "<<leftQuantaUsed<<endl;
  bool normalnoise = warmUp;
  if (newbig.get_rightBlock()->size() < 2)
    normalnoise = true;
  
  dmrginp.addnoise.start();
  double twodotnoise = 0.0;
  if (dmrginp.noise_type() == RANDOM)
    twodotnoise = additional_noise;
  
  tracedMatrix.makedensitymatrix(wave_solutions, newbig, dmrginp.weights(sweepiter), noise, twodotnoise, normalnoise, sweepiter);
  dmrginp.addnoise.stop();
  if (dmrginp.outputlevel() != 0)
    mcheck("after density matrix before rotation matrix");
  if (!mpigetrank())
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

      if (dmrginp.outputlevel() != 0)
	pout << "\t\t\t total states using dm and quanta " << totalstatesbydm << " " << totalstatesbyquanta << endl;

      /*
      for (int i=0; i<totalstatesbydm; i++)
	cout << newbig.leftBlock->get_stateInfo().quanta[inorderwts[i].first]<<" "<<eigenMatrix[inorderwts[i].first].element(inorderwts[i].second, inorderwts[i].second)<<endl;
      */

      error = assign_matrix_by_dm(rotateMatrix, eigenMatrix, transformmatrix, inorderwts, wtsbyquanta, totalstatesbydm, 
				  totalstatesbyquanta, newbig.get_leftBlock()->size(), newbig.get_rightBlock()->size());
    }
#ifndef SERIAL
  mpi::communicator world;
  broadcast(world, rotateMatrix, 0);
#endif

  SaveRotationMatrix (newbig.leftBlock->sites, rotateMatrix);
  for(int i=0;i<nroots;++i)
    wave_solutions[i].SaveWavefunctionInfo (newbig.stateInfo, newbig.leftBlock->sites, i);
  dmrginp.rotmatrixT.stop();
  if (dmrginp.outputlevel() != 0)
    mcheck("after noise and calulation of density matrix");
}
}
