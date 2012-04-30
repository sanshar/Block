#include "diis.h"
#include "global.h"
#include "solver.h"
#include "initblocks.h"
#include "MatrixBLAS.h"
#include "density.h"
#include "rotationmat.h"
#include "davidson.h"
#include "linear.h"
#include "guess_wavefunction.h"
#include <include/sortutils.h>

#ifndef SERIAL
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif

using namespace boost;
using namespace std;



double SpinAdapted::DIIS::do_one(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize)
{

  updateWavefunction(sweepParams);
  //if(!(currentIndex == 0 && buildup))
  updateHamiltonian(sweepParams, false);

  double error = updateErrors(sweepParams, true);
  updateB(sweepParams);

  ++sweepParams.set_sweep_iter();
  ++sweepParams.set_sweep_iter();
  currentIndex = (currentIndex == keepStates-1) ? 0 : currentIndex+1; 
  if (buildup && currentIndex == 0)
    buildup = false;



  return error;
}


void SpinAdapted::DIIS::initialize(SweepParams& sweepParams)
{
  buildup = true;
  b.ReSize(keepStates+1);
  x.ReSize(keepStates+1);

  B.ReSize(keepStates+1, keepStates+1);
  B = 0.0;

  for (int i=0; i<keepStates; i++) {
    b(i+1) = 0.0;
    B(keepStates+1, i+1) = -1.0;
    B(i+1, keepStates+1) = -1.0;
    
    for (int j=0; j<keepStates; j++) {
      B(i+1,j+1) = double(rand())/RAND_MAX*1.0e6;
      B(j+1,i+1) = B(i+1,j+1);
    }
  }

  b(keepStates+1) = -1.0;
  currentIndex = 0;
  cout << B<<endl;
  //cout << B<<endl;
}

void SpinAdapted::DIIS::updateB(SweepParams& sweepParams)
{
  Wavefunction w1, w2;

  for (int j=0; j<keepStates; j++) {
    B(currentIndex+1, j+1) = 0.0;
    B(j+1, currentIndex+1) = 0.0;
  }
  
  for (int i=0; i< sweepParams.get_n_iters() ; i++) {
    LoadDIISError(w1, currentIndex, i);
    for (int j=0; j<keepStates; j++) { 
      if(buildup && j> currentIndex)
	continue;
      LoadDIISError(w2, j, i);
      B(currentIndex+1, j+1) += DotProduct(w1, w2);
      B(j+1, currentIndex+1)  = B(currentIndex+1, j+1);
    }
  }

  /*
  cout <<"  B= "<<endl<< B<<endl;
  Matrix U, V;
  DiagonalMatrix d;
  svd(B, d, U, V);
  */
}
 
void SpinAdapted::DIIS::updateWavefunction(SweepParams &sweepParams)
{
  xsolve_AxeqB(B, b, x);
  cout << B<<endl;
  cout << x<<endl;
  //remove these lines
  int ind = currentIndex == 0 ?  keepStates-1 : currentIndex-1;
  x = 0.0; x(ind + 1) = 1.0;
  //cout << x<<endl;
  //till here
  
  Wavefunction w, w2;
  
  std::vector<int> sites(1,0); sites.push_back(1);
  
  Wavefunction prevW;
    
  for (int i=0; i< sweepParams.get_n_iters() ; i++) {
    if(buildup && currentIndex == 0) {
      return;
      //StateInfo newinfo;
      //w.LoadWavefunctionInfo(newinfo, sites, 0);
      //sites.push_back(i+2);
      //SaveInterpolatedWavefunction(w, i); 
      //continue;
    }
    
    for (int j=0; j<keepStates; j++) {
      if(buildup && j>=currentIndex)
	break;
      LoadDIISWavefunction(w2, j, i);
      
      if ( j == 0) {
	w = w2;
	Scale( x(j+1), w);
      }
      else 
	ScaleAdd(x(j+1), w2, w);
      
    }
    Normalise(w);
    //cout << w <<endl;
    SaveInterpolatedWavefunction(w, i); 
  }

}


void SpinAdapted::DIIS::SaveDIISError(const Wavefunction& sigma, int currentIndex, int iter)
{
  char file [30];
  sprintf (file, "%s%s%d%s%d%s", dmrginp.save_prefix().c_str(), "Error-", currentIndex, ".", iter , ".tmp");
  //pout << "Save Error "<<file<<endl;
  if (mpigetrank() == 0)
    {
      std::ofstream ofs(file, std::ios::binary);
      boost::archive::binary_oarchive save_wave(ofs);
      save_wave << sigma;
      ofs.close();
    }
}

void SpinAdapted::DIIS::LoadDIISError(Wavefunction& sigma, int currentIndex, int iter)
{
  char file [30];
  sprintf (file, "%s%s%d%s%d%s", dmrginp.save_prefix().c_str(), "Error-", currentIndex, ".", iter , ".tmp");
  //pout << "Load Error "<<file<<endl;
  if (mpigetrank() == 0)
    {
      std::ifstream ifs(file, std::ios::binary);
      boost::archive::binary_iarchive save_wave(ifs);
      save_wave >> sigma;
      ifs.close();
    }
}

void SpinAdapted::DIIS::SaveDIISWavefunction(const Wavefunction& sigma, int currentIndex, int iter)
{
  char file [30];
  sprintf (file, "%s%s%d%s%d%s", dmrginp.save_prefix().c_str(), "DIISWave-", currentIndex, ".", iter , ".tmp");
  pout << "Save Wavefunction "<<file<<endl;
  if (mpigetrank() == 0)
    {
      std::ofstream ofs(file, std::ios::binary);
      boost::archive::binary_oarchive save_wave(ofs);
      save_wave << sigma;
      ofs.close();
    }
}

void SpinAdapted::DIIS::LoadDIISWavefunction(Wavefunction& sigma, int currentIndex, int iter)
{
  char file [30];
  sprintf (file, "%s%s%d%s%d%s", dmrginp.save_prefix().c_str(), "DIISWave-", currentIndex, ".", iter , ".tmp");
  //pout << "Load Wavefunction "<<file<<endl;
  if (mpigetrank() == 0)
    {
      std::ifstream ifs(file, std::ios::binary);
      boost::archive::binary_iarchive save_wave(ifs);
      save_wave >> sigma;
      ifs.close();
    }
}

void SpinAdapted::DIIS::SaveInterpolatedWavefunction(const Wavefunction& sigma, int iter)
{
  char file [30];
  sprintf (file, "%s%s%d%s", dmrginp.save_prefix().c_str(), "InterpolatedWave-", iter , ".tmp");
  pout << "Saving wavefunction "<<file<<endl;
  if (mpigetrank() == 0)
    {
      std::ofstream ofs(file, std::ios::binary);
      boost::archive::binary_oarchive save_wave(ofs);
      save_wave << sigma;
      ofs.close();
    }
}

void SpinAdapted::DIIS::LoadInterpolatedWavefunction(Wavefunction& sigma, int iter)
{
  char file [30];
  sprintf (file, "%s%s%d%s", dmrginp.save_prefix().c_str(), "InterpolatedWave-", iter , ".tmp");
  pout << "Loading wavefunction "<<file<<endl;
  if (mpigetrank() == 0)
    {
      std::ifstream ifs(file, std::ios::binary);
      boost::archive::binary_iarchive save_wave(ifs);
      save_wave >> sigma;
      ifs.close();
    }
}

