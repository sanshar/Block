/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

// Written by N.N. for DMRG-LRT

#include "guess_wavefunction.h"
#include "modules/lrt/lrt_transform_gauge.h"
#ifndef SERIAL
#include <boost/mpi.hpp>
#endif
#include "pario.h"

namespace SpinAdapted {

//
// GAUGE-TRANSFORMATION
//

void GuessWave::LRT::transform_gauge
(std::vector<Wavefunction>& Wavefnc, const int& nroots, const SpinBlock &big,
 const guessWaveTypes& guesswavetype, const bool &onedot, const bool &transpose_guess_wave)
{
#ifndef SERIAL
  mpi::communicator world;
#endif
  assert(nroots <= Wavefnc.size());
  for(int i = 0; i < nroots; ++i) {
    Wavefnc[i].initialise(dmrginp.effective_molecule_quantum(), &big, onedot);
  }

  if (!mpigetrank()) {
    switch(guesswavetype)
    {
    case TRANSFORM: 
      transform_previous_wavefunction_deriv(Wavefnc, nroots, big, onedot, transpose_guess_wave);
      break;
    case TRANSPOSE: 
      transpose_previous_wavefunction_deriv(Wavefnc, nroots, big, onedot, transpose_guess_wave);
      break;
    defalt:
      pout << "\t\t\t Invalid guessWaveTypes for gauge-transformation was specified" << endl;
      abort();
    }

//  // DEBUG: begin check gauge-condition
//  std::vector<Wavefunction> transWavefnc(nroots);
//  for(int i = 0; i < nroots; ++i) {
//    transWavefnc[i].initialise(dmrginp.effective_molecule_quantum(), &big, onedot);
//  }
//  transpose_previous_wavefunction_deriv(transWavefnc, nroots, big, onedot, transpose_guess_wave);

//  for(int i = 0; i < nroots; ++i) {
//    Wavefunction subWavefnc(Wavefnc[i]);
//    ScaleAdd(-1.0, transWavefnc[i], subWavefnc);
//    double subnorm = DotProduct(subWavefnc, subWavefnc);
//    pout << "\t\t\t DEBUG @ GuessWave::LRT::transform_gauge: checking gauge of Wavefnc[" << i << "]: diff = " << subnorm << endl;
//  }
//  // DEBUG: end   check gauge-condition

    double norm = DotProduct(Wavefnc[0], Wavefnc[0]);
    if(dmrginp.outputlevel() > 0)
      pout << "\t\t\t norm of 0-th wavefunction = " << fixed << setprecision(24) << norm << endl;
  }
}

// transform gauge of 1-st order wavefunction
void GuessWave::LRT::transform_previous_wavefunction_deriv
(std::vector<Wavefunction>& Wavefnc, const int& nroots, const SpinBlock &big, const bool &onedot, const bool& transpose_guess_wave)
{
  assert(nroots <= Wavefnc.size());

  if (dmrginp.outputlevel() > 0) 
    pout << "\t\t\t Transforming gauge of previous wavefunction " << endl;

  if(!onedot) {
    pout << "\t\t\t Gauge-transformation has implemented only for onedot wavefunction" << endl;
    abort();
  }

  // 0-th order wavefuncton & rotation matrices
  StateInfo    oldStateInfo;
  Wavefunction oldWavefnc;
  std::vector<Matrix> leftRotationMatrix;
  std::vector<Matrix> rightRotationMatrix;

  if (transpose_guess_wave) {
    oldWavefnc.LoadWavefunctionInfo (oldStateInfo, big.get_leftBlock()->get_leftBlock()->get_sites(), 0);
    LoadRotationMatrix (big.get_leftBlock()->get_leftBlock()->get_sites(), leftRotationMatrix, 0);
  }
  else{
    oldWavefnc.LoadWavefunctionInfo (oldStateInfo, big.get_leftBlock()->get_sites(), 0);
    LoadRotationMatrix (big.get_leftBlock()->get_sites(), leftRotationMatrix, 0);
  }

  // transpose leftRotationMatrix
  for (int q = 0; q < leftRotationMatrix.size (); ++q) {
    if (leftRotationMatrix [q].Nrows () > 0) {
      leftRotationMatrix[q] = leftRotationMatrix[q].t();
    }
  }

  std::vector<int> rotsites;
  if (transpose_guess_wave) {
    rotsites = big.get_rightBlock()->get_sites();
    rotsites.insert(rotsites.end(), big.get_leftBlock()->get_rightBlock()->get_sites().begin(),
                                    big.get_leftBlock()->get_rightBlock()->get_sites().end());
    sort(rotsites.begin(), rotsites.end());
  }
  else{
    rotsites = big.get_rightBlock()->get_sites();
  }

  // should be load 0-th wavefunction rather than transform
  // since follows may introduce round-off-error which breaks normality of 0-th wavefunction

  // old right rotation matrix
  LoadRotationMatrix(rotsites, rightRotationMatrix, 0);
  // transform 0-th gauge
  onedot_transform_wavefunction(oldStateInfo, big.get_stateInfo(), oldWavefnc,
                                leftRotationMatrix, rightRotationMatrix, Wavefnc[0], transpose_guess_wave);

  for (int i = 1; i < nroots; ++i) {
    // 1-st order wavefuncton & rotation matrices
    StateInfo    oldStateInfo_deriv;
    Wavefunction oldWavefnc_deriv;
    std::vector<Matrix> rightRotationMatrix_deriv; // had to be rotated by Ritz vector and saved when making newEnvironment

    if (transpose_guess_wave) {
      // old wavefunction
      oldWavefnc_deriv.LoadWavefunctionInfo (oldStateInfo_deriv, big.get_leftBlock()->get_leftBlock()->get_sites(), i);
      // old right rotation matrix
      LoadRotationMatrix(rotsites, rightRotationMatrix_deriv, i);
    }
    else{
      // old wavefunction
      oldWavefnc_deriv.LoadWavefunctionInfo (oldStateInfo_deriv, big.get_leftBlock()->get_sites(), i);
      // old right rotation matrix
      LoadRotationMatrix(rotsites, rightRotationMatrix_deriv, i);
    }

    Wavefunction tmpWavefnc = Wavefnc[i];
    // L(0)^(t) * C(I) * R(0)
    onedot_transform_wavefunction(oldStateInfo_deriv, big.get_stateInfo(), oldWavefnc_deriv,
                                  leftRotationMatrix, rightRotationMatrix, tmpWavefnc, transpose_guess_wave);
    // L(0)^(t) * C(0) * R(I)
    onedot_transform_wavefunction(oldStateInfo, big.get_stateInfo(), oldWavefnc,
                                  leftRotationMatrix, rightRotationMatrix_deriv, Wavefnc[i], transpose_guess_wave);
    Wavefnc[i] += tmpWavefnc;
  }
}

// transpose 1-st order wavefunction ( block iter == 0 )
void GuessWave::LRT::transpose_previous_wavefunction_deriv
(std::vector<Wavefunction>& Wavefnc, const int& nroots, const SpinBlock &big, const bool &onedot, const bool& transpose_guess_wave)
{
  assert(nroots <= Wavefnc.size());
  for(int i = 0; i < nroots; ++i) {
    transpose_previous_wavefunction(Wavefnc[i], big, i, onedot, transpose_guess_wave);
  }
}

void GuessWave::LRT::rotate_previous_wavefunction
(const SpinBlock& big, const Matrix& alpha, const int& nroots,
 const guessWaveTypes& guesswavetype, const bool& onedot, const bool& transpose_guess_wave)
{
  const int mroots = alpha.Nrows()+1;

  if(!onedot) {
    pout << "\t\t\t GuessWave::LRT::rotate_previous_wavefunction  could not perform with twodot" << endl;
    abort();
  }

  if(guesswavetype == TRANSPOSE) {
    std::vector<StateInfo>    oldStateInfo(mroots);
    std::vector<Wavefunction> oldWavefnc(mroots);

    std::vector<int> wfsites = big.get_rightBlock()->get_sites();
    wfsites.insert(wfsites.end(), big.get_leftBlock()->get_rightBlock()->get_sites().begin(), big.get_leftBlock()->get_rightBlock()->get_sites().end());
    sort(wfsites.begin(), wfsites.end());
    for(int i = 1; i < mroots; ++i) {
      oldWavefnc[i].LoadWavefunctionInfo(oldStateInfo[i], wfsites, i);
    }
    std::vector<Wavefunction> oldWavefncSave;
    oldWavefncSave = oldWavefnc;
    for(int i = 1; i < nroots; ++i) {
      Scale(alpha(i, i), oldWavefnc[i]);
    }
    for(int i = 1; i < nroots; ++i) {
      for(int j = 1; j < mroots; ++j) {
        if(i != j) {
          ScaleAdd(alpha(j, i), oldWavefncSave[j], oldWavefnc[i]);
        }
      }
    }
    for(int i = 1; i < nroots; ++i) {
      oldWavefnc[i].SaveWavefunctionInfo(oldStateInfo[i], wfsites, i);
    }
  }
  else if(guesswavetype == TRANSFORM) {
    std::vector<std::vector<Matrix> > rightRotationMatrices(mroots);

    std::vector<int> rotsites;
    if (transpose_guess_wave) {
      rotsites = big.get_rightBlock()->get_sites();
      rotsites.insert(rotsites.end(), big.get_leftBlock()->get_rightBlock()->get_sites().begin(), big.get_leftBlock()->get_rightBlock()->get_sites().end());
      sort(rotsites.begin(), rotsites.end());
      for(int i = 1; i < mroots; ++i) {
        LoadRotationMatrix(rotsites, rightRotationMatrices[i], i);
      }
    }
    else {
      rotsites = big.get_rightBlock()->get_sites();
      for(int i = 1; i < mroots; ++i) {
        LoadRotationMatrix(rotsites, rightRotationMatrices[i], i);
      }
    }
    std::vector<std::vector<Matrix> > rightRotationMatricesSave;
    rightRotationMatricesSave = rightRotationMatrices;
    for(int i = 1; i < nroots; ++i) {
      for(int q = 0; q < rightRotationMatrices[i].size(); ++q) {
        MatrixScale(alpha(i, i), rightRotationMatrices[i][q]);
      }
    }
    for(int i = 1; i < nroots; ++i) {
      for(int j = 1; j < mroots; ++j) {
        if(i != j) {
          for(int q = 0; q < rightRotationMatrices[i].size(); ++q) {
            MatrixScaleAdd(alpha(j, i), rightRotationMatricesSave[j][q], rightRotationMatrices[i][q]);
          }
        }
      }
    }
    for(int i = 1; i < nroots; ++i) {
      SaveRotationMatrix(rotsites, rightRotationMatrices[i], i);
    }
  }
}

}
