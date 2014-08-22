#ifndef SPIN_GUESS_WAVEFUNCTION_HEADER
#define SPIN_GUESS_WAVEFUNCTION_HEADER
#include "spinblock.h"
#include "StateInfo.h"
#include "ObjectMatrix.h"
#include "operatorfunctions.h"
#include "rotationmat.h"
#include "wavefunction.h"
#include <vector>

namespace SpinAdapted{

namespace GuessWave
{
  void TransformLeftBlock(Wavefunction& oldwavefunction, const StateInfo& newstateinfo, const std::vector<Matrix>& RotationMatrix, Wavefunction& tempoldWave);
  void TransformRightBlock(const Wavefunction& tempnewWave, const StateInfo& tempoldStateInfo, const std::vector<Matrix>& RotationMatrix, Wavefunction& trial);
  void guess_wavefunctions(std::vector<Wavefunction>& solution, DiagonalMatrix& e, const SpinBlock &big, 
			   const guessWaveTypes &guesswavetype, const bool &onedot,  const bool& transpose_guess_wave, double additional_noise=0.0, int currentState=0);
  void guess_wavefunctions(Wavefunction& solution, DiagonalMatrix& e, const SpinBlock &big,
			   const guessWaveTypes &guesswavetype, const bool &onedot, const int &state, const bool& transpose_guess_wave, double additional_noise=0.0);

  void guess_wavefunctions(Wavefunction& solution, DiagonalMatrix& e, const SpinBlock &big,
			   const guessWaveTypes &guesswavetype, const bool &onedot, const int &state, const bool& transpose_guess_wave, double additional_noise, bool ket);

  //onedot transpose wave guess
  void transpose_previous_wavefunction(Wavefunction& trial, const StateInfo& stateInfo, const std::vector<int>& rightsites, const std::vector<int> &dotsites, const int state, const bool &onedot, const bool& transpose_guess_wave);
  void transpose_previous_wavefunction(Wavefunction& trial, const SpinBlock &big, const int state, const bool &onedot, const bool& transpose_guess_wave, bool ket=true);
  void onedot_transpose_wavefunction(const StateInfo& guessstateinfo, const StateInfo& transposestateinfo,
                                      const Wavefunction& guesswf, Wavefunction& transposewf);
  void onedot_threeindex_to_twoindex_wavefunction(const StateInfo& twostateinfo, const ObjectMatrix3D< std::vector<Matrix> >& 
						  threewavefunction, Wavefunction& twowavefunction, const StateInfo& guessstateinfo);
  void onedot_twoindex_to_threeindex_wavefunction(const StateInfo& stateinfo, const Wavefunction& twowavefunction, 
						  ObjectMatrix3D< std::vector<Matrix> > & threewavefunction);

  void onedot_shufflesysdot(const StateInfo& guessstateinfo, const StateInfo& transposestateinfo,
			    const Wavefunction& guesswf, Wavefunction& transposewf);

  void onedot_twoindex_to_threeindex_shufflesysdot(const StateInfo& stateinfo, const Wavefunction& twowavefunction, 
						   ObjectMatrix3D< vector<Matrix> >& threewavefunction);
  void transform_previous_wavefunction(Wavefunction& trial, const SpinBlock &big, const int state, const bool &onedt, 
				       const bool& transpose_guess_wave);

  void transform_previous_wavefunction(Wavefunction& trial, const SpinBlock &big, const int state, const bool &onedt, 
				       const bool& transpose_guess_wave, bool ket);
  void transform_previous_wavefunction(Wavefunction& trial, const StateInfo& stateInfo, const std::vector<int> &leftsites, const std::vector<int> &rightsites, const int state, const bool &onedot, const bool& transpose_guess_wave);

  //ondeot transform wave guess                                                                                  
  void onedot_transform_wavefunction(const StateInfo& oldstateinfo, const StateInfo& newstateinfo, const Wavefunction& oldwavefunction,
				     const std::vector<Matrix>& inverseLeftRotationMatrix,
                                     const std::vector<Matrix>& rightRotationMatrix, Wavefunction& newwavefunction, 
				     const bool& transpose_guess_wave, bool ket = true);
  void basic_guess_wavefunction(DiagonalMatrix& e, Wavefunction& trial, const StateInfo *stateinfo, const int state);
  void transform_previous_twodot_to_onedot_wavefunction(Wavefunction& trial, const SpinBlock &big, const int state);

}
}

#endif
