#ifndef LRT_SPIN_TRANSFORM_GAUGE_HEADER
#define LRT_SPIN_TRANSFORM_GAUGE_HEADER
#include "spinblock.h"
#include "StateInfo.h"
#include "ObjectMatrix.h"
#include "operatorfunctions.h"
#include "rotationmat.h"
#include "wavefunction.h"
#include <vector>

namespace SpinAdapted {

namespace GuessWave {

namespace LRT {

void transform_gauge
(std::vector<Wavefunction>& Wavefnc, const int& nroots, const SpinBlock& big,
 const guessWaveTypes& guesswavetype, const bool& onedot, const bool& transpose_guess_wave);

void transform_previous_wavefunction_deriv
(std::vector<Wavefunction>& Wavefnc, const int& nroots, const SpinBlock& big, const bool& onedot, const bool& transpose_guess_wave);

void transpose_previous_wavefunction_deriv
(std::vector<Wavefunction>& Wavefnc, const int& nroots, const SpinBlock& big, const bool& onedot, const bool& transpose_guess_wave);

void rotate_previous_wavefunction
(const SpinBlock& big, const Matrix& alpha, const int& nroots,
 const guessWaveTypes& guesswavetype, const bool& onedot, const bool& transpose_guess_wave);

};

};

};

#endif
