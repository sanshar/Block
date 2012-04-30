#ifndef SPIN_SOLVER_HEADER_H
#define SPIN_SOLVER_HEADER_H
#include "spinblock.h"
#include "wavefunction.h"

namespace SpinAdapted{
namespace Solver
{
  void solve_wavefunction(vector<Wavefunction>& solution, std::vector<double>& energies, SpinBlock& big, const double tol, 
			  const guessWaveTypes& guesswavetype, const bool &onedot, const bool& dot_with_sys, const bool& warmUp, DiagonalMatrix& e, double additional_noise=0.0);
  //double compute_spin(const Wavefunction &wave, const SpinBlock& big);
  //void compute_spin_wave(Wavefunction &spin_wave, const Wavefunction &wave, const SpinBlock& big);
};
}
#endif
