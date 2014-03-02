/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_SOLVER_HEADER_H
#define SPIN_SOLVER_HEADER_H
#include "spinblock.h"
#include "wavefunction.h"

namespace SpinAdapted{
namespace Solver
{
  void solve_wavefunction(vector<Wavefunction>& solution, std::vector<double>& energies, SpinBlock& big, const double tol, 
			  const guessWaveTypes& guesswavetype, const bool &onedot, const bool& dot_with_sys, const bool& warmUp, double additional_noise, int currentRoot, std::vector<Wavefunction>& lowerStates);
};
}
#endif
