/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

// written by N.N. for DMRG-LRT

#ifndef LRT_SPIN_SOLVER_HEADER_H
#define LRT_SPIN_SOLVER_HEADER_H

#include "spinblock.h"
#include "wavefunction.h"

namespace SpinAdapted {

namespace LRT {

void solve_wavefunction
(vector<Wavefunction>& psix, const vector<double>& eigv, vector<double>& rnorm, SpinBlock& big, int nroots, int mroots, int kroots);

namespace TDA {

void solve_correction_equation
(vector<Wavefunction>& psix, const vector<double>& eigv, vector<double>& rnorm, const DiagonalMatrix& e, Davidson_functor& h_mult, bool useprecond, int nroots, int mroots, int kroots);

void compute_matrix_elements
(const vector<Wavefunction>& psix, Davidson_functor& h_multx, Davidson_functor& h_mult0, Matrix& h_subspace, Matrix& s_subspace, int mroots);

};

namespace RPA {

// Is anyone going to implement here the RPA solver?

};

};

};

#endif
