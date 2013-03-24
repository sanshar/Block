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
#include "davidson.h"

namespace SpinAdapted {

namespace LRT {

int compute_orthogonal_transform(const Matrix& s, Matrix& u, bool is_rpa_metric);
void transform_matrix(const Matrix& a, const Matrix& u, Matrix& b);
int compute_eigenvalues(const Matrix& a_subspace, const Matrix& s_subspace, DiagonalMatrix& freq, Matrix& alpha);
int compute_eigenvalues(const Matrix& a_subspace, const Matrix& b_subspace, const Matrix& s_subspace, const Matrix& d_subspace, DiagonalMatrix& freq, Matrix& alpha);

void solve_wavefunction
(vector<Wavefunction>& psix1st, vector<Wavefunction>& psix2nd, const vector<double>& eigv, vector<double>& rnorm, SpinBlock& big,
 const guessWaveTypes& guesswavetype, const bool& onedot, const bool& dot_with_sys, const double& noise,
 const bool& rpa_sweep, const bool& rpa_sweep_2nd, int nroots, int mroots, int kroots);

void compute_guess_solutions
(vector<Wavefunction>& psix, DiagonalMatrix& h_diag, Davidson_functor& h_mult, bool rpa_sweep, bool useprecond, int nroots, const bool& dot_with_sys);

namespace TDA {

void solve_correction_equation
(vector<Wavefunction>& psix, const vector<double>& eigv, vector<double>& rnorm,
 DiagonalMatrix& h_diag, Davidson_functor& h_mult, bool useprecond, int nroots, int mroots, int kroots, double noise = 0.0);

void compute_matrix_elements
(vector<Wavefunction>& psix, const vector<double>& eigv, Davidson_functor& h_mult, Matrix& a_subspace, Matrix& s_subspace, int mroots);

};

namespace RPA {

void solve_correction_equation
(vector<Wavefunction>& psix, const vector<double>& eigv, vector<double>& rnorm,
 DiagonalMatrix& h_diag, Davidson_functor& h_mult, bool useprecond, int nroots, int mroots, int kroots, double noise = 0.0);

void compute_matrix_elements
(vector<Wavefunction>& psix, const vector<double>& eigv, vector<double>& ynorm, Davidson_functor& h_mult,
 Matrix& a_subspace, Matrix& b_subspace, Matrix& s_subspace, Matrix& d_subspace, const bool& rpa_sweep_2nd, int mroots);

};

};

};

#endif
