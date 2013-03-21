/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

// Written by N.N. for DMRG-LRT

#include "linear.h"
#include "davidson.h"
#include "blas_calls.h"
#ifndef SERIAL
#include <boost/mpi.hpp>
#endif
#include "pario.h"

#include "modules/lrt/lrt_solver.h"
#include "modules/lrt/lrt_davidson.h"
#include "modules/lrt/lrt_transform_gauge.h"

//
// subroutines for generalized eigenvalue problem
//

//
// compute orthogonal transformation from metric / returns # of orthogonal basis
//
int SpinAdapted::LRT::compute_orthogonal_transform(const Matrix& s, Matrix& u, bool is_rpa_metric)
{
  int nrows = s.Nrows();

  Matrix& s_ref = const_cast<Matrix&>(s);
  DiagonalMatrix d;
  Matrix p;
  diagonalise(s_ref, d, p);

  // ignore small and/or negative eigenvalues
  int nignore = 0;
  if(is_rpa_metric) {
    for(; nignore < nrows; nignore += 2)
      if(d.element(nignore, nignore) >= 1.0e-16) break;
  }
  else {
    for(; nignore < nrows; ++nignore)
      if(d.element(nignore, nignore) >= 1.0e-16) break;
  }

  int ncols = nrows - nignore;
  u.ReSize(nrows, ncols);

  for(int i = 0; i < ncols; ++i) {
    int ix = nignore + i;
    double di = 1.0/sqrt(d.element(ix, ix));
    for(int j = 0; j < nrows; ++j) {
      u.element(j, i) = di * p.element(j, ix);
    }
  }

  return ncols;
}

//
// Transforming matrix with orthogonal basis
//
void SpinAdapted::LRT::transform_matrix(const Matrix& a, const Matrix& u, Matrix& b)
{
  int m = u.Nrows();
  int n = u.Ncols();
  Matrix btmp(m, n); btmp = 0.0;
  MatrixMultiply(a, 'n', u, 'n', btmp, 1.0);
  b.ReSize(n, n); b = 0.0;
  MatrixMultiply(u, 't', btmp, 'n', b, 1.0);
}

//
// solve eigenvalue problem for TDA:
//  [ A ][X] = h[ S ][X]
//
int SpinAdapted::LRT::compute_eigenvalues
(const Matrix& a_subspace, const Matrix& s_subspace, DiagonalMatrix& freq, Matrix& alpha)
{
  int nrows = a_subspace.Nrows();

  Matrix u_matrix;
  int ncols = compute_orthogonal_transform(s_subspace, u_matrix, false);

  Matrix a_transformed;
  transform_matrix(a_subspace, u_matrix, a_transformed);

  DiagonalMatrix eigenvalues;
  Matrix eigenvectors;
  diagonalise(a_transformed, eigenvalues, eigenvectors);

  // ignore small and/or negative eigenvalues
  int nignore = 0;
  for(; nignore < nrows; ++nignore)
    if(eigenvalues.element(nignore, nignore) >= 1.0e-8) break;

  if(nignore != 0)
    pout << "\t\t\t Warning: found and ignored zero frequency modes in TDA" << endl;

  ncols -= nignore;

  freq.ReSize(ncols);
  Matrix eigvectmp(eigenvectors.Nrows(), ncols);

  for(int i = 0; i < ncols; ++i) {
    int ix = nignore + i;
    freq.element(i, i) = eigenvalues.element(ix, ix);
    for(int j = 0; j < eigenvectors.Nrows(); ++j) {
      eigvectmp.element(j, i) = eigenvectors.element(j, ix);
    }
  }

  alpha.ReSize(nrows, ncols); alpha = 0.0;
  MatrixMultiply(u_matrix, 'n', eigvectmp, 'n', alpha, 1.0);

  return ncols;
}

//
// solve eigenvalue problem for RPA:
//  [ A  B ][X] = h [ S  D ][X]
//  [ B  A ][Y]     [-D -S ][Y]
//
int SpinAdapted::LRT::compute_eigenvalues
(const Matrix& a_subspace, const Matrix& b_subspace,
 const Matrix& s_subspace, const Matrix& d_subspace, DiagonalMatrix& freq, Matrix& alpha)
{
  int nrows = a_subspace.Nrows();
  int ncols = a_subspace.Ncols();

  // construct [ A  B ]
  //           [ B  A ]
  Matrix h_matrix(2*nrows, 2*ncols);
  for(int i = 0; i < nrows; ++i) {
    for(int j = 0; j < ncols; ++j) {
      h_matrix.element(i,       j      ) = a_subspace.element(i, j);
      h_matrix.element(i,       j+ncols) = b_subspace.element(i, j);
      h_matrix.element(i+nrows, j      ) = b_subspace.element(i, j);
      h_matrix.element(i+nrows, j+ncols) = a_subspace.element(i, j);
    }
  }
  // construct [ S  D ]
  //           [-D -S ]
  Matrix s_matrix(2*nrows, 2*ncols);
  for(int i = 0; i < nrows; ++i) {
    for(int j = 0; j < ncols; ++j) {
      s_matrix.element(i,       j      ) = s_subspace.element(i, j);
      s_matrix.element(i,       j+ncols) = d_subspace.element(i, j);
      s_matrix.element(i+nrows, j      ) =-d_subspace.element(i, j);
      s_matrix.element(i+nrows, j+ncols) =-s_subspace.element(i, j);
    }
  }

  // to solve non-hermitian problem, h_matrix is used as a metric
  //  1/ [ A  B ][X'] = [ S  D ][X']
  //  /h [ B  A ][Y']   [-D -S ][Y']
  Matrix u_matrix;
  ncols = compute_orthogonal_transform(h_matrix, u_matrix, true);

  Matrix s_transformed;
  transform_matrix(s_matrix, u_matrix, s_transformed);

  DiagonalMatrix eigenvalues;
  Matrix eigenvectors;
  diagonalise(s_transformed, eigenvalues, eigenvectors);

  // ignore small and/or negative eigenvalues
  // first half contains imaginary part
  int imax = ncols-1;
  for(; imax >= ncols/2; --imax)
    if(eigenvalues.element(imax, imax) < 1.0e+8) break; // checking e = 1/h (h > 1.0e-8)
  int imin = ncols/2;
  for(; imin <  ncols;   ++imin)
    if(eigenvalues.element(imin, imin) >= 1.0e-16) break;

  if((imax - imin + 1) < ncols/2)
    pout << "\t\t\t Warning: found and ignored zero frequency modes in RPA" << endl;

  ncols = imax - imin + 1;

  freq.ReSize(ncols);
  Matrix eigvectmp(eigenvectors.Nrows(), ncols);

  int icol = 0;
  for(int i = imax; i >= imin; --i, ++icol) {
    freq.element(icol, icol) = 1.0/eigenvalues.element(i, i);
    for(int j = 0; j < eigenvectors.Nrows(); ++j) {
      eigvectmp.element(j, icol) = eigenvectors.element(j, i)/sqrt(eigenvalues.element(i, i));
    }
  }

  Matrix alphaRe(2*nrows, ncols); alphaRe = 0.0;
  MatrixMultiply(u_matrix, 'n', eigvectmp, 'n', alphaRe, 1.0);

  alpha.ReSize(2*nrows, 2*ncols); alpha = 0.0;
  for(int i = 0; i < nrows; ++i) {
    int ix = 2*i;
    int iy = 2*i+1;
    for(int j = 0; j < ncols; ++j) {
      int jx = 2*j;
      int jy = 2*j+1;
      // Xnew(j) = sum{i} alphaRe(i,j) [ Xold(i)  Yold(i) ] ... (X index appears even #)
      alpha.element(ix, jx) = alphaRe.element(i,       j);
      alpha.element(iy, jx) = alphaRe.element(i+nrows, j);
      // Ynew(j) = sum{i} alphaRe(i,j) [ Yold(i)  Xold(i) ] ... (Y index appears odd  #)
      alpha.element(ix, jy) = alphaRe.element(i+nrows, j);
      alpha.element(iy, jy) = alphaRe.element(i,       j);
    }
  }

// CHECK DIAGONALIZE
//{
//  pout << "State Rotation Matrix (real): " << endl << alphaRe << endl;
//  pout << "State Rotation Matrix (real/imag): " << endl << alpha << endl;
//
//  Matrix check_h_mat(2*ncols, 2*ncols); check_h_mat = 0.0;
//  Matrix check_h_tmp(2*nrows, 2*ncols); check_h_tmp = 0.0;
//  MatrixMultiply(h_matrix, 'n', alpha, 'n', check_h_tmp, 1.0);
//  MatrixMultiply(alpha, 't', check_h_tmp, 'n', check_h_mat, 1.0);
//  pout << "Diagonalized H-matrix: " << endl << check_h_mat << endl;
//
//  Matrix check_s_mat(2*ncols, 2*ncols); check_s_mat = 0.0;
//  Matrix check_s_tmp(2*nrows, 2*ncols); check_s_tmp = 0.0;
//  MatrixMultiply(s_matrix, 'n', alpha, 'n', check_s_tmp, 1.0);
//  MatrixMultiply(alpha, 't', check_s_tmp, 'n', check_s_mat, 1.0);
//  pout << "Diagonalized S-matrix: " << endl << check_s_mat << endl;
//}
// DEBUG END

  return ncols;
}

void SpinAdapted::LRT::solve_wavefunction
(vector<Wavefunction>& psix1st, vector<Wavefunction>& psix2nd, const vector<double>& eigv, vector<double>& rnorm, SpinBlock& big,
 const guessWaveTypes& guesswavetype, const bool& onedot, const bool& dot_with_sys,
 const bool& rpa_sweep, const bool& rpa_sweep_2nd, int nroots, int mroots, int kroots)
{
//pout << "DEBUG @ LRT::solve_wavefunction: nroots = " << nroots << ", mroots = " << mroots << ", kroots = " << kroots << endl;

  int Nroots = (rpa_sweep ? 2 * nroots - 1 : nroots);
  int Mroots = (rpa_sweep ? 2 * mroots - 1 : mroots);
  int Kroots = (rpa_sweep ? 2 * kroots - 1 : kroots);

  DiagonalMatrix h_diag;
  bool useprecond = true;

  h_diag.ReSize(big.get_stateInfo().totalStates); h_diag = 0;
  if (dmrginp.outputlevel() > 0)
    pout << "\t\t\t Building Diagonal Hamiltonian " << endl;
  big.diagonalH (h_diag);
  if (dmrginp.outputlevel() > 0)
    pout << "\t\t\t Done building diagonal hamiltonian "<<endl;
  FORTINT m, n = 1, nsize = h_diag.Storage();
  pout << "\t\t\t Number of elements in wavefunction :: " << h_diag.Ncols() << endl;
  if (mpigetrank()==0) {
    m = idamax_(nsize,h_diag.Store(), n); 
    if (dmrginp.outputlevel() > 0)
      pout << "highest diagonal value " << m << " " << h_diag(m) << endl;
  }
  else {
    h_diag.ReSize(0);
  }
  
  psix1st.resize(Mroots+Nroots-Kroots);
  psix2nd.resize(Mroots+Nroots-Kroots);

  if (!rpa_sweep_2nd)
    GuessWave::LRT::transform_gauge(psix1st, psix2nd, Mroots,               big, guesswavetype, onedot, dot_with_sys);
  else
    GuessWave::LRT::transform_gauge(psix1st, psix2nd, Mroots+Nroots-Kroots, big, guesswavetype, onedot, dot_with_sys);

  multiply_h_total davidson_f(big, onedot);

  if (!rpa_sweep_2nd && mroots == 1)
    compute_guess_solutions(psix1st, h_diag, davidson_f, rpa_sweep, useprecond, nroots);

  double ezero = eigv[0];
  for(int i = 1; i <= h_diag.Ncols(); ++i)
    h_diag(i) -= ezero;

  if (!rpa_sweep_2nd) {
    if (rpa_sweep)
      LRT::RPA::solve_correction_equation(psix1st, eigv, rnorm, h_diag, davidson_f, useprecond, nroots, mroots, kroots);
    else
      LRT::TDA::solve_correction_equation(psix1st, eigv, rnorm, h_diag, davidson_f, useprecond, nroots, mroots, kroots);
  }
}

//
// compute initial guesses for excited states
//
void SpinAdapted::LRT::compute_guess_solutions
(vector<Wavefunction>& psix, DiagonalMatrix& h_diag, Davidson_functor& h_mult, bool rpa_sweep, bool useprecond, int nroots)
{
  pout.precision(12);
#ifndef SERIAL
  mpi::communicator world;
#endif

  dmrginp.hmultiply -> start();

  vector<Wavefunction> psix_tmp(nroots);
  Wavefunction  psi0;
  Wavefunction *psi0_ptr = &psi0;

  if(mpigetrank() == 0) {
    psix_tmp[0] = psix[0];
    psi0_ptr = &psix_tmp[0];
  }

#ifndef SERIAL
  mpi::broadcast(world, *psi0_ptr, 0);
#endif

  // guess wavefunctions from Krylov subspace

  for(int i = 1; i < nroots; ++i) {
    Wavefunction  psix_i;
    Wavefunction *psix_ptr = &psix_i;

    if(mpigetrank() == 0) {
      psix_tmp[i] = psix_tmp[0];
      psix_tmp[i].Clear();
      psix_ptr = &psix_tmp[i];
    }
    else {
      psix_i = psi0;
      psix_i.Clear();
    }

    h_mult(*psi0_ptr, *psix_ptr, 0);

    if(mpigetrank() == 0) {
      int success = 0;
      for(int j = 0; j < i; ++j) {
        Normalise(psix_tmp[i], &success);
        double overlap = DotProduct(psix_tmp[j], psix_tmp[i]);
        ScaleAdd(-overlap, psix_tmp[j], psix_tmp[i]);
      }
      Normalise(psix_tmp[i], &success);

      psi0_ptr = &psix_tmp[i];
    }
#ifndef SERIAL
    mpi::broadcast(world, *psi0_ptr, 0);
#endif
  }

  int nvals;
  if(mpigetrank() == 0) nvals = h_diag.Ncols();
#ifndef SERIAL
  mpi::broadcast(world, nvals, 0);
#endif

  if(nvals > 2*nroots)
    Linear::block_davidson(psix_tmp, h_diag, 1.0e-6, false, h_mult, useprecond);

  if(mpigetrank() == 0) {
    for(int i = 1; i < nroots; ++i) {
      int ix = (rpa_sweep ? 2*i-1 : i);
      psix[ix] = psix_tmp[i];
      double overlap = DotProduct(psix[0], psix[ix]);
      ScaleAdd(-overlap, psix[0], psix[ix]);

      if(rpa_sweep) {
        psix[ix+1] = psix[0];
        psix[ix+1].Clear();
      }
    }
  }

}

void SpinAdapted::LRT::TDA::solve_correction_equation
(vector<Wavefunction>& psix, const vector<double>& eigv, vector<double>& rnorm, DiagonalMatrix& h_diag, Davidson_functor& h_mult,
 bool useprecond, int nroots, int mroots, int kroots)
{
//pout << "DEBUG @ LRT::TDA::solve_correction_equation: nroots = " << nroots << ", mroots = " << mroots << ", kroots = " << kroots << endl;
  pout.precision(12);
#ifndef SERIAL
  mpi::communicator world;
#endif
//const int lroots = mroots + nroots - kroots;
  double levelshift = 0.0;

  dmrginp.hmultiply -> start();

  if(mpigetrank() == 0) {
    assert(psix.size() == mroots+nroots-kroots);
  }

  // Multiply all trial vectors with Hamiltonian: V(J) = H(00) C(J) + H(0J) C(0)
  Wavefunction  psi0;
  Wavefunction *psi0_ptr = &psi0;

  if(mpigetrank() == 0) {
    psi0_ptr = &psix[0];
  }

#ifndef SERIAL
  mpi::broadcast(world, *psi0_ptr, 0);
#endif

  // guess wavefunctions from Krylov subspace

  if(mroots == 1) {
//  for(int i = 1; i < nroots; ++i) {
//    Wavefunction  psix_i;
//    Wavefunction *psix_ptr = &psix_i;

//    if(mpigetrank() == 0) {
//      psix[i] = psix[0];
//      psix[i].Clear();
//      psix_ptr = &psix[i];
//    }
//    else {
//      psix_i = psi0;
//      psix_i.Clear();
//    }

//    h_mult(*psi0_ptr, *psix_ptr, 0);

//    if(mpigetrank() == 0) {
//      ScaleAdd(-eigv[0], psix[0], psix[i]); // FIXME: (H - E0) * X
//      int success = 0;
//      for(int j = 0; j < i; ++j) {
//        Normalise(psix[i], &success);
//        double overlap = DotProduct(psix[j], psix[i]);
//        ScaleAdd(-overlap, psix[j], psix[i]);
//      }
//      Normalise(psix[i], &success);

//      ScaleAdd(-eigv[0], psix[0], psix[i]); // FIXME: (H - E0) * X
//      double overlap = DotProduct(psix[0], psix[i]);
//      ScaleAdd(-overlap, psix[0], psix[i]);
//      int success = 0;
//      Normalise(psix[i], &success);


//      psi0_ptr = &psix[i];
//    }

#ifndef SERIAL
//    mpi::broadcast(world, *psi0_ptr, 0);
#endif
//  }
    return;
  }

  // compute sigma vectors

  vector<Wavefunction> sgvx(nroots);
  vector<Wavefunction> sgv0(nroots);

  for(int i = 1; i < nroots; ++i) {
    Wavefunction  psix_i;
    Wavefunction *psix_ptr = &psix_i;

    Wavefunction  sgvx_i;
    Wavefunction *sgvx_ptr = &sgvx_i;

    Wavefunction  sgv0_i;
    Wavefunction *sgv0_ptr = &sgv0_i;

    if(mpigetrank() == 0) {
      sgvx[i] = psix[i];
      sgvx[i].Clear();
      sgvx_ptr = &sgvx[i];

      sgv0[i] = psix[i];
      sgv0[i].Clear();
      sgv0_ptr = &sgv0[i];

      psix_ptr = &psix[i];
    }

#ifndef SERIAL
    mpi::broadcast(world, *psix_ptr, 0);
#endif
    if(mpigetrank() != 0) {
      sgvx_i = psix_i;
      sgvx_i.Clear();

      sgv0_i = psix_i;
      sgv0_i.Clear();
    }

    h_mult(*psix_ptr, *sgvx_ptr, 0);
    h_mult(*psi0_ptr, *sgv0_ptr, i);

    if(mpigetrank() == 0) {
      ScaleAdd(-eigv[0], psix[i], sgvx[i]); // FIXME: (H - E0) * X
      sgvx[i] += sgv0[i];
    }
  }

  dmrginp.hmultiply -> stop();

  // solve correction
  if(mpigetrank() == 0) {
    Wavefunction r;
    for(int i = 1; i < nroots; ++i) {
      r = sgvx[i];
      ScaleAdd(-eigv[i], psix[i], r);
      rnorm[i] += DotProduct(r, r);

      if(i < kroots) continue;

      if(useprecond)
        SpinAdapted::Linear::olsenPrecondition(r, psix[i], eigv[i], h_diag, levelshift);

//    r += psix[mroots+i-kroots]; // this comes from previous block (FIXME: not exact)
      double overlap = DotProduct(r, psix[0]);
      ScaleAdd(-overlap, psix[0], r);

      int success = 0;
      Normalise(r, &success);
      Scale(1.0/dmrginp.last_site(), r); // scaled by 1/k (not necessary)
      psix[mroots+i-kroots] = r;
    }
  }

}

void SpinAdapted::LRT::TDA::compute_matrix_elements
(vector<Wavefunction>& psix, const vector<double>& eigv, Davidson_functor& h_mult, Matrix& a_subspace, Matrix& s_subspace, int mroots)
{
//pout << "DEBUG @ LRT::solve_correction_equation: mroots = " << mroots << endl;
  pout.precision(12);
#ifndef SERIAL
  mpi::communicator world;
#endif
  double levelshift = 0.0;

  vector<Wavefunction> sgvx(mroots);
  vector<Wavefunction> sgv0(mroots);

  dmrginp.hmultiply -> start();

  if(mpigetrank() == 0) {
    assert(psix.size() == mroots);
  }

  // Multiply all trial vectors with Hamiltonian: V(J) = H(00) C(J) + H(0J) C(0)
  Wavefunction  psi0;
  Wavefunction *psi0_ptr = &psi0;

  if(mpigetrank() == 0) {
    psi0_ptr = &psix[0];
  }

#ifndef SERIAL
  mpi::broadcast(world, *psi0_ptr, 0);
#endif

//pout << "DEBUG @ LRT::solve_correction_equation: check point 1" << endl;
  for(int i = 1; i < mroots; ++i) {
    Wavefunction  psix_i;
    Wavefunction *psix_ptr = &psix_i;

    Wavefunction  sgvx_i;
    Wavefunction *sgvx_ptr = &sgvx_i;

    Wavefunction  sgv0_i;
    Wavefunction *sgv0_ptr = &sgv0_i;

    if(mpigetrank() == 0) {
      sgvx[i] = psix[i];
      sgvx[i].Clear();
      sgvx_ptr = &sgvx[i];

      sgv0[i] = psix[i];
      sgv0[i].Clear();
      sgv0_ptr = &sgv0[i];

      psix_ptr = &psix[i];
    }

#ifndef SERIAL
    mpi::broadcast(world, *psix_ptr, 0);
#endif
    if(mpigetrank() != 0) {
      sgvx_i = psix_i;
      sgvx_i.Clear();

      sgv0_i = psix_i;
      sgv0_i.Clear();
    }

    h_mult(*psix_ptr, *sgvx_ptr, 0);
    h_mult(*psi0_ptr, *sgv0_ptr, i);

    if(mpigetrank() == 0)
      ScaleAdd(-eigv[0], psix[i], sgvx[i]); // FIXME: (H - E0) * X
  }

  dmrginp.hmultiply -> stop();

  if(mpigetrank() == 0) {
    // compute Xi S Xj
    for(int i = 1; i < mroots; ++i) {
      double sii = DotProduct(psix[i], psix[i]);
      s_subspace(i, i) += sii;
      for(int j = 1; j < i; ++j) {
        double sij = DotProduct(psix[i], psix[j]);
        s_subspace(i, j) += sij;
        s_subspace(j, i) += sij;
      }
    }

    // compute Xi A Xj
    for(int i = 1; i < mroots; ++i) {
      double aii = DotProduct(psix[i], sgvx[i])
                 + DotProduct(psix[i], sgv0[i]) * 2.0;
      a_subspace(i, i) += aii;
      for(int j = 1; j < i; ++j) {
        double aij = DotProduct(psix[i], sgvx[j])
                   + DotProduct(psix[i], sgv0[j])
                   + DotProduct(psix[j], sgv0[i]);
        a_subspace(i, j) += aij;
        a_subspace(j, i) += aij;
      }
    }
  }
}

//
// for RPA
//

//
// NOTE:
//   psix[ 0] contains 0-th order wavefunction
//   psix[ix] contains 1-st order wavefunction X (real) { ix = 1,3,5,7,... }
//   psix[iy] contains 1-st order wavefunction Y (imag) { iy = 2,4,6,8,... }
//

void SpinAdapted::LRT::RPA::solve_correction_equation
(vector<Wavefunction>& psix, const vector<double>& eigv, vector<double>& rnorm,
 DiagonalMatrix& h_diag, Davidson_functor& h_mult, bool useprecond, int nroots, int mroots, int kroots)
{
//pout << "DEBUG @ LRT::solve_correction_equation: nroots = " << nroots << ", mroots = " << mroots << ", kroots = " << kroots << endl;
  pout.precision(12);
#ifndef SERIAL
  mpi::communicator world;
#endif
  double levelshift = 0.0;

  int lroots = mroots + nroots - kroots;

  int Nroots = 2 * nroots - 1;
  int Mroots = 2 * mroots - 1;
  int Kroots = 2 * kroots - 1;
  int Lroots = 2 * lroots - 1;

  dmrginp.hmultiply -> start();

  if(mpigetrank() == 0) {
    assert(psix.size() == Lroots);
  }

  // Multiply all trial vectors with Hamiltonian: V(J) = H(00) C(J) + H(0J) C(0)
  Wavefunction  psi0;
  Wavefunction *psi0_ptr = &psi0;

  if(mpigetrank() == 0) {
    psi0_ptr = &psix[0];
  }

#ifndef SERIAL
  mpi::broadcast(world, *psi0_ptr, 0);
#endif

  // guess wavefunctions from Krylov subspace

  if(mroots == 1) {

//  for(int i = 1; i < nroots; ++i) {

//    int ix = 2*i-1;
//    int iy = 2*i;

//    Wavefunction  psix_ix;
//    Wavefunction *psix_ptr = &psix_ix;

//    if(mpigetrank() == 0) {
//      psix[ix] = psix[0];
//      psix[ix].Clear();
//      psix_ptr = &psix[ix];

//      // set Y[i] to zero
//      psix[iy] = psix[0];
//      psix[iy].Clear();
//    }
//    else {
//      psix_ix = psi0;
//      psix_ix.Clear();
//    }

//    h_mult(*psi0_ptr, *psix_ptr, 0);

//    if(mpigetrank() == 0) {
//      ScaleAdd(-eigv[0], psix[0], psix[ix]); // FIXME: (H - E0) * X
//      double overlap = DotProduct(psix[0], psix[ix]);
//      ScaleAdd(-overlap, psix[0], psix[ix]);
//      int success = 0;
//      Normalise(psix[ix], &success);

//      psix[iy] = psix[ix];
//      Scale(0.0, psix[iy]);

//      psi0_ptr = &psix[ix];
//    }

#ifndef SERIAL
//    mpi::broadcast(world, *psi0_ptr, 0);
#endif
//  }
    return;
  }

  // compute sigma vectors

  vector<Wavefunction> sgvx(Nroots);
  vector<Wavefunction> sgvm(Nroots);

  for(int i = 1; i < Nroots; ++i) {

    Wavefunction  psix_i;
    Wavefunction *psix_ptr = &psix_i;

    // A0 * X1
    Wavefunction  sgvx_i;
    Wavefunction *sgvx_ptr = &sgvx_i;

    // AX * C0
    Wavefunction  sgv0_i;
    Wavefunction *sgv0_ptr = &sgv0_i;

    // BX * C0
    Wavefunction  sgvm_i;
    Wavefunction *sgvm_ptr = &sgvm_i;

    if(mpigetrank() == 0) {
      sgvx[i] = psix[i];
      sgvx[i].Clear();
      sgvx_ptr = &sgvx[i];

      sgv0_i = psix[i];
      sgv0_i.Clear();

      sgvm[i] = psix[i];
      sgvm[i].Clear();
      sgvm_ptr = &sgvm[i];

      psix_ptr = &psix[i];
    }

#ifndef SERIAL
    mpi::broadcast(world, *psix_ptr, 0);
#endif

    if(mpigetrank() != 0) {
      sgvx_i = psix_i;
      sgvx_i.Clear();

      sgv0_i = psix_i;
      sgv0_i.Clear();

      sgvm_i = psix_i;
      sgvm_i.Clear();
    }

    h_mult(*psix_ptr, *sgvx_ptr, 0);
    h_mult(*psi0_ptr, *sgv0_ptr, i);
    h_mult(*psi0_ptr, *sgvm_ptr, i, true);

    if(mpigetrank() == 0) {
      // diag.(A0 - E0) * X1
      ScaleAdd(-eigv[0], psix[i], sgvx[i]);
      sgvx[i] += sgv0_i;
    }
  }

  dmrginp.hmultiply -> stop();

  // build residual
  if(mpigetrank() == 0) {

    Wavefunction rx;
    Wavefunction ry;

    for(int i = 1; i < nroots; ++i) {

      int ix = 2*i-1;
      int iy = 2*i;

      rx  = sgvx[ix];
      rx += sgvm[iy]; // AX + BY
      ScaleAdd(-eigv[i], psix[ix], rx); // (H - E0) Xi - Ei S Xi

      ry  = sgvx[iy];
      ry += sgvm[ix]; // AY + BX
      ScaleAdd(+eigv[i], psix[iy], ry); // (H - E0) Yi + Ei S Yi

      rnorm[i] += DotProduct(rx, rx);
      rnorm[i] += DotProduct(ry, ry);

      if(i < kroots) continue;

      if(useprecond) {
        SpinAdapted::Linear::olsenPrecondition(rx, psix[ix], eigv[i], h_diag, levelshift);
        SpinAdapted::Linear::olsenPrecondition(ry, psix[iy],-eigv[i], h_diag, levelshift);
      }

      int success;
      double overlap;

//    rx += psix[Mroots+ix-Kroots]; // this comes from previous block (FIXME: not exact)
      overlap = DotProduct(rx, psix[0]);
      ScaleAdd(-overlap, psix[0], rx);
//    success = 0;
//    Normalise(rx, &success);

//    ry += psix[Mroots+iy-Kroots]; // this comes from previous block (FIXME: not exact)
      overlap = DotProduct(ry, psix[0]);
      ScaleAdd(-overlap, psix[0], ry);
//    success = 0;
//    Normalise(ry, &success);

      double xnorm = DotProduct(rx, rx);
//    double ynorm = DotProduct(ry, ry);

//    double scale = 1.0/dmrginp.last_site();
      if(xnorm > 1.0) {
        double scale = 1.0/sqrt(xnorm);
        Scale(scale, rx);
        Scale(scale, ry);
      }

      psix[Mroots+ix-Kroots] = rx;
      psix[Mroots+iy-Kroots] = ry;
    }
  }

}

void SpinAdapted::LRT::RPA::compute_matrix_elements
(vector<Wavefunction>& psix, const vector<double>& eigv, vector<double>& ynorm, Davidson_functor& h_mult,
 Matrix& a_subspace, Matrix& b_subspace, Matrix& s_subspace, Matrix& d_subspace, const bool& rpa_sweep_2nd, int mroots)
{
//pout << "DEBUG @ LRT::RPA::compute_matrix_elements: mroots = " << mroots << (rpa_sweep_2nd ? "  2-nd components" : "  1-st components") << endl;
  pout.precision(12);
#ifndef SERIAL
  mpi::communicator world;
#endif
  double levelshift = 0.0;

  int Mroots = 2*mroots-1;

  vector<Wavefunction> sgvx(Mroots);
  vector<Wavefunction> sgv0(Mroots);

  dmrginp.hmultiply -> start();

  if(mpigetrank() == 0) {
    assert(psix.size()   == Mroots);
    assert(ynorm.size()  == mroots);
  }

  // Multiply all trial vectors with Hamiltonian: V(J) = H(00) C(J) + H(0J) C(0)
  Wavefunction  psi0;
  Wavefunction *psi0_ptr = &psi0;

  if(mpigetrank() == 0) {
    psi0_ptr = &psix[0];
  }

#ifndef SERIAL
  mpi::broadcast(world, *psi0_ptr, 0);
#endif

  for(int i = 1; i < Mroots; ++i) {
    Wavefunction  psix_i;
    Wavefunction *psix_ptr = &psix_i;

    Wavefunction  sgvx_i;
    Wavefunction *sgvx_ptr = &sgvx_i;

    Wavefunction  sgv0_i;
    Wavefunction *sgv0_ptr = &sgv0_i;

    if(mpigetrank() == 0) {
      sgvx[i] = psix[i];
      sgvx[i].Clear();
      sgvx_ptr = &sgvx[i];

      sgv0[i] = psix[i];
      sgv0[i].Clear();
      sgv0_ptr = &sgv0[i];

      psix_ptr = &psix[i];
    }

#ifndef SERIAL
    mpi::broadcast(world, *psix_ptr, 0);
#endif
    if(mpigetrank() != 0) {
      sgvx_i = psix_i;
      sgvx_i.Clear();

      sgv0_i = psix_i;
      sgv0_i.Clear();
    }

    if(!rpa_sweep_2nd) {
      h_mult(*psix_ptr, *sgvx_ptr, 0); // A0 * X1
      h_mult(*psi0_ptr, *sgv0_ptr, i); // AX * C0
      if(mpigetrank() == 0)
        ScaleAdd(-eigv[0], psix[i], sgvx[i]); // diag.(A - E0) * X1
    }
    else {
      h_mult(*psi0_ptr, *sgv0_ptr, i, true); // BX * C0
    }
  }

  dmrginp.hmultiply -> stop();

  // compute matrix elements
  if(mpigetrank() == 0) {

  if(!rpa_sweep_2nd) {

    // [X Y][ S  D ][X]
    //      [-D -S ][Y]
    for(int i = 1; i < mroots; ++i) {

      int ix = 2*i-1;
      int iy = 2*i;

      double sxixi = DotProduct(psix[ix], psix[ix]);
      double syiyi = DotProduct(psix[iy], psix[iy]);

      s_subspace(i, i) += sxixi;
      s_subspace(i, i) -= syiyi;

      ynorm[i] += syiyi;

      for(int j = 1; j < i; ++j) {

        int jx = 2*j-1;
        int jy = 2*j;

        double sxixj = DotProduct(psix[ix], psix[jx]);
        double syiyj = DotProduct(psix[iy], psix[jy]);

        double sxiyj = DotProduct(psix[ix], psix[jy]);
        double syixj = DotProduct(psix[iy], psix[jx]);

        s_subspace(i, j) += sxixj;
        s_subspace(j, i) += sxixj;

        s_subspace(i, j) -= syiyj;
        s_subspace(j, i) -= syiyj;

        d_subspace(i, j) += (sxiyj - syixj);
        d_subspace(j, i) += (syixj - sxiyj);
      }
    }

    // [X Y][ A  0 ][X]
    //      [ 0  A ][Y]
    for(int i = 1; i < mroots; ++i) {

      int ix = 2*i-1;
      int iy = 2*i;

      double axixi = DotProduct(psix[ix], sgvx[ix])
                   + DotProduct(psix[ix], sgv0[ix]) * 2.0;
      double ayiyi = DotProduct(psix[iy], sgvx[iy])
                   + DotProduct(psix[iy], sgv0[iy]) * 2.0;

      double axiyi = DotProduct(psix[ix], sgvx[iy])
                   + DotProduct(psix[ix], sgv0[iy])
                   + DotProduct(psix[iy], sgv0[ix]);

      a_subspace(i, i) += (axixi + ayiyi);

      b_subspace(i, i) += (2.0   * axiyi);

      for(int j = 1; j < i; ++j) {

        int jx = 2*j-1;
        int jy = 2*j;

        double axixj = DotProduct(psix[ix], sgvx[jx])
                     + DotProduct(psix[ix], sgv0[jx])
                     + DotProduct(psix[jx], sgv0[ix]);
        double ayiyj = DotProduct(psix[iy], sgvx[jy])
                     + DotProduct(psix[iy], sgv0[jy])
                     + DotProduct(psix[jy], sgv0[iy]);

        double axiyj = DotProduct(psix[ix], sgvx[jy])
                     + DotProduct(psix[ix], sgv0[jy])
                     + DotProduct(psix[jx], sgv0[iy]);
        double ayixj = DotProduct(psix[iy], sgvx[jx])
                     + DotProduct(psix[iy], sgv0[jx])
                     + DotProduct(psix[jy], sgv0[ix]);

        a_subspace(i, j) += (axixj + ayiyj);
        a_subspace(j, i) += (axixj + ayiyj);

        b_subspace(i, j) += (axiyj + ayixj);
        b_subspace(j, i) += (axiyj + ayixj);
      }
    }
  } // end if (rpa_sweep_2nd)
  else 
  {
    // [X Y][ 0  B ][X]
    //      [ B  0 ][Y]
    for(int i = 1; i < mroots; ++i) {

      int ix = 2*i-1;
      int iy = 2*i;

      double bxiyi = DotProduct(psix[ix], sgv0[iy]) * 2.0;
      double byixi = DotProduct(psix[iy], sgv0[ix]) * 2.0;

      double bxixi = DotProduct(psix[ix], sgv0[ix]) * 2.0;
      double byiyi = DotProduct(psix[iy], sgv0[iy]) * 2.0;

      a_subspace(i, i) += (bxiyi + byixi);

      b_subspace(i, i) += (bxixi + byiyi);

      for(int j = 1; j < i; ++j) {

        int jx = 2*j-1;
        int jy = 2*j;

        double bxiyj = DotProduct(psix[ix], sgv0[jy])
                     + DotProduct(psix[jx], sgv0[iy]);
        double byixj = DotProduct(psix[iy], sgv0[jx])
                     + DotProduct(psix[jy], sgv0[ix]);

        double bxixj = DotProduct(psix[ix], sgv0[jx])
                     + DotProduct(psix[jx], sgv0[ix]);
        double byiyj = DotProduct(psix[iy], sgv0[jy])
                     + DotProduct(psix[jy], sgv0[iy]);

        a_subspace(i, j) += (bxiyj + byixj);
        a_subspace(j, i) += (bxiyj + byixj);

        b_subspace(i, j) += (bxixj + byiyj);
        b_subspace(j, i) += (bxixj + byiyj);
      }
    }
  } // end if (rpa_sweep_2nd) else

  } // end if (mpigetrank() == 0)
}

