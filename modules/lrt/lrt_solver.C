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

#include "lrt_solver.h"
#include "lrt_davidson.h"
#include "gauge_transform.h"

void SpinAdapted::LRT::solve_wavefunction
(vector<Wavefunction>& psix, const vector<double>& eigv, vector<double>& rnorm, SpinBlock& big, int nroots, int mroots, int kroots)
{
  DiagonalMatrix h_diag;
  bool useprecond = true;

  h_diag.ReSize(big.get_stateInfo().totalStates); h_diag = 0;
  if (dmrginp.outputlevel() > 0)
    pout << "\t\t\t Building Diagonal Hamiltonian " << endl;
  big.diagonalH (h_diag);
  if (dmrginp.outputlevel() > 0)
    pout << "\t\t\t Done building diagonal hamiltonian "<<endl;
  FORTINT m, n = 1, nsize = h_diag.Storage();
  pout << "\t\t\t Number of elements in wavefunction :: " << e.Ncols() << endl;
  if (mpigetrank()==0) {
    m = idamax_(nsize,h_diag.Store(), n); 
    if (dmrginp.outputlevel() > 0)
      pout << "highest diagonal value " << m << " " << h_diag(m) << endl;
  }
  else 
    h_diag.ReSize(0);
  
  // mroots: # of trial vecs, kroots: # of roots to be updated
  psix.resize(mroots+kroots);
  TDA::multiply_h_lr1 davidson_f(big, true);
  GaugeTransform::gauge_fixed_wave(psix, big, mroots); // FIXME
  TDA::solve_correction_equation(psix, eigv, h_diag, davidson_f, useprecond, nroots, mroots, kroots); // FIXME
}

void SpinAdapted::LRT::TDA::solve_correction_equation
(vector<Wavefunction>& psix, const vector<double>& eigv, vector<double>& rnorm, const DiagonalMatrix& h_diag, Davidson_functor& h_mult, bool useprecond, int nroots, int mroots, int kroots)
{
  pout.precision(12);
#ifndef SERIAL
  mpi::communicator world;
#endif
  const int lroots = mroots + nroots - kroots;
  double levelshift = 0.0;

  vector<Wavefunction> sgvx;

  dmrginp.hmultiply -> start();

  if(mpigetrank() == 0) {
    assert(psix.size() == mroots)
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

  for(int i = 1; i < mroots; ++i) {
    Wavefunction  psix_i;
    Wavefunction *psix_ptr = &psix_i;

    Wavefunction  sgvx_i;
    Wavefunction *sgvx_ptr = &sgvx_i;

    if(mpigetrank() == 0) {
      sgvx.push_back(psix[i]);
      sgvx[i].Clear();
      sgvx_ptr = &sgvx[i];

      psix_ptr = &psix[i];
    }

#ifndef SERIAL
    mpi::broadcast(world, *psix_ptr, 0);
#endif
    if(mpigetrank() != 0) {
      sgvx_i = psix_i;
      sgvx_i.Clear();
    }

    h_mult(*psix_ptr, *sgvx_ptr, 0);
    h_mult(*psi0_ptr, *sgvx_ptr, i);
  }

  dmrginp.hmultiply -> stop();

  // build residual
  if(mpigetrank() == 0) {
    Wavefunction r;
    for(int i = 1; i < nroots; ++i) {
      r = sgvx[i];
      ScaleAdd(-eigv[i], psix[i], r);
      rnorm[i] += DotProduct(r, r);

      if(i < kroots) continue;

      if(useprecond)
        olsenPrecondition(r, psix[i], eigv[i], h_diag, levelshift);

      Normalise(r);
      double overlap = DotProduct(r, psix[0]);
      ScaleAdd(-overlap, psix[0], r);

      Normalise(r);
      Scale(1.0/dmrginput.last_site(), r); // scaled by 1/k (not necessary)
      psix.push_back(r);
    }
  }
}

void SpinAdapted::LRT::TDA::compute_matrix_elements
(const vector<Wavefunction>& psix, Davidson_functor& h_mult, Matrix& h_subspace, Matrix& s_subspace, int mroots)
{
  pout.precision(12);
#ifndef SERIAL
  mpi::communicator world;
#endif
  double levelshift = 0.0;

  vector<Wavefunction> sgvx;
  vector<Wavefunction> sgv0;

  dmrginp.hmultiply -> start();

  if(mpigetrank() == 0) {
    assert(psix.size() >= mroots)
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

  for(int i = 1; i < mroots; ++i) {
    Wavefunction  psix_i;
    Wavefunction *psix_ptr = &psix_i;

    Wavefunction  sgvx_i;
    Wavefunction *sgvx_ptr = &sgvx_i;

    Wavefunction  sgv0_i;
    Wavefunction *sgv0_ptr = &sgv0_i;

    if(mpigetrank() == 0) {
      sgvx.push_back(psix[i]);
      sgvx[i].Clear();
      sgvx_ptr = &sgvx[i];

      sgv0.push_back(psix[i]);
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
  }

  dmrginp.hmultiply -> stop();

  // build residual
  if(mpigetrank() == 0) {
    for(int i = 1; i < nroots; ++i) {
      double sii = DotProduct(psix[i], psix[i]);
      s_subspace(i, i) += sii;
      for(int j = 1; j < i; ++j) {
        double sij = DotProduct(psix[i], psix[j]);
        s_subspace(i, j) += sij;
        s_subspace(j, i) += sij;
      }
    }

    for(int i = 1; i < nroots; ++i) {
      double hii = DotProduct(psix[i], sgvx[i])
                 + DotProduct(psix[i], sgv0[i]) * 2.0;
      h_subspace(i, i) += hii;
      for(int j = 1; j < i; ++j) {
        double hij = DotProduct(psix[i], sgvx[j])
                   + DotProduct(psix[i], sgv0[j])
                   + DotProduct(psix[j], sgv0[i]);
        h_subspace(i, j) += hij;
        h_subspace(j, i) += hij;
      }
    }
  }
}



