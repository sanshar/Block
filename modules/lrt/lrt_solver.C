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

void SpinAdapted::LRT::solve_wavefunction
(vector<Wavefunction>& psix, const vector<double>& eigv, vector<double>& rnorm, SpinBlock& big,
 const guessWaveTypes& guesswavetype, const bool& onedot, const bool& dot_with_sys, int nroots, int mroots, int kroots)
{
//pout << "DEBUG @ LRT::solve_wavefunction: nroots = " << nroots << ", mroots = " << mroots << ", kroots = " << kroots << endl;
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
  else 
    h_diag.ReSize(0);
  
  // mroots: # of trial vecs, kroots: # of roots to be updated
  psix.resize(mroots+nroots-kroots);
  multiply_h_total davidson_f(big, onedot);
  GuessWave::LRT::transform_gauge(psix, mroots, big, guesswavetype, onedot, dot_with_sys); // FIXME
  TDA::solve_correction_equation(psix, eigv, rnorm, h_diag, davidson_f, useprecond, nroots, mroots, kroots); // FIXME
}

void SpinAdapted::LRT::TDA::solve_correction_equation
(vector<Wavefunction>& psix, const vector<double>& eigv, vector<double>& rnorm, DiagonalMatrix& h_diag, Davidson_functor& h_mult,
 bool useprecond, int nroots, int mroots, int kroots)
{
//pout << "DEBUG @ LRT::solve_correction_equation: nroots = " << nroots << ", mroots = " << mroots << ", kroots = " << kroots << endl;
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
    for(int i = 1; i < nroots; ++i) {
      Wavefunction  psix_i;
      Wavefunction *psix_ptr = &psix_i;

      if(mpigetrank() == 0) {
        psix[i] = psix[0];
        psix[i].Clear();
        psix_ptr = &psix[i];
      }
      else {
        psix_i = psi0;
        psix_i.Clear();
      }

      h_mult(*psi0_ptr, *psix_ptr, 0);

      if(mpigetrank() == 0) {
        int success = 0;
        for(int j = 0; j < i; ++j) {
          Normalise(psix[i], &success);
          double overlap = DotProduct(psix[j], psix[i]);
          ScaleAdd(-overlap, psix[j], psix[i]);
        }
        Normalise(psix[i], &success);

        psi0_ptr = &psix[i];
      }

#ifndef SERIAL
      mpi::broadcast(world, *psi0_ptr, 0);
#endif
    }
    return;
  }

  // compute sigma vectors

  vector<Wavefunction> sgvx(mroots);

  for(int i = 1; i < mroots; ++i) {
    Wavefunction  psix_i;
    Wavefunction *psix_ptr = &psix_i;

    Wavefunction  sgvx_i;
    Wavefunction *sgvx_ptr = &sgvx_i;

    if(mpigetrank() == 0) {
      sgvx[i] = psix[i];
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
        SpinAdapted::Linear::olsenPrecondition(r, psix[i], eigv[i], h_diag, levelshift);

      int success = 0;
      Normalise(r, &success);
      double overlap = DotProduct(r, psix[0]);
      ScaleAdd(-overlap, psix[0], r);

      Normalise(r, &success);
      Scale(1.0/dmrginp.last_site(), r); // scaled by 1/k (not necessary)
      psix[mroots+i-kroots] = r;
    }
  }

//  // broadcast psix
//  for(int i = mroots; i < lroots; ++i) {
//    Wavefunction  psix_i;
//    Wavefunction *psix_ptr = &psix_i;
//    if(mpigetrank() == 0) {
//      psix_ptr = &psix[i];
//    }
//#ifndef SERIAL
//    mpi::broadcast(world, *psix_ptr, 0);
//#endif
//  }
}

void SpinAdapted::LRT::TDA::compute_matrix_elements
(vector<Wavefunction>& psix, Davidson_functor& h_mult, Matrix& h_subspace, Matrix& s_subspace, int mroots)
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
  }
//pout << "DEBUG @ LRT::solve_correction_equation: check point 1 - passed" << endl;

  dmrginp.hmultiply -> stop();

  // build residual
  if(mpigetrank() == 0) {
//pout << "DEBUG @ LRT::solve_correction_equation: check point 2" << endl;
    for(int i = 1; i < mroots; ++i) {
      double sii = DotProduct(psix[i], psix[i]);
      s_subspace(i, i) += sii;
      for(int j = 1; j < i; ++j) {
        double sij = DotProduct(psix[i], psix[j]);
        s_subspace(i, j) += sij;
        s_subspace(j, i) += sij;
      }
    }

//pout << "DEBUG @ LRT::solve_correction_equation: check point 3" << endl;
    for(int i = 1; i < mroots; ++i) {
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
//pout << "DEBUG @ LRT::solve_correction_equation: check point 4 ( done )" << endl;
}



