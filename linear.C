/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "linear.h"
#include "pario.h"
#include "global.h"
#include "MatrixBLAS.h"
#include "wavefunction.h"
#include "solver.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#ifndef SERIAL
#include <boost/mpi.hpp>
#endif
#include "pario.h"
using namespace boost;

const double EPS=1.e-20;


void SpinAdapted::Linear::precondition(Wavefunction& op, double e, DiagonalMatrix& diagonal, double levelshift)
{
  if (!mpigetrank())
  {
    int index = 1;
    for (int lQ = 0; lQ < op.nrows (); ++lQ)
      for (int rQ = 0; rQ < op.ncols (); ++rQ)
	if (op.allowed(lQ, rQ))
	  for (int lQState = 0; lQState < op.operator_element(lQ, rQ).Nrows (); ++lQState)
	    for (int rQState = 0; rQState < op.operator_element(lQ, rQ).Ncols (); ++rQState)
	      {
		if (fabs(e - diagonal(index)) > 1.e-12)		
		  op.operator_element(lQ, rQ).element(lQState, rQState) /= (e - diagonal(index)+levelshift);
		++index;
	      }
  }
#ifndef SERIAL
  //mpi::communicator world;
  //mpi::broadcast(world, op, 0);
#endif
}


void SpinAdapted::Linear::olsenPrecondition(Wavefunction& op, Wavefunction& C0, double e, DiagonalMatrix& diagonal, double levelshift)
{

  Wavefunction C0copy = C0;
  precondition(C0copy, e, diagonal, levelshift);
  double numerator = DotProduct(C0copy, op);
  double denominator = DotProduct(C0, C0copy);
  ScaleAdd(-numerator/denominator, C0, op);
  precondition(op, e, diagonal, levelshift);

}

void SpinAdapted::Linear::Lanczos(vector<Wavefunction>& b, DiagonalMatrix& e, double normtol, Davidson_functor& h_multiply, int nroots)
{
  int iter = 0;
  pout.precision(12);


  if(mpigetrank() == 0) {
    //b[0].Randomise();
    Normalise(b[0]);
    b.resize(dmrginp.max_lanczos_dimension());
  }


  double beta_minus = 1.0;
  std::vector<double> diagonal(dmrginp.max_lanczos_dimension(), 0);
  std::vector<double> offdiagonal(dmrginp.max_lanczos_dimension(), 0);

  Wavefunction b0;
  Wavefunction* bptr = &b0;
  if (mpigetrank() == 0)
    bptr = &b[0];
  else
    bptr = &b0;
  
#ifndef SERIAL
  mpi::communicator world;
  mpi::broadcast(world, *bptr, 0);
#endif
  Wavefunction r(*bptr);
  
  bool notconverged = true;
  while (notconverged && iter < dmrginp.max_lanczos_dimension()) {
    r.Clear();
    pout << "\t\t\t Lanczos Iteration :: "<<iter<<endl;
    iter ++;

#ifndef SERIAL
    mpi::broadcast(world, *bptr, 0);
#endif

    h_multiply(*bptr, r);


    if(mpigetrank() == 0) {
      if (iter > 1)
        ScaleAdd(-beta_minus, b[iter-2], r);


      diagonal[iter-1] = DotProduct(*bptr, r);


      double alphai = diagonal[iter-1];

      Matrix alpha;
      std::vector<double> diagonal_temporary(iter,0);
      std::vector<double> offdiagonal_temporary(iter,0);
      for (int i=0; i<iter; i++) {
        diagonal_temporary[i] = diagonal[i];
        offdiagonal_temporary[i] = offdiagonal[i];
      }

      diagonalise_tridiagonal(diagonal_temporary, offdiagonal_temporary, iter, alpha);

      ScaleAdd(-alphai, *bptr, r);
      for (int i=iter-1; i>=0; i--) {
	double overlap = DotProduct(b[i], r);
	ScaleAdd(-overlap, b[i], r);
      }
      offdiagonal[iter-1]  = sqrt(DotProduct(r, r));

      if (dmrginp.outputlevel() > 0) {
        if (iter > 1)
          notconverged = false;
        for (int i = 1; i <= min(iter-1, nroots); ++i) {
          pout << "\t\t\t " << i << " ::  " << diagonal_temporary[i-1]+dmrginp.get_coreenergy() <<"  "<<pow(offdiagonal[iter-1]*alpha(iter,i),2)<<endl;
          if (pow(offdiagonal[iter-1]*alpha(iter,i),2) >= normtol) {
            notconverged = true;
            break;
          }
        }

      }
      if (iter == dmrginp.max_lanczos_dimension()) {
        pout << "Reached the maximum number of allowed iterations!!"<<endl;
        exit(0);
      }


      //cout << offdiagonal[iter-1]<<"  "<<iter<<endl;                                                                                                                       
      if(!notconverged) {
	for (int i=0; i<nroots; i++)
	  e(i+1) = diagonal_temporary[i];

	vector<Wavefunction> btmp = b;
	for (int i=0; i<nroots; i++)
	  Scale(alpha.element(i,i), b[i]);
	for (int i=0; i<nroots; i++)
	  for (int j=0; j<iter; j++) 
	    if (i != j)
	      ScaleAdd(alpha.element(j,i), btmp[j], b[i]);

	break;
      }

      Normalise(r);
      b[iter] = r;
      if(mpigetrank() == 0)
	bptr = &b[iter];
    }
#ifndef SERIAL
    mpi::broadcast(world, notconverged, 0);
#endif
    r.Clear();
  }
}


void SpinAdapted::Linear::block_davidson(vector<Wavefunction>& b, DiagonalMatrix& h_diag, double normtol, const bool &warmUp, Davidson_functor& h_multiply, bool& useprecond)
{

  pout.precision (12);
#ifndef SERIAL
  mpi::communicator world;
#endif
  int iter = 0;
  double levelshift = 0.0;
  int nroots = b.size();
  //normalise all the guess roots
  if(mpigetrank() == 0) {
    for(int i=0;i<nroots;++i) {
      for(int j=0;j<i;++j) {
        double overlap = DotProduct(b[j], b[i]);
        ScaleAdd(-overlap, b[j], b[i]);
      }
      Normalise(b[i]);
    }
  }
  vector<Wavefunction> sigma;
  int converged_roots = 0;
  while(true) {
    if (dmrginp.outputlevel() > 0)
	  pout << "\t\t\t Davidson Iteration :: " << iter << endl;

    ++iter;
    dmrginp.hmultiply -> start();

    int sigmasize, bsize;

    if (mpigetrank() == 0) {
	  sigmasize = sigma.size();
	  bsize = b.size();
    }
#ifndef SERIAL
    mpi::broadcast(world, sigmasize, 0);
    mpi::broadcast(world, bsize, 0);
#endif
    //multiply all guess vectors with hamiltonian c = Hv

    for(int i=sigmasize;i<bsize;++i) {
	  Wavefunction sigmai, bi;
	  Wavefunction* sigmaptr=&sigmai, *bptr = &bi;

	  if (mpigetrank() == 0) {
	    sigma.push_back(b[i]);
	    sigma[i].Clear();
	    sigmaptr = &sigma[i];
	    bptr = &b[i];
	  }

#ifndef SERIAL
	  mpi::broadcast(world, *bptr, 0);
#endif
	  if (mpigetrank() != 0) {
	    sigmai = bi;
	    sigmai.Clear();
	  }

	  h_multiply(*bptr, *sigmaptr);
	}
    dmrginp.hmultiply -> stop();

    Wavefunction r;
    DiagonalMatrix subspace_eigenvalues;

    if (mpigetrank() == 0) {
	  Matrix subspace_h(b.size(), b.size());
	  for (int i = 0; i < b.size(); ++i)
	    for (int j = 0; j <= i; ++j) {
	      subspace_h.element(i, j) = DotProduct(b[i], sigma[j]);
	      subspace_h.element(j, i) = subspace_h.element(i, j);
	    }

	Matrix alpha;
	diagonalise(subspace_h, subspace_eigenvalues, alpha);
	
	if (dmrginp.outputlevel() > 0) {
	  for (int i = 1; i <= subspace_eigenvalues.Ncols (); ++i)
	    pout << "\t\t\t " << i << " ::  " << subspace_eigenvalues(i,i)+dmrginp.get_coreenergy() << endl;
	}	
	
	//now calculate the ritz vectors which are approximate eigenvectors
	vector<Wavefunction> btmp = b;
	vector<Wavefunction> sigmatmp = sigma;
	for (int i = 0; i < b.size(); ++i)
	  {
	    Scale(alpha.element(i, i), b[i]);
	    Scale(alpha.element(i, i), sigma[i]);
	  }
	for (int i = 0; i < b.size(); ++i)
	  for (int j = 0; j < b.size(); ++j)
	    {
	      if (i != j)
		{
		  ScaleAdd(alpha.element(i, j), btmp[i], b[j]);
		  ScaleAdd(alpha.element(i, j), sigmatmp[i], sigma[j]);
		}
	    }
	
	// build residual                                                                                                              
        for (int i=0; i<converged_roots; i++) {
          r = sigma[i];
          ScaleAdd(-subspace_eigenvalues(i+1), b[i], r);
          double rnorm = DotProduct(r,r);
          if (rnorm > normtol) {
            converged_roots = i;
            cout << "\t\t\t going back to converged root "<<i<<"  "<<rnorm<<" > "<<normtol<<endl;
            continue;
          }
        }
        r = sigma[converged_roots];
        ScaleAdd(-subspace_eigenvalues(converged_roots+1), b[converged_roots], r);
      }


      double rnorm;
      if (mpigetrank() == 0)
	rnorm = DotProduct(r,r);  

#ifndef SERIAL
      mpi::broadcast(world, converged_roots, 0);
      mpi::broadcast(world, rnorm, 0);
#endif

      if (useprecond && mpigetrank() == 0)
	olsenPrecondition(r, b[converged_roots], subspace_eigenvalues(converged_roots+1), h_diag, levelshift);


      if (dmrginp.outputlevel() > 0)
	pout << "\t \t \t residual :: " << rnorm << endl;
      if (rnorm < normtol)
	{
	  if (dmrginp.outputlevel() > 0)
	    pout << "\t\t\t Converged root " << converged_roots << endl;

	  ++converged_roots;
	  if (converged_roots == nroots)
	    {
	      if (mpigetrank() == 0) {

		for (int i = 0; i < min((int)(b.size()), h_diag.Ncols()); ++i)
		  h_diag.element(i) = subspace_eigenvalues.element(i);
	      }
	      break;
	    }
	}
      else if (mpigetrank() == 0)
	{
	  if(b.size() >= dmrginp.deflation_max_size())
	    {
	      if (dmrginp.outputlevel() > 0)
		pout << "\t\t\t Deflating block Davidson...\n";
	      b.resize(dmrginp.deflation_min_size());
	      sigma.resize(dmrginp.deflation_min_size());
	    }
	  for (int j = 0; j < b.size(); ++j)
	    {
	      //Normalize
	      double normalization = DotProduct(r, r);
	      Scale(1./sqrt(normalization), r);

	      double overlap = DotProduct(r, b[j]);
	      ScaleAdd(-overlap, b[j], r);
	    }
	  //double tau2 = DotProduct(r,r);

	  Normalise(r);
	  b.push_back(r);

	}


    }
}

