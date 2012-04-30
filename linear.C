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


void SpinAdapted::Linear::block_davidson(vector<Wavefunction>& b, DiagonalMatrix& h_diag, double normtol, const bool &warmUp, Davidson_functor& h_multiply, bool& useprecond, bool& solved)
{
  //mcheck("************* Memory in davidson routine ");

  pout.precision (12);
#ifndef SERIAL
  mpi::communicator world;
#endif
  int iter = 0;
  double levelshift = 0.0;
  int nroots = b.size();
  //normalise all the guess roots
  for(int i=0;i<nroots;++i)
    {
      for(int j=0;j<i;++j)
	{
	  double overlap = DotProduct(b[j], b[i]);
	  ScaleAdd(-overlap, b[j], b[i]);
	}
      Normalise(b[i]);
    }
  
  if(warmUp)
    {
      for(int i=0;i<nroots;++i)
	{
	  b.push_back(b[i]);
	  b[i+nroots].Randomise();
	  for(int j=0;j<i+nroots;++j)
	    {
	      Normalise(b[i+nroots]);
	      double overlap = DotProduct(b[j], b[i+nroots]);
	      ScaleAdd(-overlap, b[j], b[i+nroots]);
	    }
	  int success = 0;
	  Normalise(b[i+nroots], &success);
	}
    }

  vector<Wavefunction> sigma;
  int converged_roots = 0;
  while(true)
    {
      if (dmrginp.outputlevel() != 0)
	pout << "\t\t\t Davidson Iteration :: " << iter << endl;
      //mcheck(" ******************* ");
      ++iter;
      dmrginp.hmultiply.start();

      int sigmasize, bsize;

      if (mpigetrank() == 0) {
	sigmasize = sigma.size();
	bsize = b.size();
      }

      mpi::broadcast(world, sigmasize, 0);
      mpi::broadcast(world, bsize, 0);
      //multiply all guess vectors with hamiltonian c = Hv

      for(int i=sigmasize;i<bsize;++i)
	{
	  Wavefunction* sigmaptr, *bptr;
	  Wavefunction sigmai;
	  if (mpigetrank() == 0) {
	    sigma.push_back(b[i]);
	    sigma[i].Clear();
	    sigmaptr = &sigma[i];
	    bptr = &b[i];
	  }
	  else {
	    sigmai = b[0];
	    sigmai.Clear();
	    sigmaptr = &sigmai;
	    bptr = &b[converged_roots]; //b vector is not doing anything in procs with non-zero rank
	  }
	  mpi::broadcast(world, *bptr, 0);
	   
	  h_multiply(*bptr, *sigmaptr);
	  
	}
      dmrginp.hmultiply.stop();

      Wavefunction r;
      DiagonalMatrix subspace_eigenvalues;
      if (mpigetrank() == 0) {
	//generate the hamiltonian in b space h(i,j) = v_i H v_j
	Matrix subspace_h(b.size(), b.size());
	for (int i = 0; i < b.size(); ++i)
	  for (int j = 0; j <= i; ++j)
	    {
	      subspace_h.element(i, j) = DotProduct(b[i], sigma[j]);
	      subspace_h.element(j, i) = subspace_h.element(i, j);
	    }
	//diagonalise the h
	//DiagonalMatrix subspace_eigenvalues;
	Matrix alpha;
	diagonalise(subspace_h, subspace_eigenvalues, alpha);
	
	if (dmrginp.outputlevel() != 0) {
	  for (int i = 1; i <= subspace_eigenvalues.Ncols (); ++i)
	    pout << "\t\t\t " << i << " ::  " << subspace_eigenvalues(i,i) << endl;
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
	r = sigma[converged_roots];
	ScaleAdd(-subspace_eigenvalues(converged_roots+1), b[converged_roots], r);
      }

      double rnorm;
      if (mpigetrank() == 0)
	rnorm = DotProduct(r,r);
  
      mpi::broadcast(world, rnorm, 0);

      if (useprecond && mpigetrank() == 0)
	olsenPrecondition(r, b[converged_roots], subspace_eigenvalues(converged_roots+1), h_diag, levelshift);


      if (dmrginp.outputlevel() != 0)
	pout << "\t \t \t residual :: " << rnorm << endl;
      if (rnorm < normtol)
	{
	  if (dmrginp.outputlevel() != 0)
	    pout << "\t\t\t Converged root " << converged_roots << endl;
	  mpi::broadcast(world, b[converged_roots], 0);
	  ++converged_roots;	  
	  if (converged_roots == nroots)
	    {
	      if (mpigetrank() == 0) {
		for (int i = 0; i < min((int)(b.size()), h_diag.Ncols()); ++i)
		  h_diag.element(i) = subspace_eigenvalues.element(i);
	      }
	      mpi::broadcast(world, h_diag, 0);
	      solved = true;
	      break;
	    }
	}
      else if (mpigetrank() == 0)
	{
	  if(b.size() >= dmrginp.deflation_max_size())
	    {
	      if (dmrginp.outputlevel() != 0)
		pout << "\t\t\t Deflating block davidson...\n";
	      b.resize(dmrginp.deflation_min_size());
	      sigma.resize(dmrginp.deflation_min_size());
	    }
	  for (int j = 0; j < b.size(); ++j)
	    {
	      Normalise(r);
	      double overlap = DotProduct(r, b[j]);
	      ScaleAdd(-overlap, b[j], r);
	    }
	  //double tau2 = DotProduct(r,r);

	  Normalise(r);
	  b.push_back(r);

	}


    }
}

