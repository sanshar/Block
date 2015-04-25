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
  p3out.precision(12);


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
    p3out << "\t\t\t Lanczos Iteration :: "<<iter<<endl;
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
          p3out << "\t\t\t " << i << " ::  " << diagonal_temporary[i-1] <<"  "<<pow(offdiagonal[iter-1]*alpha(iter,i),2)<<endl;
          if (pow(offdiagonal[iter-1]*alpha(iter,i),2) >= normtol) {
            notconverged = true;
            break;
          }
        }

      }
      if (iter == dmrginp.max_lanczos_dimension()) {
        p3out << "Reached the maximum number of allowed iterations!!"<<endl;
        exit(0);
      }


      //p3out << offdiagonal[iter-1]<<"  "<<iter<<endl;                                                                                                                       
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


void SpinAdapted::Linear::block_davidson(vector<Wavefunction>& b, DiagonalMatrix& h_diag, double normtol, const bool &warmUp, Davidson_functor& h_multiply, bool& useprecond, int currentRoot, std::vector<Wavefunction> &lowerStates)
{

  p3out.precision (12);
#ifndef SERIAL
  mpi::communicator world;
#endif
  int iter = 0;
  double levelshift = 0.0;
  int nroots = b.size();
  //normalise all the guess roots
  if(mpigetrank() == 0) {
    for(int i=0;i<nroots;++i)
      {
	for(int j=0;j<i;++j)
	  {
	    double overlap = DotProduct(b[j], b[i]);
	    ScaleAdd(-overlap, b[j], b[i]);
	  }
	Normalise(b[i]);
      }
  
  
    //if we are doing state specific, lowerstates has lower energy states
    if (lowerStates.size() != 0) {
      for (int i=0; i<lowerStates.size(); i++) {
	double overlap = DotProduct(b[0], lowerStates[i]);
	ScaleAdd(-overlap/DotProduct(lowerStates[i], lowerStates[i]), lowerStates[i], b[0]);
      }
      Normalise(b[0]);
    }

  }
  vector<Wavefunction> sigma;
  int converged_roots = 0;
  int maxiter = h_diag.Ncols() - lowerStates.size();
  while(true)
    {
      p3out << "\t\t\t Davidson Iteration :: " << iter << endl;

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
      //if (mpigetrank() == 0) {
      //  p3out << *bptr << endl;
      //  p3out << *sigmaptr << endl;
      //}
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
	
	for (int i = 1; i <= subspace_eigenvalues.Ncols (); ++i)
	  p3out << "\t\t\t " << i << " ::  " << subspace_eigenvalues(i,i) << endl;

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
	    p3out << "\t\t\t going back to converged root "<<i<<"  "<<rnorm<<" > "<<normtol<<endl;
            continue;
          }
        }
        r = sigma[converged_roots];
        ScaleAdd(-subspace_eigenvalues(converged_roots+1), b[converged_roots], r);

	if (lowerStates.size() != 0) {
	  for (int i=0; i<lowerStates.size(); i++) {
	    double overlap = DotProduct(r, lowerStates[i]);
	    ScaleAdd(-overlap/DotProduct(lowerStates[i], lowerStates[i]), lowerStates[i], r);
	    //ScaleAdd(-overlap, lowerStates[i], r);
	  }
	}
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


      p3out << "\t \t \t residual :: " << rnorm << endl;
      if (rnorm < normtol)
	{
	  p3out << "\t\t\t Converged root " << converged_roots << endl;

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
	      p3out << "\t\t\t Deflating block Davidson...\n";
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

	  //if we are doing state specific, lowerstates has lower energy states
	  if (lowerStates.size() != 0) {
	    for (int i=0; i<lowerStates.size(); i++) {
	      double overlap = DotProduct(r, lowerStates[i]);
	      ScaleAdd(-overlap/DotProduct(lowerStates[i], lowerStates[i]), lowerStates[i], r);
	      //ScaleAdd(-overlap, lowerStates[i], r);
	    }
	  }
	  //double tau2 = DotProduct(r,r);

	  Normalise(r);
	  b.push_back(r);

	}


    }
}

//solves the equation (H-E)^T(H-E)|psi_1> = -(H-E)^TQV|\Psi_0>  where lowerState[1] contains |\Psi_0> to enforce orthogonality (Q), and lowerState[0] contains V|\Psi_0> so we can calculate its projection in the krylov space 
//the algorithm is taken from wikipedia page
//the algorithm is taken from wikipedia page
double SpinAdapted::Linear::ConjugateGradient(Wavefunction& xi, double normtol, Davidson_functor& h_multiply, std::vector<Wavefunction>& lowerStates)
{
  setbuf(stdout, NULL);
  p3out.precision (12);
  int iter = 0, maxIter = 100;
  double levelshift = 0.0, overlap2 = 0.0, oldError=0.0, functional=0.0, Error=0.0;

  Wavefunction& targetState = lowerStates[0];
  if (mpigetrank() == 0) {
    for (int i=1; i<lowerStates.size(); i++) {
      overlap2 = pow(DotProduct(lowerStates[i], lowerStates[i]), 0.5);
      if (fabs(overlap2) > NUMERICAL_ZERO) { 
	ScaleAdd(-DotProduct(targetState, lowerStates[i])/overlap2, 
		 lowerStates[i], targetState);
	ScaleAdd(-DotProduct(xi, lowerStates[i])/overlap2, 
		 lowerStates[i], xi);
      }
    }
  }

#ifndef SERIAL
    mpi::communicator world;
    mpi::broadcast(world, xi, 0);
#endif

  Wavefunction pi, ri; 
  ri=xi; ri.Clear();
  h_multiply(xi, ri);  

  //Check if we should even perform CG or just exit with a zero vector.
  bool doCG = true;
  if (mpigetrank() == 0) {
    Wavefunction ricopy = ri; ricopy.Clear(); ricopy.Randomise();
    Wavefunction ricopy2 = ricopy;

    for (int i=1; i<lowerStates.size(); i++) {
      overlap2 = pow(DotProduct(lowerStates[i], lowerStates[i]), 0.5);
      if (fabs(overlap2) > NUMERICAL_ZERO) { 
	ScaleAdd(-DotProduct(ricopy, lowerStates[i])/overlap2, 
		 lowerStates[i], ricopy2);
      }
    }

    if (abs(DotProduct(ricopy2, targetState)) < NUMERICAL_ZERO) {
      pout << "The problem is ill posed or the initial guess is very bad "<<DotProduct(ricopy, targetState)<<endl;
      doCG = false;
    }
  }
#ifndef SERIAL
    mpi::broadcast(world, doCG, 0);
#endif
  if (!doCG) {
    xi.Clear();
    int success = 0;

    functional = 0.0;
    return functional;
  }

  if (mpigetrank() == 0) {
    ScaleAdd(-1.0, targetState, ri);
    Scale(-1.0, ri);
    
    for (int i=1; i<lowerStates.size(); i++) {
      overlap2 = pow(DotProduct(lowerStates[i], lowerStates[i]),0.5);
      if (fabs(overlap2) > NUMERICAL_ZERO) { 
	ScaleAdd(-DotProduct(ri, lowerStates[i])/overlap2, 
		 lowerStates[i], ri);
      }
    }
    
    pi = ri;

    oldError = DotProduct(ri, ri);
    printf("\t\t\t %15s  %15s  %15s\n", "iter", "Functional", "Error");
  }

#ifndef SERIAL
    mpi::broadcast(world, Error, 0);
    mpi::broadcast(world, oldError, 0);
    mpi::broadcast(world, functional, 0);
#endif

    if (oldError < normtol) {
      if (mpigetrank() == 0) {
	functional = -DotProduct(xi, ri) - DotProduct(xi, targetState);
	printf("\t\t\t %15i  %15.8e  %15.8e\n", 0, functional, oldError);
      }
#ifndef SERIAL
    mpi::broadcast(world, functional, 0);
#endif
      return functional;
    }

  while(true) {
#ifndef SERIAL
    mpi::communicator world;
    mpi::broadcast(world, pi, 0);
#endif

    Wavefunction Hp = pi; Hp.Clear();


    h_multiply(pi, Hp);

    if (mpigetrank() == 0) {

      for (int i=1; i<lowerStates.size(); i++) {
	overlap2 = pow(DotProduct(lowerStates[i], lowerStates[i]),0.5);
	if (fabs(overlap2) > NUMERICAL_ZERO) { 
	  ScaleAdd(-DotProduct(Hp, lowerStates[i])/overlap2, 
		   lowerStates[i], Hp);
	}
      }

      double alpha = oldError/DotProduct(pi, Hp);
      
      ScaleAdd(alpha, pi, xi);
      ScaleAdd(-alpha, Hp, ri);
      
      Error = DotProduct(ri, ri);
      functional = -DotProduct(xi, ri) - DotProduct(xi, targetState);
      printf("\t\t\t %15i  %15.8e  %15.8e\n", iter, functional, Error);
    }

#ifndef SERIAL
    mpi::broadcast(world, Error, 0);
    mpi::broadcast(world, functional, 0);
#endif

    if (Error < normtol || iter >maxIter) {
      return functional;
    }
    else {      
      if (mpigetrank() == 0) {
	double beta = Error/oldError;
	oldError = Error;
	ScaleAdd(1.0/beta, ri, pi);
	Scale(beta, pi);
	
	for (int i=1; i<lowerStates.size(); i++) {
	  overlap2 = pow(DotProduct(lowerStates[i], lowerStates[i]),0.5);
	  if (fabs(overlap2) > NUMERICAL_ZERO) { 
	    ScaleAdd(-DotProduct(pi, lowerStates[i])/overlap2, 
		     lowerStates[i], pi);
	  }
	}
      }
      iter ++;
    }
  }
}


void makeOrthogonalToLowerStates(Wavefunction& targetState, std::vector<Wavefunction>& lowerStates) {
  for (int i=1; i<lowerStates.size(); i++) {
    double overlap2 = pow(DotProduct(lowerStates[i], lowerStates[i]), 0.5);
    if (fabs(overlap2) > NUMERICAL_ZERO) { 
      ScaleAdd(-DotProduct(targetState, lowerStates[i])/overlap2, 
	       lowerStates[i], targetState);
    }
  }
}

double SpinAdapted::Linear::MinResMethod(Wavefunction& xi, double normtol, Davidson_functor& h_multiply, std::vector<Wavefunction>& lowerStates)
{
  setbuf(stdout, NULL);
  p3out.precision (12);
  int iter = 0, maxIter = 100;
  double levelshift = 0.0, overlap2 = 0.0, oldError=0.0, functional=0.0, Error=0.0;

  Wavefunction& targetState = lowerStates[0];
  if (mpigetrank() == 0) {
    makeOrthogonalToLowerStates(targetState, lowerStates);
    makeOrthogonalToLowerStates(xi, lowerStates);
  }

#ifndef SERIAL
    mpi::communicator world;
    mpi::broadcast(world, xi, 0);
#endif

    Wavefunction pi, ri; 
    ri=xi; ri.Clear();
    h_multiply(xi, ri);  

  //Check if we should even perform CG or just exit with a zero vector.
  bool doCG = true;
  if (mpigetrank() == 0) {
    Wavefunction ricopy = ri; ricopy.Clear(); ricopy.Randomise();
    Wavefunction ricopy2 = ricopy;

    makeOrthogonalToLowerStates(ricopy2, lowerStates);

    if (abs(DotProduct(ricopy2, targetState)) < NUMERICAL_ZERO) {
      pout << "The problem is ill posed or the initial guess is very bad "<<DotProduct(ricopy, targetState)<<endl;
      doCG = false;
    }
  }
#ifndef SERIAL
    mpi::broadcast(world, doCG, 0);
#endif
  if (!doCG) {
    xi.Clear();
    int success = 0;

    functional = 0.0;
    return functional;
  }

  if (mpigetrank() == 0) {
    ScaleAdd(-1.0, targetState, ri);
    Scale(-1.0, ri);
    
    makeOrthogonalToLowerStates(ri, lowerStates);
    
    pi = ri;

    oldError = DotProduct(ri, ri);
    printf("\t\t\t %15s  %15s  %15s\n", "iter", "Functional", "Error");
  }

#ifndef SERIAL
  mpi::broadcast(world, oldError, 0);
#endif
  
  if (oldError < normtol) {
    if (mpigetrank() == 0) {
      functional = -DotProduct(xi, ri) - DotProduct(xi, targetState);
      printf("\t\t\t %15i  %15.8e  %15.8e\n", 0, functional, oldError);
    }
#ifndef SERIAL
    mpi::broadcast(world, functional, 0);
#endif
    return functional;
  }

#ifndef SERIAL
    mpi::broadcast(world, ri, 0);
#endif
  
  double betaNumerator = 0, betaDenominator = 0;
  Wavefunction Hr = ri; Hr.Clear();
  h_multiply(ri, Hr);
  betaDenominator = DotProduct(ri, Hr);

  Wavefunction Hp = Hr;

  if (mpigetrank() == 0) {
    makeOrthogonalToLowerStates(Hp, lowerStates);
    makeOrthogonalToLowerStates(Hr, lowerStates);
  }

  while(true) {

    if (mpigetrank() == 0) {
      double alpha = DotProduct(ri, Hr)/DotProduct(Hp, Hp);
      
      ScaleAdd(alpha, pi, xi);
      ScaleAdd(-alpha, Hp, ri);
      
      Error = DotProduct(ri, ri);

      functional = -DotProduct(xi, targetState);
      printf("\t\t\t %15i  %15.8e  %15.8e \n", iter, functional, Error);
    }

#ifndef SERIAL
    mpi::broadcast(world, Error, 0);
    mpi::broadcast(world, functional, 0);
    mpi::broadcast(world, ri, 0);
#endif

    if (Error < normtol || iter >maxIter) {
      return functional;
    }
    else {      
      Hr.Clear();
      h_multiply(ri, Hr);
      if (mpigetrank() == 0) {
	makeOrthogonalToLowerStates(Hr, lowerStates);

	betaNumerator = DotProduct(ri, Hr);
	double beta = betaNumerator/betaDenominator;
	betaDenominator = betaNumerator;

	ScaleAdd(1./beta, ri, pi);
	Scale(beta, pi);

	ScaleAdd(1./beta, Hr, Hp);
	Scale(beta, Hp);
	
      }
      iter ++;
    }
  }
}

