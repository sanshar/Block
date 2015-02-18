/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "solver.h"
#include "linear.h"
#include "davidson.h"
#include "guess_wavefunction.h"
#include "blas_calls.h"
#ifndef SERIAL
#include <boost/mpi.hpp>
#endif
#include "pario.h"


void SpinAdapted::Solver::solve_wavefunction(vector<Wavefunction>& solution, vector<double>& energies, SpinBlock& big, const double tol, 
					     const guessWaveTypes& guesswavetype, const bool &onedot, const bool& dot_with_sys, const bool& warmUp,
					     double additional_noise, int currentRoot, std::vector<Wavefunction>& lowerStates)
{
  const int nroots = solution.size();

  DiagonalMatrix e;
  bool useprecond = true;

  e.ReSize(big.get_stateInfo().totalStates); e= 0;
  p1out << "\t\t\t Building Diagonal Hamiltonian " << endl;


  big.diagonalH(e);
  p1out << "\t\t\t Done building diagonal hamiltonian "<<endl;
  FORTINT m, n=1, nsize=e.Storage();
  p2out << "\t\t\t Number of elements in wavefunction :: " << e.Ncols() << endl;
  if (mpigetrank()==0) {
    m = idamax_(nsize,e.Store(), n); 
    p3out << "\t\t\t highest diagonal value "<<m<<" "<<e(m)<<endl;
  }
  else 
    e.ReSize(0);


  bool haveEnoughStates = (e.Ncols()< nroots) ? false : true;
#ifndef SERIAL
  mpi::communicator world;
  broadcast(world, haveEnoughStates, 0);
#endif


  if (!haveEnoughStates) {
    //sometimes when you need many roots and at the start of the sweep the hilbert space is not big
    //enough to support all the roots
    
    solution.resize(nroots);

    if (dmrginp.calc_type() != RESPONSE) {
      for (int i=0; i<nroots; i++) {
	solution[i].initialise(dmrginp.effective_molecule_quantum_vec(), &big, onedot);
	solution[i].Randomise();
	Normalise(solution[i]);
      }
    }
    else {
      solution.resize(1);
      GuessWave::guess_wavefunctions(solution[0], e, big, guesswavetype, onedot, currentRoot, 
				     dot_with_sys, 0.0); 
    }
  
  }
  else {
    if(dmrginp.solve_method() == DAVIDSON) {
      solution.resize(nroots);
      multiply_h davidson_f(big, onedot);
      GuessWave::guess_wavefunctions(solution, e, big, guesswavetype, onedot, dot_with_sys, additional_noise, currentRoot); 

      if (mpigetrank() == 0) {
	for (int istate=0; istate<lowerStates.size(); istate++)  {
	  for (int jstate=istate+1; jstate<lowerStates.size(); jstate++) {
	    double overlap = DotProduct(lowerStates[istate], lowerStates[jstate]);
	    ScaleAdd(-overlap/DotProduct(lowerStates[istate], lowerStates[istate]), lowerStates[istate], lowerStates[jstate]);
	  }
	}
      }

      if (nroots == 1 && currentRoot >= e.Ncols()) //state specific calculation
	lowerStates.resize(0);

      Linear::block_davidson(solution, e, tol, warmUp, davidson_f, useprecond, currentRoot, lowerStates);

    }
    else if (dmrginp.solve_method() == CONJUGATE_GRADIENT) {
      solution.resize(1);

      //SS PLEASE REMOVE THIS
      //lowerStates.resize(2);

      multiply_h davidson_f(big, onedot);

      if (mpigetrank()!=0) 
	e.ReSize(0);

      GuessWave::guess_wavefunctions(solution[0], e, big, guesswavetype, onedot, currentRoot, 
				     dot_with_sys, 0.0); 

      if (guesswavetype == BASIC)
	solution[0].Clear();

      double functional = Linear::MinResMethod(solution[0], tol, davidson_f, lowerStates);
      if (mpigetrank() == 0)
	e(1) = functional;

    }
    else {
      solution.resize(1);
      multiply_h davidson_f(big, onedot);
      GuessWave::guess_wavefunctions(solution, e, big, guesswavetype, onedot, dot_with_sys, additional_noise, currentRoot); 
      Linear::Lanczos(solution, e, tol, davidson_f, nroots);
    }
  }

  solution.resize(nroots);
  energies.resize(nroots);
  if (haveEnoughStates) {
    for (int i=0; i<nroots&& mpigetrank() == 0;i++) {
      energies[i] = e(i+1);
      //pout << "\t\t\t Energy of wavefunction "<<i<<"  =  "<<e(i+1)<<endl;
    }
  }
  else {
    for (int i=0; i<nroots&& mpigetrank() == 0;i++) {
      if (dmrginp.calc_type() == RESPONSE)
	energies[i] = 1.e10;
      else
	energies[i] = e(1);
    }
  }
#ifndef SERIAL
  broadcast(world, energies, 0);
#endif
  pout<<endl;
}

