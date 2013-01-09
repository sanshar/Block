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
					     double additional_noise)
{
  const int nroots = solution.size();

  DiagonalMatrix e;
  bool useprecond = true;
  int iter = 0;
  bool solved = false;

  while (iter <= 1 && !solved) {
    solution.resize(nroots);
    e.ReSize(big.get_stateInfo().totalStates); e= 0;
    if (dmrginp.outputlevel() > 0)
      pout << "\t\t\t Building Diagonal Hamiltonian " << endl;
    big.diagonalH ( e);
    if (dmrginp.outputlevel() > 0)
      pout << "\t\t\t Done building diagonal hamiltonian "<<endl;
    FORTINT m, n=1, nsize=e.Storage();
    pout << "\t\t\t Number of elements in wavefunction :: " << e.Ncols() << endl;
    if (mpigetrank()==0) {
      m = idamax_(nsize,e.Store(), n); 
      if (dmrginp.outputlevel() > 0)
         pout << "highest diagonal value "<<m<<" "<<e(m)<<endl;
    }


    multiply_h davidson_f(big, onedot);
    //etemp = e;
#ifndef SERIAL
    if (mpigetrank() != 0) e.ReSize(0);
#endif
    GuessWave::guess_wavefunctions(solution, e, big, guesswavetype, onedot, dot_with_sys, additional_noise); 
    Linear::block_davidson(solution, e, tol, warmUp, davidson_f, useprecond, solved);

 
    iter++;
  }
  if (iter == 2 && !solved)
  {
    cout << "Unable to solve the eigenvalue equation and bailing!!"<<endl;
    exit(2);
  }

  solution.resize(nroots);
  energies.resize(nroots);
  for (int i=0; i<nroots&& mpigetrank() == 0;i++) {
    energies[i] = e(i+1);
    //pout << "\t\t\t Energy of wavefunction "<<i<<"  =  "<<e(i+1)<<endl;
  }
#ifndef SERIAL
  mpi::communicator world;
  broadcast(world, energies, 0);
#endif
  pout<<endl;
}

