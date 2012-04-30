#include "solver.h"
#include "linear.h"
#include "davidson.h"
#include "guess_wavefunction.h"

void SpinAdapted::Solver::solve_wavefunction(vector<Wavefunction>& solution, vector<double>& energies, SpinBlock& big, const double tol, 
				const guessWaveTypes& guesswavetype, const bool &onedot, const bool& dot_with_sys, const bool& warmUp,
					     DiagonalMatrix& e, double additional_noise)
{
  const int nroots = solution.size();

  //if (mpigetrank() && !BIG_DAVIDSON) // release diagonal memory    
      //e.ReSize(1); // save a spot for the lowest eigenvector only



  DiagonalMatrix etemp;

  bool useprecond = true;
  int iter = 0;
  bool solved = false;

  while (iter <= 1 && !solved) {
    solution.resize(nroots);
    e.ReSize(big.get_stateInfo().totalStates); e= 0;
    if (dmrginp.outputlevel() != 0)
      pout << "\t\t\t Building Diagonal Hamiltonian " << endl;
    big.diagonalH ( e);
    if (dmrginp.outputlevel() != 0)
      pout << "\t\t\t Done building diagonal hamiltonian "<<endl;
    int m, n=1, nsize=e.Storage();
    if (mpigetrank()==0) {
      m = idamax_(nsize,e.Store(), n); 
      if (dmrginp.outputlevel() != 0)
	pout << "highest diagonal value "<<m<<" "<<e(m)<<endl;
    }


    pout << "\t\t\t Number of elements in wavefunction :: " << e.Ncols() << endl;
    multiply_h davidson_f(big, onedot);
    etemp = e;
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
  for (int i=0; i<nroots;i++) {
    energies[i] = e(i+1);
    //pout << "\t\t\t Energy of wavefunction "<<i<<"  =  "<<e(i+1)<<endl;
  }
  pout<<endl;
  e = etemp;
}

