/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_LINEAR_HEADER_H
#define SPIN_LINEAR_HEADER_H
#define WANT_MATH
#define WANT_STREAM
#include <newmat.h>
#include <newmatap.h>
#include <newmatio.h>
#include <iostream>
#include <cmath>
#include <vector>
#include "wavefunction.h"
#include "davidson.h"

namespace SpinAdapted{
namespace Linear
{
  void precondition(Wavefunction& op, double e, DiagonalMatrix& diagonal, double levelshift=0.0);
  void olsenPrecondition(Wavefunction& op, Wavefunction& C0, double e, DiagonalMatrix& diagonal, double levelshift=0.0);
  void block_davidson(vector<Wavefunction>& b, DiagonalMatrix& e, double normtol, const bool &warmUp, Davidson_functor& h_mult, bool& useprecond, int currentRoot, vector<Wavefunction>& lowerStates);
  void Lanczos(vector<Wavefunction>& b, DiagonalMatrix& e, double normtol, Davidson_functor& h_multiply, int nroots);
  double ConjugateGradient(Wavefunction& xi, double normtol, Davidson_functor& h_multiply, std::vector<Wavefunction> &lowerStates);
  double MinResMethod(Wavefunction& xi, double normtol, Davidson_functor& h_multiply, std::vector<Wavefunction> &lowerStates);
};
}
#endif

