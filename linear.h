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
  //void small_davidson(Wavefunction& c1, DiagonalMatrix& diagonal, double normtol, Davidson_functor& h_multiply);
  void precondition(Wavefunction& op, double e, DiagonalMatrix& diagonal, double levelshift=0.0);
  void olsenPrecondition(Wavefunction& op, Wavefunction& C0, double e, DiagonalMatrix& diagonal, double levelshift=0.0);
  void block_davidson(vector<Wavefunction>& b, DiagonalMatrix& e, double normtol, const bool &warmUp, Davidson_functor& h_mult, bool& useprecond, bool& solved);
};
}
#endif

