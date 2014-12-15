/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef NPDM_EXPECTATIONS_ENGINE_H
#define NPDM_EXPECTATIONS_ENGINE_H

#include "spinblock.h"
#include "wavefunction.h"
#include "BaseOperator.h"

namespace SpinAdapted{
namespace Npdm{

  void FormLeftOp(const SpinBlock* leftBlock, const SparseMatrix& leftOp, const SparseMatrix& dotOp, SparseMatrix& Aop, int totalspin);
  double DotProduct(const Wavefunction& w1, const Wavefunction& w2, const SpinBlock& big);
  double spinExpectation(Wavefunction& wave1, Wavefunction& wave2, SparseMatrix &leftOp, SparseMatrix& dotOp, SparseMatrix& rightOp, const SpinBlock& big);
}
}

#endif

