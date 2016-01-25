/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_OPERATORFUNCTIONS_HEADER_H
#define SPIN_OPERATORFUNCTIONS_HEADER_H
#include "Operators.h"
#include "timer.h"
#include "spinblock.h"
#include "MatrixBLAS.h"
#include <math.h>
#include "global.h"
#define TINY 1.e-20

namespace SpinAdapted{

namespace operatorfunctions
{
  
  void TensorTrace(const SpinBlock *ablock, const Baseoperator<Matrix>& a, const SpinBlock *cblock, const StateInfo *cstateinfo, Baseoperator<Matrix>& c, double scale= 1.0, int num_thrds = 1);

  void TensorTrace(const SpinBlock *ablock, const Baseoperator<Matrix>& a, const SpinBlock *cblock, const StateInfo *cstateinfo, DiagonalMatrix& c, double scale);

  void TensorTraceElement(const SpinBlock *ablock, const Baseoperator<Matrix>& a, const SpinBlock *cblock, const StateInfo *cstateinfo, Baseoperator<Matrix>& c, Matrix& cel, int cq, int cqprime);

  void TensorTraceElement(const SpinBlock *ablock, const Baseoperator<Matrix>& a, const SpinBlock *cblock, const StateInfo *cstateinfo, Baseoperator<Matrix>& c, Matrix& cel, int cq, int cqprime, double scale);


  void TensorProduct (const SpinBlock *ablock, const Baseoperator<Matrix>& a, const Baseoperator<Matrix>& b, const SpinBlock *cblock, const StateInfo *cstateinfo, Baseoperator<Matrix>& c, double scale, int num_thrds=1);

  void Product (const SpinBlock *ablock, const Baseoperator<Matrix>& a, const Baseoperator<Matrix>& b, Baseoperator<Matrix>& c, double scale);

  void TensorProduct (const SpinBlock *ablock, const Baseoperator<Matrix>& a, const Baseoperator<Matrix>& b, const SpinBlock *cblock, const StateInfo *cstateinfo, DiagonalMatrix& c, double scale);

  void TensorProductElement(const SpinBlock *ablock, const Baseoperator<Matrix>& a, const Baseoperator<Matrix>& b, const SpinBlock *cblock, const StateInfo *cstateinfo, Baseoperator<Matrix>& c, Matrix& cel, int cq, int cqprime, double scale);

  void TensorProduct (const Baseoperator<Matrix>& a, const Baseoperator<Matrix>& b, const StateInfo *brastateinfo, const StateInfo *ketstateinfo, Baseoperator<Matrix>& c, double scale, bool aIsLeftOp, int num_thrds=1);

  void TensorProductElement(const Baseoperator<Matrix>& a, const Baseoperator<Matrix>& b, const StateInfo *brastateinfo, const StateInfo *ketstateinfo, Baseoperator<Matrix>& c, Matrix& cel, int cq, int cqprime, bool aIsLeftOp, double scale);

  void TensorMultiply(const SpinBlock *ablock, const Baseoperator<Matrix>& a, const SpinBlock *cblock, Wavefunction& c, Wavefunction& v, const SpinQuantum dQ, double scale, int num_thrds = 1);
  void TensorMultiply(const SpinBlock *ablock, const Baseoperator<Matrix>& a, const Baseoperator<Matrix>& b, const SpinBlock *cblock, Wavefunction& c, Wavefunction& v, const SpinQuantum opQ, double scale);
  void TensorMultiply(const Baseoperator<Matrix>& a, const Baseoperator<Matrix>& b, const StateInfo *brastateinfo, const StateInfo *ketstateinfo, const Wavefunction& c, Wavefunction& v, const SpinQuantum opQ, bool aIsLeft, double scale);
  void TensorMultiply(const Baseoperator<Matrix>& a, const StateInfo *brastateinfo, const StateInfo *ketstateinfo, const Wavefunction& c, Wavefunction& v, const SpinQuantum dQ, bool left, double scale, int num_thrds = 1);
  
void braTensorMultiply(const SpinBlock *ablock, const Baseoperator<Matrix>& a, const SpinBlock *cblock, Wavefunction& c, Wavefunction& v, double scale, int num_thrds = 1);

  void OperatorScaleAdd(double scaleV, const SpinBlock& b, const Baseoperator<Matrix>& op1, Baseoperator<Matrix>& op2);
  void MultiplyProduct(const Baseoperator<Matrix>& a, const Baseoperator<Matrix>& b, Baseoperator<Matrix>& c, Real scale);
}
}
#endif
