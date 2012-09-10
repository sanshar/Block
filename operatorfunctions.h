/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is free software: you can redistribute it and/or modify         
it under the terms of the GNU General Public License as published by         
the Free Software Foundation, either version 3 of the License, or            
(at your option) any later version.                                          
                                                                             
This program is distributed in the hope that it will be useful,              
but WITHOUT ANY WARRANTY; without even the implied warranty of               
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
GNU General Public License for more details.                                 
                                                                             
You should have received a copy of the GNU General Public License            
along with this program.  If not, see <http://www.gnu.org/licenses/>.        

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

  void TensorMultiply(const SpinBlock *ablock, const Baseoperator<Matrix>& a, const SpinBlock *cblock, Wavefunction& c, Wavefunction& v, const SpinQuantum dQ, double scale, int num_thrds = 1);
  void TensorMultiply(const SpinBlock *ablock, const Baseoperator<Matrix>& a, const Baseoperator<Matrix>& b, const SpinBlock *cblock, Wavefunction& c, Wavefunction& v, const SpinQuantum opQ, double scale);
  
  void OperatorScaleAdd(double scaleV, const SpinBlock& b, const Baseoperator<Matrix>& op1, Baseoperator<Matrix>& op2);
  void MultiplyProduct(const Baseoperator<Matrix>& a, const Baseoperator<Matrix>& b, Baseoperator<Matrix>& c, Real scale);
}
}
#endif
