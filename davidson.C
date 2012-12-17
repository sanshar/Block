/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "davidson.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "BaseOperator.h"
#include "MatrixBLAS.h"


SpinAdapted::multiply_h::multiply_h(const SpinBlock& b, const bool &onedot_) : block(b){}

void SpinAdapted::multiply_h::operator()(Wavefunction& c, Wavefunction& v)
{
  block.multiplyH( c, &v, MAX_THRD);
}


