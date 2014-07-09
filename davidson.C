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

SpinAdapted::multiply_h_e::multiply_h_e(const SpinBlock& b, const bool &onedot_, double E0_) : block(b), E0(E0_){}

void SpinAdapted::multiply_h_e::operator()(Wavefunction& c, Wavefunction& v)
{
  Wavefunction oi = v; oi.Clear();
  block.multiplyH( c, &v, MAX_THRD);
  block.multiplyOverlap(c, &oi, MAX_THRD);
  ScaleAdd(-E0, oi, v);
}




