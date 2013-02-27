/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "modules/lrt/lrt_davidson.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "BaseOperator.h"
#include "MatrixBLAS.h"


SpinAdapted::LRT::multiply_h_left::multiply_h_left(const SpinBlock& b, const bool &onedot_) : block(b){}
SpinAdapted::LRT::multiply_h_total::multiply_h_total(const SpinBlock& b, const bool &onedot_) : block(b){}

void SpinAdapted::LRT::multiply_h_left::operator()(Wavefunction& c, Wavefunction& v, int iState, int jState)
{
  block.multiplyH_lrt_left( c, &v, iState, MAX_THRD);
}

void SpinAdapted::LRT::multiply_h_left::operator()(Wavefunction& c, Wavefunction& v, int iState, int jState)
{
  block.multiplyH_lrt_total( c, &v, iState, MAX_THRD);
}


