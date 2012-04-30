#include "davidson.h"
#include <omp.h>
#include "BaseOperator.h"
#include "MatrixBLAS.h"


SpinAdapted::multiply_h::multiply_h(const SpinBlock& b, const bool &onedot_) : block(b){}

void SpinAdapted::multiply_h::operator()(Wavefunction& c, Wavefunction& v)
{
  block.multiplyH( c, &v, MAX_THRD);
}


