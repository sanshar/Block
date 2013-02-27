#ifndef SPIN_ROTATE_OPERATORS_HEADER
#define SPIN_ROTATE_OPERATORS_HEADER
 
#include "spinblock.h"
#include "Operators.h"

namespace SpinAdapted {

namespace RotateOps {

void rotate_operator_arrays(const Matrix& alpha, SpinBlock& b, const opTypes& ot, int nroots);

void rotate_operators(const Matrix& alpha, SpinBlock* b, std::vector<std::vector<boost::shared_ptr<SparseMatrix> > >& ops, int nroots);

};

};

#endif 
