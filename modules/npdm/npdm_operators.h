#ifndef NPDM_OPERATORS_HEADER
#define NPDM_OPERATORS_HEADER
#include "BaseOperator.h"
//#include <boost/function.hpp>
//#include <boost/functional.hpp>

//typedef boost::function<void (std::vector<boost::shared_ptr<SpinAdapted::SparseMatrix> >&)> Functor;

namespace SpinAdapted{

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class CreCreCre: public SpinAdapted::SparseMatrix
{
 public:
  CreCreCre() { orbs.resize(3); fermion = false;}
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

}

#endif
