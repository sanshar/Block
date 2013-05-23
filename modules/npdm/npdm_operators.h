#ifndef NPDM_OPERATORS_HEADER
#define NPDM_OPERATORS_HEADER
#include "BaseOperator.h"
//#include <boost/function.hpp>
//#include <boost/functional.hpp>

//typedef boost::function<void (std::vector<boost::shared_ptr<SpinAdapted::SparseMatrix> >&)> Functor;

namespace SpinAdapted{

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class CreCreDes: public SpinAdapted::SparseMatrix
{
  public:
    CreCreDes() { orbs.resize(3); fermion = false;}
//MAW
    void build_in_csf_space(const SpinBlock& b) {assert(false);}
    void build(const SpinBlock& b) ;
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class CreDesDes: public SpinAdapted::SparseMatrix
{
  public:
    CreDesDes() { orbs.resize(3); fermion = false;}
//MAW
    void build_in_csf_space(const SpinBlock& b) {assert(false);}
    void build(const SpinBlock& b) ;
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

}

#endif
