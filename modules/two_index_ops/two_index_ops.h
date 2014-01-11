#ifndef TWO_INDEX_OPERATORS_H
#define TWO_INDEX_OPERATORS_H
#include "BaseOperator.h"

namespace SpinAdapted{

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
// 3PDM operators
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class DesCre: public SpinAdapted::SparseMatrix
{
  public:
    DesCre() { orbs.resize(2); fermion = false; build_pattern = "(DC)";} // default build_pattern
    void build_in_csf_space(const SpinBlock& b);
    void build(const SpinBlock& b) ;
    void build_from_disk(SpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs) { assert(false); }
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

}

#endif
