#ifndef FOURPDM_OPERATORS_HEADER
#define FOURPDM_OPERATORS_HEADER
#include "BaseOperator.h"

namespace SpinAdapted{

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
// 4PDM operators
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class DesCreDes: public SpinAdapted::SparseMatrix
{
  public:
    DesCreDes() { orbs.resize(3); fermion = true; build_pattern = "((DC)D)";} // default build_pattern
    void build(const SpinBlock& b) ;
    void build_from_disk(SpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs);
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class DesDesCre: public SpinAdapted::SparseMatrix
{
  public:
    DesDesCre() { orbs.resize(3); fermion = true; build_pattern = "((DD)C)";} // default build_pattern
    void build(const SpinBlock& b) ;
    void build_from_disk(SpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs);
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class CreCreDesDes: public SpinAdapted::SparseMatrix
{
  public:
    CreCreDesDes() { orbs.resize(4); fermion = true; build_pattern = "((CC)(DD))";} // default build_pattern
    void build(const SpinBlock& b) ;
    void build_from_disk(SpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs) { assert(false); }
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

}

#endif
