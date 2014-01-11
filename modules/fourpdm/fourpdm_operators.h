#ifndef FOURPDM_OPERATORS_HEADER
#define FOURPDM_OPERATORS_HEADER
#include "BaseOperator.h"

namespace SpinAdapted{

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
// 4PDM operators
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
// This is a skeleton class only used to set up infrastructure, and the operator contents are never stored or built during the standard sweep
class RI4index: public SpinAdapted::SparseMatrix
{
  public:
    RI4index() { orbs.resize(4); fermion = true; build_pattern = "((XX)(XX))";} // default build_pattern
//    void build_in_csf_space(const SpinBlock& b) { assert(false); };
    void build(const SpinBlock& b) { assert(false); };
    void build_from_disk(SpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs) { assert(false); };
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block) { assert(false); };
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b) { assert(false); };
};

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
