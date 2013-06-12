#ifndef NPDM_OPERATORS_HEADER
#define NPDM_OPERATORS_HEADER
#include "BaseOperator.h"
//#include <boost/function.hpp>
//#include <boost/functional.hpp>

//typedef boost::function<void (std::vector<boost::shared_ptr<SpinAdapted::SparseMatrix> >&)> Functor;

namespace SpinAdapted{

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class DesCre: public SpinAdapted::SparseMatrix
{
  public:
    DesCre() { orbs.resize(2); fermion = false; build_pattern = "(DC)";} // default build_pattern
    void build_in_csf_space(const SpinBlock& b);
    void build(const SpinBlock& b) ;
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class CreCreDes: public SpinAdapted::SparseMatrix
{
  public:
    CreCreDes() { orbs.resize(3); fermion = true; build_pattern = "((CC)D)";} // default build_pattern
    void build_in_csf_space(const SpinBlock& b);
    void build(const SpinBlock& b) ;
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class CreDesDes: public SpinAdapted::SparseMatrix
{
  public:
    CreDesDes() { orbs.resize(3); fermion = true; build_pattern = "((CD)D)";} // default build_pattern
    void build_in_csf_space(const SpinBlock& b) {assert(false);}
    void build(const SpinBlock& b) ;
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class CreDesCre: public SpinAdapted::SparseMatrix
{
  public:
    CreDesCre() { orbs.resize(3); fermion = true; build_pattern = "((CD)C)";} // default build_pattern
    void build_in_csf_space(const SpinBlock& b);
    void build(const SpinBlock& b) ;
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class CreCreCre: public SpinAdapted::SparseMatrix
{
  public:
    CreCreCre() { orbs.resize(3); fermion = true; build_pattern = "((CC)C)";} // default build_pattern
    void build(const SpinBlock& b) ;
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class DesCreDes: public SpinAdapted::SparseMatrix
{
  public:
    DesCreDes() { orbs.resize(3); fermion = true; build_pattern = "((DC)D)";} // default build_pattern
    void build(const SpinBlock& b) ;
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class DesDesCre: public SpinAdapted::SparseMatrix
{
  public:
    DesDesCre() { orbs.resize(3); fermion = true; build_pattern = "((DD)C)";} // default build_pattern
    void build(const SpinBlock& b) ;
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

}

#endif
